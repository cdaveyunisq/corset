// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <MakeClusters.h>
using namespace std;


shared_ptr<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> > > MakeClusters::getMapElement(
    shared_ptr<Transcript> trans) {
    //  StringSet<Cluster>::iterator it=transMap.find(trans);
    auto it = transMap.find(trans);
    if (it != transMap.end())
        return make_shared<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> > >(*it);
    shared_ptr<Cluster> clust = make_shared<Cluster>();
    clust->add_tran(trans);
    clusterList.push_back(clust);

    auto newTrans = transMap.insert(
        pair<shared_ptr<Transcript> const, shared_ptr<Cluster> >(trans, clust));
    return make_shared<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> > >(*newTrans.first);
};

void MakeClusters::checkAgainstCurrentCluster(shared_ptr<Transcript> trans) {
    //compare the current transcript's cluster to the "current" cluster.
    shared_ptr<Cluster> this_cluster = getMapElement(trans)->second;
    //"this_cluster" should match "current_cluster"
    if (this_cluster == current_cluster) return;

    //otherwise, change the cluster of all transcripts in "this_cluster" to "current_cluster"
    for (int i = 0; i < this_cluster->n_trans(); i++) {
        transMap[this_cluster->get_tran(i)] = current_cluster;
        current_cluster->add_tran(this_cluster->get_tran(i));
    }
    //and copy all the reads across too
    for (int i = 0; i < this_cluster->n_reads(); i++) {
        current_cluster->add_read(this_cluster->get_read(i));
    }
    //clean up
    clusterList.erase(find(clusterList.begin(), clusterList.end(), this_cluster));
};

// this method uses the unordered map data structures to improve
// iteration slightly
void MakeClusters::makeSuperClustersSingleThreaded(const vector<shared_ptr<ReadList>> &readLists) {
    //loop through each read:
    int i = 0;
    unordered_map<string, shared_ptr<Transcript>> transCache {};
    
    for (int sample = 0; sample < readLists.size(); sample++) {
        const shared_ptr<ReadList>& reads = readLists.at(sample);
        const vector<shared_ptr<Read>> readsVec = reads->getReads();
        for (auto rIt = readsVec.begin(); rIt != readsVec.end(); rIt++) {
            const shared_ptr<Read>& r = *rIt;
            // continue of the Read is already processed
            if (!r->try_mark_processed()) {
                continue;
            }
            //     if(r->get_weight() >= 0 && distance(tIt,tIt_end) <= 50){
            for (auto tIt = r->align_begin(); tIt != r->align_end(); tIt++) {
                const string& name = *tIt;

                if (!transCache.count(name)) {
                    const shared_ptr<Transcript> trans = reads->getTranscript(name);
                    if (trans == nullptr) {
                        continue;
                    }
                    transCache[name] = trans;
                }
                const auto& item = transCache.find(name);
                if (item == transCache.end()) {
                    continue;
                }
                const shared_ptr<Transcript>& trans = item->second;

                if (tIt == r->align_begin()) {
                    setCurrentCluster(trans);
                    current_cluster->add_read(r);
                } else
                    checkAgainstCurrentCluster(trans);
            }
            //	}
            if (i % 100000 == 0)
                cout << float(i) / float(1000000) << " million equivalence classes read" << endl;
            i++;
        }
    }
    //clean up
    transMap.clear();
}

void MakeClusters::makeSuperClusters(const vector<shared_ptr<ReadList> > &readLists) {
    // ── Flatten all reads across all samples into a single indexed vector ────
    // Each entry records the read and which ReadList it came from (for
    // transcript lookup via the per-sample transCache).
    struct ReadEntry {
        const shared_ptr<Read>& read;
        const unordered_map<string, shared_ptr<Transcript>>& cache {};
    };

    struct Edge {
        shared_ptr<Transcript> a;
        shared_ptr<Transcript> b;
    };
    struct ReadAssignment {
        const shared_ptr<Read>& read;
        shared_ptr<Transcript> root; // first transcript of this read
    };

    // Build a per-sample transcript name → pointer cache once, serially.
    // This is cheap (one pass over the transcript map per sample).
    vector<unordered_map<string, shared_ptr<Transcript> > > perSampleCache(readLists.size());
    for (int s = 0; s < (int) readLists.size(); s++) {
        const auto &tmap = readLists[s]->get_transcript_map();
        for (const auto &[name, ptr]: tmap)
            perSampleCache[s][name] = ptr;
    }

    vector<ReadEntry> allReads;
    allReads.reserve(1 << 20); // rough pre-allocation

    for (int s = 0; s < (int) readLists.size(); s++) {
        const auto &rv = readLists[s]->getReads();
        for (const auto &r: rv)
            allReads.push_back({r, perSampleCache[s]});
    }
    const int totalReads = (int) allReads.size();

    // ── Phase 1 (parallel): collect edges ───────────────────────────────────
    // Each thread processes a contiguous slice of allReads and records
    // (Transcript*, Transcript*) pairs for every pair of transcripts that
    // co-occur on a single read.  No shared state is written.
    //
    // We also record every read->cluster assignment: the first transcript of
    // each read nominates the "root" transcript, and the read belongs to
    // whichever cluster that root ends up in after DSU.

    const int hw = (int) std::thread::hardware_concurrency();
    const int nthreads = hw > 1 ? hw : 1;

    vector<vector<Edge> > perThreadEdges(nthreads);
    vector<vector<ReadAssignment> > perThreadAssignments(nthreads);
    vector<unordered_map<shared_ptr<Transcript>, int, TransPtrHash, TransPtrEqual>> perThreadTranscriptCounts(nthreads);

    std::atomic<int> next{0};

    auto worker = [&](int tid) {
        auto& myEdges       = perThreadEdges[tid];
        auto& myAssignments = perThreadAssignments[tid];
        auto& myCounts      = perThreadTranscriptCounts[tid]; 

        int idx;
        while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < totalReads) {
            const auto& entry = allReads[idx];
            const shared_ptr<Read> r     = entry.read;
            const auto& cache = entry.cache;

            shared_ptr<Transcript> root;
            auto & alignments = r->getAlignments();

            // set the reader as processed and ignore it if it was already processed
            // only one thread will successfully process this read.
            if (!r->try_mark_processed()) {
                continue;
            }
            
            for (auto tIt = alignments.begin(); tIt != alignments.end(); ++tIt) {
                auto cIt = cache.find(*tIt);
                if (cIt == cache.end()) continue;
                // deliberately need to copy the shared pointer at this stage
                // otherwise the root object ends up empty in later iter
                shared_ptr<Transcript> t = cIt->second;
                // record the transcript count for this thread
                myCounts[t]++;
                // make it more obvious that the very first alignment determines the root transcript
                if (tIt == alignments.begin()) {
                    root = t;
                    myAssignments.push_back({entry.read, root});
                } else {
                    if (root != nullptr) {
                        myEdges.push_back({root, t});
                    }
                }
                
            }
        }
    };

    vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; t++)
        threads.emplace_back(worker, t);
    for (auto& t : threads) t.join();

    // combine transcript counts from all threads
    unordered_map<shared_ptr<Transcript>, int, TransPtrHash, TransPtrEqual> globalTranscriptCounts;
    for (const auto& threadCounts : perThreadTranscriptCounts) {
        for (const auto& [t, count] : threadCounts) {
            globalTranscriptCounts[t] += count;
        }
    }


    // ── Phase 2 (serial): build DSU and apply all edges ─────────────────────
    //CustomDSU dsu;
    
    DSU dsu;

    for (const auto& cache : perSampleCache) {
        for (const auto& [name, ptr] : cache) {   
            // check that we have at least counted the transcript.
            // some transcripts may exist with 0 reads we skip those.
            auto countIt = globalTranscriptCounts.find(ptr);
            if (countIt != globalTranscriptCounts.end()) {
                dsu.add(ptr);
            }
        }
    }

    // for any connected component - an edge between node a and node b
    // unite the nodes in the DSU by finding the node a, and assigning it as the parent of node b. 
    for (auto& threadEdges : perThreadEdges)
        for (auto& e : threadEdges)
            dsu.unite(e.a, e.b);

    // ── Phase 3 (serial): build clusterList from DSU components ─────────────
    unordered_map<shared_ptr<Transcript>, shared_ptr<Transcript>,
                  TransPtrHash, TransPtrEqual> resolvedRoot;
    // assign each parent child pair to the resolved root collection.
    for (auto& [t, p] : dsu.parent)
        resolvedRoot[t] = dsu.find(t);

    unordered_map<shared_ptr<Transcript>, shared_ptr<Cluster>,
                  TransPtrHash, TransPtrEqual> rootToCluster;
    // add each transcript to its corresponding cluster based on the resolved root.
    for (auto& [t, root] : resolvedRoot) {
        if (!rootToCluster.count(root)) {
            rootToCluster[root] = make_shared<Cluster>();
            clusterList.push_back(rootToCluster[root]);
        }
        rootToCluster[root]->add_tran(t);
    }

    // Assign reads to their cluster via the recorded root transcript.
    for (auto& threadAssignments : perThreadAssignments) {
        for (auto& a : threadAssignments) {
            const shared_ptr<Transcript>& root = resolvedRoot.at(a.root);
            auto cIt = rootToCluster.find(root);
            if (cIt != rootToCluster.end())
                cIt->second->add_read(a.read);
        }
    }

    
}

void MakeClusters::processSuperClusters(map<float, string> &distance_thresholds, vector<int> &groups) {
    //now do the hierarchical clustering...
    cout << "Starting hierarchial clustering..." << endl;
    for (int i = 0; clusterList.size() > 0; i++) {
        if (i % 1000 == 0)
            cout << float(i) / float(1000) << " thousand clusters done" << endl;
        shared_ptr<Cluster> back = clusterList.back();
        back->set_id(i);
        back->set_sample_groups(groups);
        back->cluster(distance_thresholds);
        //    back->print_alignments();
        clusterList.pop_back();
    }
}

MakeClusters::MakeClusters(const vector<shared_ptr<ReadList> > &readLists,
                           map<float, string> &distance_thresholds,
                           vector<int> &groups) {
    //stage 1: process all the reads and looked for shared hits.
    //groups all transcripts which share at least one read
#ifdef CLUSTER_SYNC 
    makeSuperClustersSingleThreaded(readLists);
#else
    makeSuperClusters(readLists);
#endif
    cout << clusterList.size() << " super clusters formed." << endl;
    
    //stage 2: loop over each of the newly created groups (super clusters)
    //and perform the hierarchical clustering
    processSuperClusters(distance_thresholds, groups);
};
