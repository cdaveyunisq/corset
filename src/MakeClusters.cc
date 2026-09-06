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

void MakeClusters::makeSuperClusters(const vector<shared_ptr<ReadList>> &readLists) {
    //loop through each read:
    int i = 0;
    unordered_map<string, shared_ptr<Transcript>> transCache {};

    for (int sample = 0; sample < readLists.size(); sample++) {
        const shared_ptr<ReadList>& reads = readLists.at(sample);
        const vector<shared_ptr<Read>> readsVec = reads->getReads();
        for (auto rIt = readsVec.begin(); rIt != readsVec.end(); rIt++) {
            const shared_ptr<Read>& r = *rIt;
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

MakeClusters::MakeClusters(const vector<shared_ptr<ReadList>> &readLists,
                           map<float, string> &distance_thresholds,
                           vector<int> &groups) {
    //stage 1: process all the reads and looked for shared hits.
    //groups all transcripts which share at least one read
    makeSuperClusters(readLists);
    //stage 2: loop over each of the newly created groups (super clusters)
    //and perform the hierarchical clustering
    processSuperClusters(distance_thresholds, groups);
};
