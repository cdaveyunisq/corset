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
#include <future>
using namespace std;


shared_ptr<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> > > MakeClusters::getMapElement(
    shared_ptr<Transcript> trans) {
    //  StringSet<Cluster>::iterator it=transMap.find(trans);
    map<shared_ptr<Transcript>, shared_ptr<Cluster> >::iterator it = transMap.find(trans);
    if (it != transMap.end())
        return make_shared<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> > >(*it);
    shared_ptr<Cluster> clust = make_shared<Cluster>();
    clust->add_tran(trans);
    clusterList.push_back(clust);

    pair<map<shared_ptr<Transcript> const, shared_ptr<Cluster> >::iterator, bool> newTrans = transMap.insert(
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
    for (int sample = 0; sample < readLists.size(); sample++) {
        ReadList reads = *(readLists.at(sample));

        for (auto rIt = reads.begin(); rIt != reads.end(); rIt++) {
            shared_ptr<Read> r = *rIt;
            //     if(r->get_weight() >= 0 && distance(tIt,tIt_end) <= 50){
            for (auto tIt = r->align_begin(); tIt != r->align_end(); tIt++) {
                string name = *tIt;
                shared_ptr<Transcript> trans = reads.getTranscript(name);
                if (trans == nullptr) {
                    continue;
                }
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

    // Assign stable IDs before going parallel so each cluster has a unique id_
    for (int i = 0; i < (int)clusterList.size(); i++) {
        clusterList[i]->set_id(i);
        clusterList[i]->set_sample_groups(groups);
    }
    // Launch one async task per cluster — each cluster() call is self-contained.
    // output_clusters() inside each call is guarded by Cluster::output_mutex.
    vector<future<void>> futures;
    futures.reserve(clusterList.size());
    for (auto &clust : clusterList) {
        futures.push_back(std::async(std::launch::async, [&clust, &distance_thresholds]() {
            clust->cluster(distance_thresholds);
        }));
    }

    for (auto &f : futures)
        f.get(); // propagate any exceptions

    clusterList.clear();

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
