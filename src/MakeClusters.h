// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#ifndef MAKECLUSTERS_H
#define MAKECLUSTERS_H

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <Cluster.h>
#include <Read.h>
#include <Transcript.h>

using namespace std;

// ---------------------------------------------------------------------------
// Disjoint Set Union (Union-Find) over Transcript raw pointers.
// Path compression + union by rank gives near-O(1) amortised per operation.
// ---------------------------------------------------------------------------
struct DSU {
    unordered_map<shared_ptr<Transcript>, shared_ptr<Transcript>,
                  TransPtrHash, TransPtrEqual> parent;
    unordered_map<shared_ptr<Transcript>, int,
                  TransPtrHash, TransPtrEqual> rank;

    void add(const shared_ptr<Transcript>& t) {
        if (!parent.count(t)) {
            parent[t] = t;
            rank[t]   = 0;
        }
    }

    shared_ptr<Transcript> find(const shared_ptr<Transcript>& t) {
        if (parent[t] != t)
            parent[t] = find(parent[t]);
        return parent[t];
    }

    void unite(const shared_ptr<Transcript>& a, const shared_ptr<Transcript>& b) {
        auto ra = find(a);
        auto rb = find(b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) {
            parent[ra] = rb;
        } else if (rank[ra] > rank[rb]) {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra]++;
        }
    }
};

class MakeClusters {
private:
    vector<shared_ptr<Cluster>> clusterList;
    shared_ptr<Cluster> current_cluster;

public:
    void setCurrentCluster(shared_ptr<Transcript> trans) {
        current_cluster = getMapElement(trans)->second;
    };

private:
    // O(1) lookup as opposed to map
    unordered_map<shared_ptr<Transcript>, shared_ptr<Cluster>, TransPtrHash, TransPtrEqual> transMap;

    shared_ptr<pair<shared_ptr<Transcript> const, shared_ptr<Cluster> >> getMapElement(shared_ptr<Transcript> trans);

    void checkAgainstCurrentCluster(shared_ptr<Transcript> trans);

    void makeSuperClusters(const vector<shared_ptr<ReadList>> &readLists);

    void makeSuperClustersSingleThreaded(const vector<shared_ptr<ReadList>> &readLists);

    void processSuperClusters(map<float, string> &distance_thresholds, vector<int> &groups);

public:
    MakeClusters(const vector<shared_ptr<ReadList>> &readLists, map<float, string> &distance_thresholds, vector<int> &groups);
};

#endif
