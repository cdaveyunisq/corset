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


class CustomDSU {
public:
    std::unordered_map<std::shared_ptr<Transcript>, std::shared_ptr<Transcript>,
                       TransPtrHash, TransPtrEqual> parent;

    void add(const std::shared_ptr<Transcript>& t) {
        if (parent.find(t) == parent.end()) {
            parent[t] = t;
        }
    }

    std::shared_ptr<Transcript> find(const std::shared_ptr<Transcript>& t) {
        if (parent[t] == t) {
            return t;
        }
        // Path compression is still safe to use here because it only flattens the resolved paths
        return parent[t] = find(parent[t]);
    }

    // DIRECTIONAL UNION: Force 'b' to merge into 'a's tree, matching the original logic
    void unite(const std::shared_ptr<Transcript>& a, const std::shared_ptr<Transcript>& b) {
        auto root_a = find(a);
        auto root_b = find(b);

        if (root_a != root_b) {
            // In the original code, 'this_cluster' (b) is merged INTO 'current_cluster' (a)
            // Therefore, the root of b must now point directly to the root of a.
            parent[root_b] = root_a;
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
