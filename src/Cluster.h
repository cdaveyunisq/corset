// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 ** A Cluster object is responsible for holding the information 
 ** associated with a "superCluster" group of transcripts. It also
 ** contains the methods to cluster a "superCluster" into the final
 ** smaller cluster groups. ie. it contains the code for hierarcical
 ** clustering and expression testing. Cluster objects are first created
 ** in the MakeClusters class, then the Cluster::cluster method must be
 ** called to actually perform the clustering. See the description
 ** below for an outline of the Cluster::cluster algorithm.
 **
 ** Last Modified: 11th July 2014, NMD
 **/

#ifndef ONECLUSTER_H
#define ONECLUSTER_H

//#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <climits>
#include <Read.h>
#include <cstdint>
#include <mutex>

using namespace std;

typedef vector<vector<int> > group;
typedef vector<vector<vector<int> > > read_group;

#ifndef UNORDEREDMAP
typedef map<uint64_t, unsigned char> dist_map;
typedef map<uint64_t, unsigned char>::iterator dist_iterator;
#else
typedef unordered_map<uint64_t, unsigned char> dist_map;
typedef unordered_map<uint64_t, unsigned char>::iterator dist_iterator;
#endif

/** Sparse distance matrix backed by two parallel vectors:
 **   indices_ : sorted flat keys (ntrans*i + j) of non-zero cells
 **   values_  : distance values where values_[k] corresponds to indices_[k]
 **
 ** Lookup: convert (i,j) to flat key, binary-search indices_.
 **   - Not found → distance is 0 (no link)
 **   - Found at position k → values_[k] is the distance
 **
 ** Parallel write safety: indices_ is pre-allocated from the known set of
 ** non-zero pairs before any parallel work begins. Each thread owns a unique
 ** slot index and writes only to values_[k] for its assigned k — no two
 ** threads ever touch the same element.
 **/
class DistanceMatrix {
private:
    vector<uint64_t> indices_; // sorted flat keys of non-zero cells
    vector<unsigned char> values_; // values_[k] <-> indices_[k]
    uint64_t ntrans_;

public:
    void set_size(int ntrans) {
        ntrans_ = ntrans;
    }

    // Pre-allocate both vectors from the known non-zero flat keys.
    // Keys must be sorted and deduplicated before calling this.
    // Called once serially before any parallel distance computation.
    void allocate(const vector<uint64_t> &sorted_keys) {
        indices_ = sorted_keys;
        values_.assign(sorted_keys.size(), 0);
    }

    // Return the slot index for a flat key, or -1 if not present.
    // O(log P) via binary search on the sorted indices_ vector.
    int find_slot(uint64_t key) const {
        auto it = std::lower_bound(indices_.begin(), indices_.end(), key);
        if (it == indices_.end() || *it != key)
            return -1;
        return (int) std::distance(indices_.begin(), it);
    }

    unsigned char get(int i, int j) const {
        int slot = find_slot(ntrans_ * i + j);
        return (slot == -1) ? 0 : values_[slot];
    }

    // Safe to call concurrently for different (i,j) pairs:
    // each maps to a unique pre-allocated slot in values_.
    void set(int i, int j, unsigned char value) {
        int slot = find_slot(ntrans_ * i + j);
        if (slot != -1)
            values_[slot] = value;
    }

    // Direct slot write — used during parallel init when the caller
    // already holds the slot index, avoiding a second binary search.
    void set_by_slot(int slot, unsigned char value) {
        values_[slot] = value;
    }

    bool no_link(int i, int j) const {
        int slot = find_slot(ntrans_ * i + j);
        return (slot == -1) || (values_[slot] == 0);
    }


    int get_i(uint64_t key) const { return (int) (key / ntrans_); }
    int get_j(uint64_t key) const { return (int) (key % ntrans_); }


    // Iteration over non-zero entries — used by find_next_pair.
    // Walks values_ directly; skips zeros.
    struct Entry {
        uint64_t first;
        unsigned char second;
    };

    struct Iterator {
        const DistanceMatrix *mat;
        int slot;

        void advance() {
            while (slot < (int) mat->values_.size() && mat->values_[slot] == 0)
                ++slot;
        }

        Iterator &operator++() {
            ++slot;
            advance();
            return *this;
        }

        // Post-increment: return copy of current state, then advance.
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator!=(const Iterator &o) const { return slot != o.slot; }
        Entry operator*() const { return {mat->indices_[slot], mat->values_[slot]}; }

        // Allows ditr->first and ditr->second to work in find_next_pair
        // by returning a temporary Entry by value via a proxy pointer.
        struct Proxy {
            Entry e;
            const Entry *operator->() const { return &e; }
        };

        Proxy operator->() const { return Proxy{{mat->indices_[slot], mat->values_[slot]}}; }
    };


    Iterator begin() {
        Iterator it{this, 0};
        it.advance();
        return it;
    }

    Iterator end() { return Iterator{this, (int) values_.size()}; }

    int size() const { return (int) indices_.size(); }

    // Zeroes the value at slot (i,j). The slot remains in indices_ but
    // values_[slot]==0 is treated as "no link" by get(), no_link() and
    // the iterator — matching the original sparse map sentinel convention.
    // Called during merge() to invalidate rows/columns of absorbed transcripts.
    void remove(int i, int j) {
        int slot = find_slot(ntrans_ * i + j);
        if (slot != -1)
            values_[slot] = 0;
    }
};

class Cluster {
private:
    //object which store the data
    vector<shared_ptr<Transcript> > cluster_;
    // a quick lookup of transcripts by their name.
    map<string, shared_ptr<Transcript> > transcript_map;
    vector<shared_ptr<Read> > read_;
    //3D vector. Dimensions are: sample, transcript, read
    read_group read_groups;
    group groups;
    //a 2D array to hold the number of reads in read_group taking into account the weights
    group read_group_sizes;

    //represent the distances using an unsigned char array with range 0-255
    //(0 = distance of 0, 255 = distance of 1). Precision doesn't matter
    //since the dist array is updated using the reads and not from the
    //current dist array. This implementation saves lots of memory when the
    //array is large.
    DistanceMatrix dist;

    int id_;
    vector<int> sample_groups;

    // a flag for whether is the minimum distance found still zero?
    //bool zero_dist_done;
    // variable below holds the corrdinates in dist for the last minimum found
    // (only used if the last found was a zero).
    //int zero_dist_i;
    //int zero_dist_j;

    //private methods called from "cluster"

    /** get_dis is responsible for returning the distance between
     ** two clusters (the clusters have positions in the distance
     ** matrix of i and j). This first calculates the distance based on the
     ** proportion of shared reads between the clusters. It then performs
     ** a log likelihood ratio test looking for a difference in the ratio
     ** of expressions cluster i / cluster j between the different experimental
     ** groups. If this is found, then the distance is increased to the maximum.
     ** The return values will be a distance between 0 and 1. NOTE that the
     ** value returned is actually 1-distance described in our paper. So
     ** 1 = clusters share all reads, 0 = clusters share no reads.
     **/
    float get_dist(int i, int j);

    /** find_next_pair will search through the distance matrix and find
     ** the next smalled distance value. i.e. the next two pairs of
     ** clusters to be merged together (returned as max_i and max_j)
     **/
    unsigned char find_next_pair(int &max_i, int &max_j);

    /** the merge method takes two clusters (at position i and j in
     ** the distance matrix) and will merge them together. This involves
     ** merging their list of reads together, recalculating distance
     ** based on the new allocation of reads and updating the distance matrix.
     **/
    void merge(int i, int j);

    /** Before starting to cluster all transcripts are allocated to
     ** either own cluster, and the distance matrix is set-up */
    void initialise_matrix();

    /** get_counts will work out how many reads have been allocated to
     ** each cluster given the current configuration of clusters. */
    vector<int> get_counts(int s);

    /** output_cluster will write to the output "clusters.txt" and
     ** "counts.txt" files with the current configuration of clusters. */
    void output_clusters(string threshold);

public: //basic getter/setter functions
    void add_tran(shared_ptr<Transcript> trans);

    void add_read(shared_ptr<Read> read);

    shared_ptr<Transcript> get_tran(int i);

    // lookup a transcript by na
    shared_ptr<Transcript> get_tran(string name);

    shared_ptr<Read> get_read(int i);

    int n_trans();

    int n_reads();

    int get_id() { return id_; };
    void set_id(int id) { id_ = id; };
    void set_sample_groups(vector<int> &sgroups) { sample_groups = sgroups; };

    /** This method does all the work! A basic outline is:
     ** - initialise_matrix() is called to set up the distance matrix
     ** - Until we have merged all the transcripts into a single cluster:
     **     - Find the closest clusters to merge with find_next_pair()
     **     - Are we over the distance threshold now? If yes then
     **       call get_counts() followed by output_clusters().
     **     - merge the two closest cluster together and update the
     **       distance matrix with merge(). merge() will call get_dist()
     **       to recalculate the distances.
     ** - When finished, free the memory we've been using with clean_up()
     **/
    void cluster(map<float, string> &thresholds);

    //     const static unsigned char distMAX;
    static float D_cut;
    static string file_prefix;
    const static string file_counts;
    const static string file_clusters;
    const static string file_ext;
    const static string cluster_id_prefix;
    const static string cluster_id_joiner;
    const static string cluster_id_prefix_no_reads;

    // Protects concurrent appends to the shared output files during
    // parallel processSuperClusters.
    static std::mutex output_mutex;

    void print_alignments();
};

#endif
