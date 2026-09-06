// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/**
 * A class to store the names of transcripts.
 * This allows transcripts to be idenfitied by a pointer to their "pair" 
 * rather than using the full name (which is a string). This makes the code
 * run faster and require less memory.
 *
 * Author: Nadia Davidson
 * Modified: 3 May 3013
 */

#ifndef TRANSCRIPT_H
#define TRANSCRIPT_H

#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <StringSet.h>
#include <BinaryIO.h>

using namespace std;

class Read;

class Transcript {
    string name_;
    int pos_; // used later by Cluster
    // vector of read ids
    vector<int64_t> reads_; //temporary vector so we can quickly remove the alignments
    //for transcripts with less then min_count hits.
    //  bool reached_min_counts_;

public:
    Transcript() { name_ = ""; };

    Transcript(string name);

    string get_name() { return name_; } ;
    void pos(int position) { pos_ = position; };
    int pos() { return pos_; };

    void add_read(const std::shared_ptr<Read> &read);

    bool reached_min_counts(const std::vector<std::shared_ptr<Read> > &reads);

    //return reached_min_counts_; };

    void remove(const std::vector<std::shared_ptr<Read> > &reads); //remove myself from the reads lists ..

    static int samples;
    static int groups;
    static int min_counts;
    static int min_reads_for_link;
    static int max_alignments;

    vector<int64_t> get_reads() { return reads_; } ;

    // Binary serialisation — writes name and position only.
    // reads_ is rebuilt during Read deserialisation via add_read().
    void serialise(std::ostream &out) const {
        BinaryIO::write_string(out, name_);
        BinaryIO::write_pod(out, pos_);
    }

    // Binary deserialisation — restores name and position from stream.
    // Caller is responsible for inserting into the owning TranscriptList.
    static std::shared_ptr<Transcript> deserialise(std::istream &in) {
        std::string name;
        BinaryIO::read_string(in, name);
        int pos;
        BinaryIO::read_pod(in, pos);
        auto t = std::make_shared<Transcript>(name);
        t->pos_ = pos;
        return t;
    }
};

typedef StringSet<Transcript> TranscriptList;


// Hash and equality for shared_ptr<Transcript> keyed on the raw pointer address.
struct TransPtrHash {
    size_t operator()(const shared_ptr<Transcript>& p) const noexcept {
        return std::hash<Transcript*>()(p.get());
    }
};
struct TransPtrEqual {
    bool operator()(const shared_ptr<Transcript>& a, const shared_ptr<Transcript>& b) const noexcept {
        return a.get() == b.get();
    }
};


#endif
