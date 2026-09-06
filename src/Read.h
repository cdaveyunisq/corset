// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/***
 ** This file contains the declarations and definitions for two classes:
 ** Read and ReadList.
 **
 ** Author: Nadia Davidson
 ** Modified: 11th July 2014
 **/

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <Transcript.h>
#include <memory>
#include <cstdint>
#include <atomic>
#include <limits>
#include <BinaryIO.h>
#include <unordered_map>
#include <unordered_set>
#include <thread>

using namespace std;

class ReadId {
protected:
    static constexpr int64_t MAX_ID = std::numeric_limits<int64_t>::max();

    static std::atomic<int64_t> next_id;
    friend class Read;
};

// Read is a basic container object for a read. It can record the read ID,
// all the transcript alignments for that read and which sample the read belongs to.
class Read {
public:
    Read() : id(ReadId::next_id) {
        if (ReadId::next_id >= ReadId::MAX_ID) {
            throw std::runtime_error("Maximum number of reads reached");
        }
        ++ReadId::next_id;
        weight_ = 1;
    };

    Read(string name) : id(ReadId::next_id) {
        if (ReadId::next_id >= ReadId::MAX_ID) {
            throw std::runtime_error("Maximum number of reads reached");
        }
        ++ReadId::next_id;
        weight_ = 1;
    };

    Read(const Read &other) : id(other.id),
                              weight_(other.weight_),
                              sample_(other.sample_),
                              transcriptIds_(other.transcriptIds_) {
        if (id >= ReadId::MAX_ID) {
            throw std::runtime_error("Maximum read ID limit reached during copy.");
        }
        set_trans_hash();
    }


    // dummy for StringSet
    void set_sample(int sample) { sample_ = sample; };
    int get_sample() { return sample_; };
    int alignments() { return transcriptIds_.size(); };
    void set_weight(int weight) { weight_ = weight; };
    int get_weight() { return weight_; };

    void add_alignment(const std::shared_ptr<Transcript> &trans) {
        transcriptIds_.push_back(trans->get_name());
        trans->add_read(make_shared<Read>(*this));
    };

    const vector<string>& getAlignments() const { return transcriptIds_; }

    vector<std::string>::iterator align_begin() { return transcriptIds_.begin(); };
    vector<std::string>::iterator align_end()  { return transcriptIds_.end(); };

    void sort_alignments() { sort(transcriptIds_.begin(), transcriptIds_.end()); };

    // function to check whether an alignment has already been recorded.
    bool has(const std::shared_ptr<Transcript> &t) {
        return (find(align_begin(), align_end(), t->get_name()) != align_end());
    };

    bool has(const std::string &name) {
        return (find(align_begin(), align_end(), name) != align_end());
    }

    void remove(const std::shared_ptr<Transcript> &trans) {
        transcriptIds_.erase(find(align_begin(), align_end(), trans->get_name()));
    };


    void remove(const std::string &name) {
        transcriptIds_.erase(find(align_begin(), align_end(), name));
    };


    // check if two reads align to all the same transcripts
    // this assumed that "sort_alignments" has been called first.
    // it is used in ReadList::compactify_reads
    bool has_same_alignments(const std::shared_ptr<Read> &r);

    // return the hashed values of all transcript pointers that this read aligns to
    // used to quickly compare alignments in the compactify_reads function
    uintptr_t get_trans_hash() {
        return trans_hash;
    }

    // sets the hashed values of all transcript pointers that this read aligns to using xor
    // used to quickly compare alignments in the compactify_reads function
    uintptr_t set_trans_hash() {
        trans_hash = 0;
        if (transcriptIds_.size() > 0) {
            for (auto t_itr = transcriptIds_.begin(); t_itr != transcriptIds_.end(); t_itr++) {
                // this has function was derived from boost::hash_combine()
                std::size_t name_hash = std::hash<std::string>{}(*t_itr);
                trans_hash ^= hash<uintptr_t>()(name_hash + 0x9e3779b9 + (trans_hash << 6) + (
                                                    trans_hash >> 2));
            }
        }
        return trans_hash;
    }


    void print_alignments() {
        if (transcriptIds_.size() > 0) {
            for (auto t_itr = transcriptIds_.begin(); t_itr != transcriptIds_.end(); t_itr++) {
                cout << *t_itr << " ";
            }
        }
        cout << endl;
    }

    int64_t getId() {
        return id;
    }

    Read &operator=(const Read &other) {
        if (this != &other) {
            id = other.id;
            transcriptIds_ = other.transcriptIds_;
            weight_ = other.weight_;
            sample_ = other.sample_;
            this->set_trans_hash();
        }
        return *this;
    }

    // Binary serialisation — writes all fields needed for reconstruction.
    // trans_hash is not written — it is recomputed on deserialise().
    void serialise(std::ostream &out) const {
        BinaryIO::write_pod(out, id);
        int32_t w = weight_;
        int32_t s = static_cast<int32_t>(sample_);
        BinaryIO::write_pod(out, w);
        BinaryIO::write_pod(out, s);
        BinaryIO::write_string_vector(out, transcriptIds_);
    }

    // Binary deserialisation — restores a Read from a stream.
    // The original id is restored directly without incrementing next_id.
    // trans_hash is recomputed after all transcriptIds_ are loaded.
    // Caller (ReadList::deserialise) must call tran->add_read() for each
    // transcript name to rebuild Transcript::reads_ back-references.
    static std::shared_ptr<Read> deserialise(std::istream &in) {
        int64_t id;
        BinaryIO::read_pod(in, id);
        int32_t weight, sample;
        BinaryIO::read_pod(in, weight);
        BinaryIO::read_pod(in, sample);
        std::vector<std::string> tids;
        BinaryIO::read_string_vector(in, tids);

        // Use the private constructor path — create a default Read then
        // overwrite its fields. The default constructor increments next_id
        // which we correct below.
        auto r = std::make_shared<Read>();
        r->id = id;
        r->weight_ = weight;
        r->sample_ = static_cast<unsigned char>(sample);
        r->transcriptIds_ = std::move(tids);
        r->set_trans_hash();

        // Keep next_id ahead of any restored id to avoid future collisions.
        int64_t expected = r->id + 1;
        int64_t current = ReadId::next_id.load(std::memory_order_relaxed);
        while (current < expected &&
               !ReadId::next_id.compare_exchange_weak(
                   current, expected, std::memory_order_relaxed)) {
        }

        return r;
    }

private:
    vector<std::string> transcriptIds_;
    unsigned char sample_;
    int weight_;
    uintptr_t trans_hash;

    int64_t id;
};

// ReadList is a container for a set of Reads and contains functions for
// inserting and accessing reads.
class ReadList {
private:
    std::shared_ptr<StringSet<Transcript> > transcript_list{};
    std::shared_ptr<StringSet<Read> > reads_map{}; // warning this is deleted after the reads are read
    // the StringSet above may be replaced by the internal map for
    // the id value to the read this way a read can be accessed more quickly
    // than iterating over the string set.
    map<int64_t, shared_ptr<Read> > read_id_map{};

    vector<std::shared_ptr<Read> > reads_vector{};
    // housekeeping for read ids referened by this list.
    vector<int64_t> read_ids{};

public:
    // we need to know all the transcripts before we can build a ReadList.
    ReadList(const std::shared_ptr<TranscriptList> &transcripts) {
        transcript_list = transcripts;
        reads_map = std::make_shared<StringSet<Read> >();
        reads_vector = std::vector<std::shared_ptr<Read> >{};
        read_ids = vector<int64_t>{};
        read_id_map = map<int64_t, shared_ptr<Read> >{};
    };

    // add a new alignment into the list
    // this one is used when reading bam files
    void add_alignment(string read, string trans, int sample);

    // this one is used when reading corset read summary files
    void add_alignment(vector<string> trans, int sample, int weight);

    // saves memory by reducing each read into a set of
    //"compact reads" with a weight.
    // Also the map object is clear and the reads are
    // stored as a vector instead. Read IDs are cleared.
    void compactify_reads(shared_ptr<TranscriptList> trans, string outputReadsName = "");

    const vector<std::shared_ptr<Read> >& getReads() const {
        return reads_vector;
    }

    // Returns a const reference to the underlying transcript name->ptr map.
    // Used by MakeClusters DSU phase to build the per-sample transcript cache.
    const auto& get_transcript_map() const {
        return transcript_list->get_map();
    }


    vector<std::shared_ptr<Read> >::iterator begin() { return reads_vector.begin(); };
    vector<std::shared_ptr<Read> >::iterator end() { return reads_vector.end(); };

    void write(string outputReadsName);

    vector<int64_t> getReadIds();

    map<int64_t, shared_ptr<Read> > getReadIdMap();

    shared_ptr<Read> getRead(int64_t id);

    shared_ptr<Transcript> getTranscript(const string& name) const;

    // Re-point this ReadList to a merged TranscriptList after parallel reading.
    // Must be called before any downstream use when using parallel file loading.
    void rebind_transcript_list(const shared_ptr<TranscriptList> &merged);

    // Binary serialisation — writes the TranscriptList, all compactified Reads,
    // and read_ids to a single binary file. Must be called after compactify_reads().
    // Returns false if the file could not be opened or writing failed.
    bool serialise(const std::string &path) const;

    // Binary deserialisation — restores a ReadList and populates the provided
    // TranscriptList. read_id_map and trans_hash are rebuilt automatically.
    // Transcript::reads_ back-references are also rebuilt.
    // Returns nullptr if the file cannot be read or the header is invalid.
    static std::shared_ptr<ReadList> deserialise(
        const std::string &path,
        const std::shared_ptr<TranscriptList> &trans);
};

#endif
