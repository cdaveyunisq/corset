// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <Read.h>

std::atomic<int64_t> ReadId::next_id = 0;

bool Read::has_same_alignments(const std::shared_ptr<Read> &r) {
    if (r->get_trans_hash() != get_trans_hash())
        return false;
    return (r->get_sample() == get_sample());
};

void ReadList::add_alignment(string read, string trans, int sample) {
    // find the transcript id if it already exists:
    shared_ptr<Read> r = reads_map->insert(read);
    read_ids.push_back(r->getId());

    if (!read_id_map.contains(r->getId())) {
        read_id_map.insert(std::make_pair(r->getId(), r));
    }
    shared_ptr<Transcript> t = transcript_list->insert(trans);
    // don't try to insert an alignment if 1. it already exists or
    // 2. if the transcript ID is not in the TranscriptList.
    if (!r->has(t)) {
        r->add_alignment(t);
        r->set_sample(sample);
    }
};

// This is the method called when we read equivalence class files as input.
// reads are already 'compact', so we can add to the read_vector instead of the map.
void ReadList::add_alignment(vector<string> trans_names, int sample, int weight) {
    // make the new read
    shared_ptr<Read> r = make_shared<Read>();
    reads_vector.push_back(move(r));
    r = reads_vector.back();

    r->set_sample(sample);
    r->set_weight(weight); // this read represents mutliple reads in the original bam
    if (find(read_ids.begin(), read_ids.end(), r->getId()) == read_ids.end()) {
        read_ids.push_back(r->getId());
    }
    if (!read_id_map.contains(r->getId())) {
        read_id_map.insert(std::make_pair(r->getId(), r));
    }
    // loop over all the transcripts that this read aligns to
    for (auto itrTrans = trans_names.begin(); itrTrans != trans_names.end(); itrTrans++) {
        // find the transcript object with the name
        shared_ptr<Transcript> t = transcript_list->insert(*itrTrans);
        r->add_alignment(t);
    }
};

// save memory by reducing all the reads of a sample into
// a vector of "compact reads". compact reads are stored in
// regular read object, with a weight equal to the number
// of regular reads. This saves a lot of RAM if reads from multiple
// samples are processed. The map obect (StringSet) with read IDs
// is also destroyed to save memory.
void ReadList::compactify_reads(shared_ptr<TranscriptList> trans, string outputReadsName) {

    // ── Phase 1: build idMap/nameMap, sort alignments and compute hashes ─────
    // Each read is independent — parallelise using a flat vector + atomic index.
    map<int64_t, shared_ptr<Read>> idMap{};
    map<int64_t, string>           nameMap{};

    // Flatten reads_map into a vector so we can index into it from threads.
    vector<pair<string, shared_ptr<Read>>> flat_reads(
        reads_map->begin(), reads_map->end());

    for (auto &[name, read] : flat_reads) {
        int64_t id = read->getId();
        // quicker to assign
        idMap[id] = read;
        nameMap[id] = name;
    }

    // Parallel sort + hash — each element is independent.
    const int hw       = (int)std::thread::hardware_concurrency();
    const int nthreads = hw > 1 ? hw : 1;
    cout << "compact " << outputReadsName << " with " << nthreads << " threads" << endl;
    {
        std::atomic<int> next{0};
        const int total = (int)flat_reads.size();
        auto worker = [&]() {
            int idx;
            while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total) {
                flat_reads[idx].second->sort_alignments();
                flat_reads[idx].second->set_trans_hash();
            }
        };
        vector<std::thread> threads;
        threads.reserve(nthreads);
        for (int t = 0; t < nthreads; t++) threads.emplace_back(worker);
        for (auto &t : threads) t.join();
    }

    // ── Phase 2: deduplicate reads per transcript, partitioned by hash ────────
    // Group reads by trans_hash. Two reads can only be duplicates if they share
    // the same hash — so processing each bucket in a separate thread is safe:
    // no two threads ever touch the same Read object.
    unordered_map<uintptr_t, vector<int64_t>> hash_buckets;
    for (auto &[id, read] : idMap)
        hash_buckets[read->get_trans_hash()].push_back(id);

    // Flatten buckets into a vector so threads can claim them via atomic index.
    vector<vector<int64_t>> bucket_list;
    bucket_list.reserve(hash_buckets.size());
    for (auto &[hash, ids] : hash_buckets)
        bucket_list.push_back(std::move(ids));

    // Each thread accumulates its own removableReads — merged after joining.
    vector<vector<int64_t>> per_thread_removable(nthreads);
    {
        std::atomic<int> next{0};
        const int total = (int)bucket_list.size();

        auto worker = [&](int thread_idx) {
            vector<int64_t> &my_removable = per_thread_removable[thread_idx];
            int idx;
            while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total) {
                vector<int64_t> &ids = bucket_list[idx];
                int bucket_size = (int)ids.size();
                for (int i = 0; i < bucket_size - 1; i++) {
                    int64_t readId = ids[i];
                    if (!idMap.contains(readId)) continue;
                    shared_ptr<Read> read = idMap.at(readId);
                    if (read->get_weight() == 0) continue;
                    for (int j = i + 1; j < bucket_size; j++) {
                        int64_t otherId = ids[j];
                        if (!idMap.contains(otherId)) continue;
                        if (readId == otherId) continue;
                        shared_ptr<Read> otherRead = idMap.at(otherId);
                        if (otherRead->get_weight() != 0 &&
                            read->has_same_alignments(otherRead)) {
                            read->set_weight(read->get_weight() + otherRead->get_weight());
                            otherRead->set_weight(0);
                            my_removable.push_back(otherId);
                        }
                    }
                }
            }
        };

        vector<std::thread> threads;
        threads.reserve(nthreads);
        for (int t = 0; t < nthreads; t++) threads.emplace_back(worker, t);
        for (auto &t : threads) t.join();
    }

    // Merge per-thread removable lists into one set for O(1) lookup.
    unordered_set<int64_t> removable_set;
    for (auto &rv : per_thread_removable)
        for (int64_t id : rv)
            removable_set.insert(id);

    // ── Phase 3: erase zeroed reads from read_ids (serial, cheap) ────────────
    erase_if(read_ids, [&removable_set](int64_t id) {
        return removable_set.contains(id);
    });

    // Populate reads_vector with surviving reads (weight > 0).
    reads_vector.reserve(idMap.size());
    for (auto &[id, read] : idMap) {
        if (read->get_weight() > 0)
            reads_vector.push_back(read);
    }

    reads_map->clear();
    read_id_map.clear();
}

void ReadList::write(string outputReadsName) {
    // now output the read mapping to file if requested
    ofstream readFile;
    readFile.open(outputReadsName);
    for (auto read = reads_vector.begin(); read != reads_vector.end(); read++) {
        readFile << (*read)->get_weight();
        for (auto name = (*read)->align_begin(); name != (*read)->align_end(); name++) {
            readFile << "\t" << *name;
        }
        readFile << endl;
    }
    readFile.close();
    cout << "Done writing " << outputReadsName << endl;
}


vector<int64_t> ReadList::getReadIds() {
    return read_ids;
}


map<int64_t, shared_ptr<Read>> ReadList::getReadIdMap() {
    return read_id_map;
}

shared_ptr<Read> ReadList::getRead(int64_t id) {
    return read_id_map.at(id);
}


shared_ptr<Transcript> ReadList::getTranscript(string name) {
    // contains followed by map[key] causes 2 search in the lookup.
    // use an iterator for one step.
    auto it = transcript_list->get_map().find(name);
    if (it != transcript_list->get_map().end()) {
        return it->second;
    }
    return nullptr;
}

// Re-point this ReadList to a merged TranscriptList after parallel reading.
// Must be called before any downstream use when using parallel file loading.
void ReadList::rebind_transcript_list(const shared_ptr<TranscriptList> &merged) {
    transcript_list = merged;
}

bool ReadList::serialise(const std::string &path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) {
        std::cerr << "ReadList::serialise: cannot open " << path << std::endl;
        return false;
    }

    BinaryIO::write_header(out);

    // ── TranscriptList ───────────────────────────────────────────────────────
    auto tmap = transcript_list->get_map();
    uint64_t n_trans = tmap.size();
    BinaryIO::write_pod(out, n_trans);
    for (const auto &kv : tmap)
        kv.second->serialise(out);  // delegates to Transcript::serialise

    // ── reads_vector ─────────────────────────────────────────────────────────
    uint64_t n_reads = reads_vector.size();
    BinaryIO::write_pod(out, n_reads);
    for (const auto &r : reads_vector)
        r->serialise(out);          // delegates to Read::serialise

    // ── read_ids ─────────────────────────────────────────────────────────────
    BinaryIO::write_pod_vector(out, read_ids);

    out.flush();
    return out.good();
}

std::shared_ptr<ReadList> ReadList::deserialise(
    const std::string &path,
    const std::shared_ptr<TranscriptList> &trans)
{
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) {
        std::cerr << "ReadList::deserialise: cannot open " << path << std::endl;
        return nullptr;
    }

    try {
        BinaryIO::read_header(in);
    } catch (const std::exception &e) {
        std::cerr << "ReadList::deserialise: " << e.what() << std::endl;
        return nullptr;
    }

    auto rList = std::make_shared<ReadList>(trans);

    // ── TranscriptList ───────────────────────────────────────────────────────
    uint64_t n_trans;
    BinaryIO::read_pod(in, n_trans);
    for (uint64_t t = 0; t < n_trans; ++t) {
        auto transcript = Transcript::deserialise(in); // delegates to Transcript::deserialise
        trans->insert(transcript->get_name());          // no-op if already present
    }

    // ── reads_vector ─────────────────────────────────────────────────────────
    uint64_t n_reads;
    BinaryIO::read_pod(in, n_reads);
    rList->reads_vector.reserve(n_reads);

    for (uint64_t i = 0; i < n_reads; ++i) {
        auto r = Read::deserialise(in);  // delegates to Read::deserialise

        // Rebuild Transcript::reads_ back-references for each transcript
        // this read aligns to — Transcript::deserialise does not do this
        // since it has no knowledge of which reads map to it.
        for (auto it = r->align_begin(); it != r->align_end(); ++it) {
            shared_ptr<Transcript> tran = trans->insert(*it);
            tran->add_read(r);
        }

        rList->reads_vector.push_back(r);
        rList->read_id_map.insert({r->getId(), r});  // rebuild id map
    }

    // ── read_ids ─────────────────────────────────────────────────────────────
    BinaryIO::read_pod_vector(in, rList->read_ids);

    return rList;
}

