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
    if (find(read_ids.begin(), read_ids.end(), r->getId()) == read_ids.end()) {
        read_ids.push_back(r->getId());
    }
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
    reads_vector.push_back(r);
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
    // first lets sort the alignments for each read
    // and calculate a hash value to be used when comparing alignments
    // TODO: is it possible parallise this
    map<int64_t, shared_ptr<Read> > idMap{};
    // house keeping data structure to back reference keys in read_map
    map<int64_t, string> nameMap{};
    for (auto itr = reads_map->begin(); itr != reads_map->end(); itr++) {
        int64_t id = itr->second->getId();
        string name = itr->first;
        if (!idMap.contains(id)) {
            idMap.insert(std::make_pair(id, itr->second));
        }
        if (!nameMap.contains(id)) {
            nameMap.insert(std::make_pair(id, name));
        }
        itr->second->sort_alignments();
        itr->second->set_trans_hash();
    }

    // then loop over the transcripts
    vector<int64_t> removableReads{};
    for (auto transItr = trans->begin(); transItr != trans->end(); transItr++) {
        vector<int64_t> readIds = transItr->second->get_reads();
        int reads_size = readIds.size();
        for (int i = 0; i < reads_size - 1; i++) {
            int readId = readIds.at(i);
            if (idMap.find(readId) == idMap.end()) {
                // transcript references unknown read
                continue;
            }
            shared_ptr<Read> read = idMap.at(readId);
            if (read->get_weight() != 0) {
                for (int j = (i + 1); j < reads_size; j++) {
                    int otherId = readIds.at(j);
                    if (idMap.find(otherId) == idMap.end()) {
                        // unknown read
                        continue;
                    }
                    if (readId == otherId) {
                        continue;
                    }
                    shared_ptr<Read> otherRead = idMap.at(otherId);
                    if (otherRead->get_weight() != 0 &&
                        read->has_same_alignments(otherRead)) {
                        int new_weight = read->get_weight() + otherRead->get_weight();
                        read->set_weight(new_weight); // add to the weight
                        otherRead->set_weight(0); // set the weight of duplicates to zero
                        removableReads.push_back(otherId);
                    }
                }
            }
        }

        // now remove the reads with weight=0 from each transcript's list of reads
        erase_if(readIds, [&removableReads](int64_t id) {
            return (find(removableReads.begin(), removableReads.end(), id) != removableReads.end());
        });
        // remove internal references to read ids
        erase_if(read_ids, [&removableReads](int64_t id) {
            return (find(removableReads.begin(), removableReads.end(), id) != removableReads.end());
        });

    }

    // after changing the shared_ptr we can simply clear reads_map
    // instead of iterating over and delet
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
    if (transcript_list->get_map().contains(name)) {
        return transcript_list->get_map()[name];
    }
    return nullptr;
}
