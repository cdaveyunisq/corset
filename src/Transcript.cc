// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

//#include <string>
//#include <vector>
//#include <sstream>
//#include <cstdlib>
//#include <algorithm>
//#include <StringSet.h>
#include <Transcript.h>
#include <Read.h>

using namespace std;

// initialise the number of samples to 0.
// this will be set later on.
int Transcript::samples = 0;
int Transcript::groups = 0;
int Transcript::min_counts = 10;
int Transcript::min_reads_for_link = 1;
int Transcript::max_alignments = -1;

void Transcript::remove(const std::vector<std::shared_ptr<Read> > &reads) {
    for (int r = 0; r < reads.size(); r++) {
        if (reads.at(r)->has(this->get_name())) {
            reads.at(r)->remove(this->get_name());
        }
    }
};

Transcript::Transcript(string name) {
    name_ = name;
    //    reached_min_counts_=false;
};

void Transcript::add_read(const std::shared_ptr<Read> &read) {
    //  if(reads_.size()>=min_counts)//{
    //    reads_.clear();
    //   reached_min_counts_=true;
    //} else {
    reads_.push_back(read->getId());
    //  }
};

bool Transcript::reached_min_counts(const std::vector<std::shared_ptr<Read> > &reads) {
    int counts = 0;
    for (int i = 0; i < reads.size(); i++) {
        if (!reads.at(i)->has(this->get_name())) {
            continue;
        }
        counts += reads.at(i)->get_weight();
        if (counts >= min_counts) return true;
    }
    return false;
}
