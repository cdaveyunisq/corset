// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** This file contains the main function for the corset program. It
 ** controls the i/o including command line options.
 **
 ** The program runs in four main stages:
 ** 1 - All the alignments in the bam files are read. See the function
 **     read_bam_file below and add_alignment() in ReadList for more
 **     detail.
 ** 2 - Filter out any transcripts which have fewer than a certain number
 **     of reads aligning. By deafult this is 10.
 ** 3 - The read-transcript associations are parsed and the transcripts
 **     undergo an initial clustering, see MakeClusters::makeSuperClusters(),
 **     in which transcripts are grouped together if they share a single read
 **     with another transcript in the same group.
 ** 4 - For each of the clusters above, hierarchical clustering is performed,
 **     see Cluster::cluster(). The group is split into multiple smaller
 **     clusters and counts are output for each of these smaller groups.
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 25 June 2019
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <climits>
#include <map>
#include <random>
// #include <pthread.h>

#include <MakeClusters.h>
#include <Read.h>
#include <Transcript.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <thread>
#include <future>
#include <htslib/sam.h>
// #include <bam.h>

// #include <gperftools/profiler.h>

#include "version.h"

#define MAX_BAM_LINE_LENGTH 10000

using namespace std;

const string corset_extension = ".corset-reads";

// a function to parse a bam file. It required samtools
// to read it. The alignments are stored in a ReadList object.
shared_ptr<ReadList> read_bam_file(string all_file_names, const shared_ptr<TranscriptList> &trans, int sample) {
    shared_ptr<ReadList> rList = make_shared<ReadList>(trans);
    string filename;
    stringstream ss(all_file_names);

    // read each file in the comma separated list
    while (getline(ss, filename, ',')) {
        std::cout << "Reading bam file : " << filename << std::endl;
        if (!(filename.find(corset_extension) == string::npos)) {
            std::cout << "The input files look like corset read files (filenames with " << corset_extension;
            std::cout << "). Perhaps you meant to run with the -i corset option?" << std::endl;
            exit(1);
        }
        samFile *in = 0;
        if ((in = sam_open(filename.c_str(), "rb")) == 0) {
            std::cerr << "fail to open " << filename << " for reading." << std::endl;
            exit(1);
        }
        // Give HTSlib background threads for BGZF block decompression.
        // Reserve 1 core for the main reading loop; use at least 1 bg thread.
        const unsigned int hw_threads = std::thread::hardware_concurrency();
        const int bg_threads = (hw_threads > 1) ? static_cast<int>(hw_threads - 1) : 1;
        hts_set_threads(in, bg_threads);

        bam1_t *b = bam_init1();
        int r;
        int i = 0;

        // read the header first..... only really need this when we want to
        // output all the transcripts including those with no reads
        sam_hdr_t *header = sam_hdr_read(in);
        int total_targets = sam_hdr_nref(header);
        for (int tid = 0; tid < total_targets; tid++) {
            const char *target_name = sam_hdr_tid2name(header, tid);
            trans->insert(target_name);
        }

        while ((r = sam_read1(in, header, b)) >= 0) {
            // read one alignment from `in'
            // add the alignment into our data object
            string read_name = string(bam_get_qname(b));
            int tid = b->core.tid;
            if (tid != -1) {
                // unmapped reads have tid=-1
                const char *target_name = sam_hdr_tid2name(header, tid);
                if (target_name != nullptr) {
                    rList->add_alignment(read_name, string(target_name), sample);
                }
            }
            if (i % 200000 == 0) // output some info about how it's proceeding.
                std::cout << float(i) / float(1000000) << " million alignments read" << std::endl;
            //       if(i>5000000) break;
            i++;
        }
        bam_destroy1(b);
        sam_close(in);
        std::cout << "Done reading " << filename << std::endl;
    }
    rList->compactify_reads(trans); // reduce the memory usage
    return rList; // return the alignment list
}

// a function to parse a fasta file and k-mers to cluster
shared_ptr<ReadList> read_fasta_file(string all_file_names, const shared_ptr<TranscriptList> &trans, int sample) {
    shared_ptr<ReadList> rList = make_shared<ReadList>(trans);
    string filename;
    stringstream ss(all_file_names);
    ifstream file;
    string line;

    // read each file in the comma separated list
    while (getline(ss, filename, ',')) {
        std::cout << "Reading fasta file : " << filename << std::endl;

        // now open the fasta file
        file.open(filename);
        if (!(file.good())) {
            std::cout << "Unable to open file " << filename << std::endl;
            exit(1);
        }
        map<string, string> sequences;
        string id = "";
        while (getline(file, line)) {
            int start = line.find(">") + 1;
            if (start == 1) {
                // if this is the ID line...
                int end = line.find_first_of("\t\n ") - 1;
                id = line.substr(start, end);
            } else {
                sequences[id] = sequences[id] + line;
            }
        }
        std::cout << "Done reading file, getting k-mers" << std::endl;
        map<string, string>::iterator it = sequences.begin();
        map<string, string>::iterator it_end = sequences.end();
        int k = 70;
        int t_count = 1;
        for (; it != it_end; it++, t_count++) {
            if (t_count % 1000 == 0) // output some info about how it's proceeding.
                std::cout << float(t_count) / float(1000) << " thousand contigs read" << std::endl;
            string seq = it->second;
            trans->insert(it->first);
            for (int i = 0; i <= seq.length() - k; i++) {
                string kmer = seq.substr(i, k);
                rList->add_alignment(kmer, it->first, sample);
            }
        }
        std::cout << "Done reading " << filename << std::endl;
    }
    rList->compactify_reads(trans); // reduce the memory usage
    return rList; // return the alignment list
}

#define READS_REDISTRIBUTED 0
#define READS_COUNTED 1
#define READS_FILTERED 2

/** Code to add an EC into the readList. This code is common to corset and
    salmon EC file parsing **/
int add_equivalence_class(const shared_ptr<ReadList> &rList, int &sample, vector<string> &transNames, int &weight) {
    // case that the number of reads is smaller than threshold for a link b/n transcripts
    if (weight < Transcript::min_reads_for_link) {
        shuffle(transNames.begin(), transNames.end(), default_random_engine());
        for (int rr = 0; rr < transNames.size(); rr++) {
            int this_weight = weight / transNames.size();
            if (rr < (weight % transNames.size()))
                this_weight++;
            vector<string> this_tran;
            this_tran.push_back(transNames.at(rr));
            rList->add_alignment(this_tran, sample, this_weight);
        }
        return READS_REDISTRIBUTED;
    } else if (transNames.size() <= Transcript::max_alignments || Transcript::max_alignments <= 0) {
        rList->add_alignment(transNames, sample, weight);
        return READS_COUNTED;
    }
    return READS_FILTERED;
}

/** Process a corset format equivalence class file **/
shared_ptr<ReadList> read_corset_file(string all_file_names, const shared_ptr<TranscriptList> &trans, int sample) {
    shared_ptr<ReadList> rList = make_shared<ReadList>(trans);
    string filename;
    stringstream ss(all_file_names);

    // read each file in the comma separated list
    while (getline(ss, filename, ',')) {
        std::cout << "Reading corset file : " << filename << std::endl;
        // check that the filename has the correct file extension
        // i.e. that we can't accidently pass a bam file.
        if (filename.find(corset_extension) == string::npos) {
            std::cout << "The input files don't have the extension, " << corset_extension;
            std::cout << ". Please check them." << std::endl;
            exit(1);
        }
        ifstream file(filename);
        string line;
        int reads_counted = 0;
        int reads_filtered = 0;
        int reads_redistributed = 0;
        while (getline(file, line)) {
            istringstream istream(line);
            int weight;
            vector<string> transNames;
            string name;
            int alignments = -1; // start at -1 because first column is weight
            istream >> weight;
            while (istream >> name)
                transNames.push_back(name);
            int ret_flag = add_equivalence_class(rList, sample, transNames, weight);
            switch (ret_flag) {
                case READS_REDISTRIBUTED:
                    reads_redistributed += weight;
                    break;
                case READS_COUNTED:
                    reads_counted += weight;
                    break;
                case READS_FILTERED:
                    reads_filtered += weight;
                    break;
            }
        }
        std::cout << reads_counted << " reads counted, "
                << reads_filtered << " reads filtered, "
                << reads_redistributed << " reads redistributed."
                << std::endl;
    }
    return rList;
}

// Process a salom equivalence class file
// This has the same information as a corset-reads file
shared_ptr<ReadList> read_salmon_eq_classes_file(string all_file_names, const shared_ptr<TranscriptList> &trans,
                                                 int sample) {
    shared_ptr<ReadList> rList = make_shared<ReadList>(trans);
    string filename;
    stringstream ss(all_file_names);
    // read each file in the comma separated list
    while (getline(ss, filename, ',')) {
        std::cout << "Reading salmon eq_classes file : " << filename << std::endl;
        ifstream file(filename);
        string line;
        int reads_counted = 0;
        int reads_filtered = 0;
        int reads_redistributed = 0;
        // First read the header information with the transcript IDs
        getline(file, line);
        int ntranscripts = atoi(line.c_str());
        getline(file, line);
        int neqclasses = atoi(line.c_str());
        std::cout << "Reading data on " << ntranscripts << " transcripts in " << neqclasses << " equivalence classes" <<
                std::endl;
        vector<string> trans_ids;
        for (int nt = 0; nt < ntranscripts; nt++) {
            getline(file, line);
            trans_ids.push_back(line);
        }
        // Now read the equivalence class information
        for (int ne = 0; ne < neqclasses; ne++) {
            vector<string> transNames;
            getline(file, line);
            istringstream istream(line);
            int eq_size;
            istream >> eq_size;
            for (int es = 0; es < eq_size; es++) {
                int trans_pos;
                istream >> trans_pos;
                transNames.push_back(trans_ids.at(trans_pos));
            }
            int weight;
            istream >> weight;
            int ret_flag = add_equivalence_class(rList, sample, transNames, weight);
            switch (ret_flag) {
                case READS_REDISTRIBUTED:
                    reads_redistributed += weight;
                    break;
                case READS_COUNTED:
                    reads_counted += weight;
                    break;
                case READS_FILTERED:
                    reads_filtered += weight;
                    break;
            }
        }
        std::cout << reads_counted << " reads counted, "
                << reads_filtered << " reads filtered, "
                << reads_redistributed << " reads redistributed."
                << std::endl;
    }
    return rList;
}

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage() {
    std::cout << std::endl;
    std::cout << "Corset clusters contigs and counts reads from de novo assembled transcriptomes." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: corset [options] <input bam files>" << std::endl;
    std::cout << std::endl;
    std::cout << "Input bam files:" << std::endl;
    std::cout << "\t The input files should be multi-mapped bam files. They can be single, paired-end or mixed" <<
            std::endl;
    std::cout << "\t and do not need to be indexed. A space separated list should be given." << std::endl;
    std::cout << "\t e.g. corset sample1.bam sample2.bam sample3.bam" << std::endl;
    std::cout << "\t or just: corset sample*.bam" << std::endl;
    std::cout << std::endl;
    std::cout << "\t If you want to combine the results from different transcriptomes. i.e. the same reads have " <<
            std::endl;
    std::cout << "\t been mapped twice or more, you can used a comma separated list like below:" << std::endl;
    std::cout << "\t corset sample1_Trinity.bam,sample1_Oases.bam sample2_Trinity.bam,sample2_Oases.bam ..." <<
            std::endl;
    std::cout << std::endl;
    std::cout << "Options are:" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -d <double list> A comma separated list of distance thresholds. The range must be" << std::endl;
    std::cout << "\t                  between 0 and 1. e.g -d 0.4,0.5. If more than one distance threshold" <<
            std::endl;
    std::cout << "\t                  is supplied, the output filenames will be of the form:" << std::endl;
    std::cout << "\t                  counts-<threshold>.txt and clusters-<threshold>.txt " << std::endl;
    std::cout << "\t                  Default: 0.3" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -D <double>      The value used for thresholding the log likelihood ratio. The default " <<
            std::endl;
    std::cout << "\t                  value will depend on the number of degrees of freedom (which is the " <<
            std::endl;
    std::cout << "\t                  number of groups -1). By default D = 17.5 + 2.5 * ndf, which corresponds " <<
            std::endl;
    std::cout << "\t                  approximately to a p-value threshold of 10^-5, when there are fewer than" <<
            std::endl;
    std::cout << "\t                  10 groups." << std::endl;
    std::cout << std::endl;
    std::cout << "\t -I               Switches off the log likelihood ratio test and should be used " << std::endl;
    std::cout << "\t                  for downstream differential transcript usage analysis. It will prevent " <<
            std::endl;
    std::cout << "\t                  differentially expressed transcript being split into different clusters. " <<
            std::endl;
    std::cout << "\t                  This option is the equivalent of setting -D to a very large number.  " <<
            std::endl;
    std::cout << "\t                  Default: log likelihood ratio test on" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -m <int>         Filter out any transcripts with fewer than this many reads aligning." <<
            std::endl;
    std::cout << "\t                  Default: 10" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -g <list>        Specifies the grouping. i.e. which samples belong to which experimental" <<
            std::endl;
    std::cout << "\t                  groups. The parameter must be a comma separated list (no spaces), with the " <<
            std::endl;
    std::cout << "\t                  groupings given in the same order as the bam filename. For example:" << std::endl;
    std::cout << "\t                  -g Group1,Group1,Group2,Group2 etc. If this option is not used, each sample" <<
            std::endl;
    std::cout << "\t                  is treated as an independent experimental group." << std::endl;
    std::cout << std::endl;
    std::cout << "\t -p <string>      Prefix for the output filenames. The output files will be of the form" <<
            std::endl;
    std::cout << "\t                  <prefix>-counts.txt and <prefix>-clusters.txt. Default filenames are:" <<
            std::endl;
    std::cout << "\t                  counts.txt and clusters.txt" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -f <true/false>  Specifies whether the output files should be overwritten if they already exist."
            << std::endl;
    std::cout << "\t                  Default: false" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -n <string list> Specifies the sample names to be used in the header of the output count file." <<
            std::endl;
    std::cout << "\t                  This should be a comma separated list without spaces." << std::endl;
    std::cout << "\t                  e.g. -n Group1-ReplicateA,Group1-ReplicateB,Group2-ReplicateA etc." << std::endl;
    std::cout << "\t                  Default: the input filenames will be used." << std::endl;
    std::cout << std::endl;
    std::cout << "\t -r <true/true-stop/false> " << std::endl;
    std::cout << "\t                  Output a file summarising the read alignments. This may be used if you" <<
            std::endl;
    std::cout << "\t                  would like to read the bam files and run the clustering in seperate runs" <<
            std::endl;
    std::cout << "\t                  of corset. e.g. to read input bam files in parallel. The output will be the" <<
            std::endl;
    std::cout << "\t                  bam filename appended with " << corset_extension << "." << std::endl;
    std::cout << "\t                  Default: false" << std::endl;
    std::cout << std::endl;
    std::cout << "\t -i <bam/corset/salmon_eq_classes>  The input file type. Use -i corset, if you previously ran" <<
            std::endl;
    std::cout << "\t                  corset with the -r option and would like to restart using those" << std::endl;
    std::cout << "\t                  read summary files. Use salmon_eq_classes, if you aligned with salmon with" <<
            std::endl;
    std::cout << "\t                  the flag --dumpEq and are passing corset the equivalent class files." <<
            std::endl;
    std::cout << "\t                  Running with either -i corset or salmon_eq_classes will switch off the -r option."
            << std::endl;
    std::cout << "\t                  Default: bam" << std::endl;
    std::cout << std::endl;
    std::cout <<
            "\t -l <int>         If running with -i corset or salmon_eq_classes, this will filter out a link between contigs"
            << std::endl;
    std::cout <<
            "\t                  if the link is supported by less than this many reads (performed sample-wise). Reads will "
            << std::endl;
    std::cout << "\t                  be reassigned uniformly to the contigs in the equivalence class. This option will"
            << std::endl;
    std::cout <<
            "\t                  improve runtime and memory usage, but will increase the number of clusters reported."
            << std::endl;
    std::cout << "\t                  Default: 1 (no filtering)" << std::endl;
    std::cout << std::endl;
    std::cout <<
            "\t -x <int>         If running with -i corset or salmon_eq_classes, this option will filter out reads that"
            << std::endl;
    std::cout << "\t                  align to more than x contigs. Default: no filtering" << std::endl;
    std::cout << std::endl;
    std::cout << "Citation: Nadia M. Davidson and Alicia Oshlack, Corset: enabling differential gene expression " <<
            std::endl;
    std::cout << "          analysis for de novo assembled transcriptomes, Genome Biology 2014, 15:410" << std::endl;
    std::cout << std::endl;
    std::cout << "Please see https://github.com/Oshlack/Corset/wiki for more information" << std::endl;
    std::cout << std::endl;
}

// the real stuff starts here.
int main(int argc, char **argv) {
    // default parameters:
    string distance_string("0.3");
    bool force = false;
    string sample_names;
    int c;
    int params = 1;
    vector<int> groups;
    bool output_reads = false;
    bool stop_after_read = false;
    string input_type = "bam";

    // function pointer to the method to read the bam or corset input files
    shared_ptr<ReadList> (*read_input)(string, const shared_ptr<TranscriptList> &, int) = read_bam_file;

    std::cout << std::endl;
    std::cout << "Running Corset Version " << CORSET_VERSION_STRING << std::endl;

    // parse the command line options
    while ((c = getopt(argc, argv, "f:p:d:n:g:D:Im:r:i:l:x:")) != EOF) {
        switch (c) {
            case 'f': {
                // f=force output to be overwritten?
                std::string value(optarg);
                transform(value.begin(), value.end(), value.begin(), ::tolower);
                if (value.compare("true") == 0 | value.compare("t") == 0 | value.compare("1") == 0) {
                    force = true;
                    std::cout << "Setting output files to be overridden" << std::endl;
                } else if (value.compare("false") == 0 | value.compare("f") == 0 | value.compare("0") == 0)
                    force = false;
                else {
                    cerr << "Unknown argument passed with -f. Please specify true or false." << endl;
                    print_usage();
                    exit(1);
                }
                params += 2;
                break;
            }
            case 'p': // p=the prefix for the output files
                std::cout << "Setting output filename prefix to " << optarg << std::endl;
                Cluster::file_prefix = string(optarg) + "-";
                params += 2;
                break;
            case 'd': // d=the distance threshold for the clustering
                std::cout << "Setting distance threshold to: " << optarg << std::endl;
                distance_string = optarg;
                params += 2;
                break;
            case 'n': {
                // names to give the samples in the output files
                std::cout << "Setting sample names to:" << optarg << std::endl;
                sample_names = string(optarg);
                replace(sample_names.begin(), sample_names.end(), ',', '\t');
                params += 2;
                break;
            }
            case 'g': {
                // the sample groups. Used for the likelihood ratio test
                std::cout << "Setting sample groups:" << optarg;
                stringstream ss(optarg);
                string s;
                vector<string> names;
                while (getline(ss, s, ','))
                    names.push_back(s);
                int ngroups = 0;
                for (int s1 = 0; s1 < names.size(); s1++) {
                    // find the first occurance
                    int s2 = 0;
                    for (; s2 < names.size(); s2++) {
                        if (names.at(s1).compare(names.at(s2)) == 0)
                            break;
                    }
                    if (groups.size() <= s2) {
                        groups.push_back(ngroups);
                        ngroups++;
                    } else
                        groups.push_back(groups.at(s2));
                }
                Transcript::groups = ngroups;
                std::cout << ", " << Transcript::groups << " groups in total" << std::endl;
                params += 2;
                break;
            }
            case 'D': {
                std::cout << "Setting likelihood threshold at " << optarg << std::endl;
                Cluster::D_cut = atof(optarg);
                params += 2;
                break;
            }
            case 'I': {
                std::cout << "Switching likelihood test off  - will not separate differentially expressed isoforms" <<
                        std::endl;
                Cluster::D_cut = INT_MAX;
                params += 1;
                break;
            }
            case 'm': {
                std::cout << "Setting minimum counts to " << optarg << std::endl;
                Transcript::min_counts = atoi(optarg);
                params += 2;
                break;
            }
            case 'r': {
                // r=whether to output read information
                std::string value(optarg);
                transform(value.begin(), value.end(), value.begin(), ::tolower);
                if (value.compare("true-stop") == 0) {
                    output_reads = true;
                    stop_after_read = true;
                    std::cout << "Setting read alignments to be output to file." << std::endl;
                    std::cout << "Corset will exit after reading bam files." << std::endl;
                } else if (value.compare("true") == 0 | value.compare("t") == 0) {
                    output_reads = true;
                    std::cout << "Setting read alignments to be output to file." << std::endl;
                } else if (value.compare("false") != 0 & value.compare("f") != 0) {
                    cerr << "Unknown argument passed with -r. Please specify true or false." << endl;
                    print_usage();
                    exit(1);
                }
                params += 2;
                break;
            }
            case 'i': {
                std::string value(optarg);
                transform(value.begin(), value.end(), value.begin(), ::tolower);
                if (value.compare("corset") == 0) {
                    read_input = read_corset_file;
                    output_reads = false; // override the -r option
                    stop_after_read = false;
                } else if (value.compare("salmon_eq_classes") == 0) {
                    read_input = read_salmon_eq_classes_file;
                    output_reads = false;
                    stop_after_read = false;
                } else if (value.compare("fasta") == 0) {
                    read_input = read_fasta_file;
                } else if (value.compare("bam") != 0 && value.compare("fasta") != 0) {
                    cerr << "Unknown input type, " << value << ", passed with -i. Please check options." << endl;
                    print_usage();
                    exit(1);
                }
                params += 2;
                break;
            }
            case 'l': {
                std::cout << "Setting minimum reads for a link to " << optarg
                        << " (only used if -i is set to corset or salmon_eq_classes)" << std::endl;
                Transcript::min_reads_for_link = atoi(optarg);
                params += 2;
                break;
            }
            case 'x': {
                std::cout << "Setting maximum alignments for a read to " << optarg
                        << " (only used if -i is set to corset or salmon_eq_classes )" << std::endl;
                Transcript::max_alignments = atoi(optarg);
                params += 2;
                break;
            }
            case '?': // other unknown options
                std::cerr << "Unknown option.. stopping" << std::endl;
                print_usage();
                exit(1);
                break;
        }
    }

    // convert the distances provided as a string into doubles
    map<float, string> distance_thresholds;
    char *cstr = new char[distance_string.length() + 1];
    strcpy(cstr, distance_string.c_str());
    char *single_string = strtok(cstr, ",");
    while (single_string != NULL) {
        double value = atof(single_string);
        // check that the value is valid
        if (value < 0 | value > 1) {
            cerr << "The distance provided is invalid. It must be between 0 and 1 (inclusive)." << endl;
            print_usage();
            exit(1);
        }
        distance_thresholds[value] = string("-") + single_string;
        single_string = strtok(NULL, ",");
    }
    if (distance_thresholds.size() == 1)
        distance_thresholds.begin()->second = ""; // no suffix if only one distance cut-off

    // the remaining command line arguments are names of the input files
    int smpls = argc - params;
    if (smpls == 0) {
        cerr << "No input files specified" << endl;
        print_usage();
        exit(1);
    }
    if (sample_names == "") {
        // if the -n option (sample names) have not been
        // specified, then use the filenames.
        for (int s = params; s < argc; s++) {
            sample_names += string(argv[s]);
            if (s < (argc - 1))
                sample_names += '\t';
        }
    }
    // if the number of sample provided for the header of the counts file
    // is inconsistent with the number of samples print an error and exit
    int n_sample_names = count(sample_names.begin(), sample_names.end(), '\t') + 1;
    if (n_sample_names != smpls) {
        cerr << "The number of sample names passed (via the -n option), " << n_sample_names
                << ", does not match the number of samples, " << smpls
                << ". Please check how many values you have passed." << endl;
        exit(1);
    }

    if (groups.size() == 0) {
        // if no -g option (no groups information) is given
        // treat each file as a separate group
        for (int s = 0; s < smpls; s++)
            groups.push_back(s);
        Transcript::groups = smpls;
    }

    // before we go any further. Check that there are enough groups and
    // sample names for all the files
    if (groups.size() != smpls) {
        cerr << "The number of experimental groups passed (via the -g option) "
                << "does not match the number of input files. Please check how "
                << "many values you have passed." << endl;
        exit(1);
    }

    // if D was not set, then set the default value
    if (Cluster::D_cut == 0)
        Cluster::D_cut = 17.5 + 2.5 * (Transcript::groups - 1);

    // test whether the output files already exist
    for (map<float, string>::iterator it = distance_thresholds.begin();
         stop_after_read == false & it != distance_thresholds.end(); ++it) {
        string type[] = {Cluster::file_counts, Cluster::file_clusters};
        for (int t = 0; t < 2; t++) {
            string filename(Cluster::file_prefix + type[t] + it->second + Cluster::file_ext);
            ifstream file(filename.c_str());
            if (file.good()) {
                if (force & remove(filename.c_str()) != 0) {
                    std::cerr << "Could not replace the file, " << filename << std::endl;
                    std::exit(1);
                } else if (!force) {
                    std::cerr << "File already exists, " << filename << ". Use \"-f true\" option to overwrite " <<
                            std::endl;
                    std::exit(1);
                }
            }
            file.close();
            ofstream ofile(filename.c_str());
            // output the header for the counts files.
            if (type[t].compare(Cluster::file_counts) == 0) {
                ofile << '\t' + sample_names << endl;
            }
            ofile.close();
        }
    }

    // set the number of samples
    Transcript::samples = smpls;
    // set-up out data structures ready to receive the alignment information
    shared_ptr<TranscriptList> tList = make_shared<TranscriptList>();
    vector<shared_ptr<ReadList> > rList{}; // one ReadList per sample
    // a map to access relationship between read id and read list.
    map<int64_t, shared_ptr<ReadList> > read_list_map{};

    // ---------------------------------------------------------------------------
    // (parallel): read each input file into its own private TranscriptList
    // and ReadList. No shared state between threads — no locks needed.
    // ---------------------------------------------------------------------------
    using ReadResult = pair<shared_ptr<TranscriptList>, shared_ptr<ReadList> >;
    vector<future<ReadResult> > futures;
    futures.reserve(smpls);

    for (int bam_file = 0; bam_file < smpls; bam_file++) {
        string file_arg = string(argv[params + bam_file]);
        futures.push_back(std::async(std::launch::async, [=]() -> ReadResult {
            // Each thread owns its own TranscriptList — no sharing, no locks.
            shared_ptr<TranscriptList> privateTrans = make_shared<TranscriptList>();
            shared_ptr<ReadList> readList = read_input(file_arg, privateTrans, bam_file);
            return {privateTrans, readList};
        }));
    }

    // ---------------------------------------------------------------------------
    // Phase 2 (serial): collect futures, merge all private TranscriptLists into
    // tList, rebind each ReadList to the unified tList, then build read_list_map.
    // ---------------------------------------------------------------------------
    rList.resize(smpls);
    vector<shared_ptr<TranscriptList> > privateTransLists(smpls);

    for (int bam_file = 0; bam_file < smpls; bam_file++) {
        auto [privateTrans, readList] = futures[bam_file].get();
        privateTransLists[bam_file] = privateTrans;
        rList[bam_file] = readList;
    }

    // Merge transcript names from all private lists into the single shared tList.
    for (auto &privateTrans: privateTransLists) {
        for (auto it = privateTrans->begin(); it != privateTrans->end(); ++it) {
            tList->insert(it->first); // no-op if name already present
        }
    }
    // Rebind each ReadList to the unified tList and populate read_list_map.
    for (int bam_file = 0; bam_file < smpls; bam_file++) {
        rList[bam_file]->rebind_transcript_list(tList);

        for (auto idItr = rList[bam_file]->getReadIds().begin();
             idItr != rList[bam_file]->getReadIds().end(); ++idItr) {
            if (!read_list_map.contains(*idItr)) {
                read_list_map.insert(make_pair(*idItr, rList[bam_file]));
            }
        }

        // This is where we output the read alignment summary file for future runs of corset
        if (output_reads) {
            string bam_filename = string(argv[params + bam_file]);
            string base_filename = bam_filename.substr(bam_filename.find_last_of("/\\,") + 1);
            rList[bam_file]->write(base_filename + corset_extension);
        }
    }

    std::cout << "Done reading all files. " << std::endl;

    if (stop_after_read == true)
        std::exit(0);

    // output clusters with no reads in this special case
    int n = 0;
    for (auto it = tList->begin(); it != tList->end(); it++) {
        // extract the reads for this transcript.
        vector<int64_t> read_ids = it->second->get_reads();
        vector<shared_ptr<Read> > read_ref{};
        for (auto idItr = read_ids.begin(); idItr != read_ids.end(); idItr++) {
            if (read_list_map.contains(*idItr)) {
                read_ref.push_back(read_list_map.at(*idItr)->getRead(*idItr));
            }
        }

        if (!(*it).second->reached_min_counts(read_ref)) {
            if (Transcript::min_counts == 0) {
                n++;
                // output clusters
                for (map<float, string>::iterator dis = distance_thresholds.begin(); dis != distance_thresholds.end();
                     ++dis) {
                    ofstream clusterFile;
                    string filename(Cluster::file_prefix + Cluster::file_clusters + dis->second + Cluster::file_ext);
                    clusterFile.open(filename.c_str(), ios_base::app);
                    clusterFile << (*it).second->get_name() << "\t" << Cluster::cluster_id_prefix_no_reads << n << endl;
                    clusterFile.close();
                }
            }
            it->second->remove(read_ref);
        }
    };


    std::cout << "Start to cluster the reads" << std::endl;

    // Now the reads are parsed and the clustering and counting is performed.
    //   ProfilerStart("gprof.prof");
    // TODO: The make clusters constructor triggers the clustering on the rList.
    // Need to investigate how to parallelise rList based on parallel reads.
    // process may need to incorporate reads from a subset and clustering on that subset.
    // then have a separate step to combine clusters.
    MakeClusters cList(rList, distance_thresholds, groups);
    //  ProfilerStop();

    std::cout << "Finished" << std::endl;

    return (0);
}
