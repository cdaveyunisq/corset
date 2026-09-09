Corset is a command-line program to go from a de novo transcriptome assembly to gene-level counts. This software takes a set of reads that have been multi-mapped to the transcriptome (where multiple alignments per read were reported) and hierarchically clusters the transcripts based on the proportion of shared reads and expression patterns. It will report the clusters and gene-level counts for each sample, which are easily tested for differential expression with count based tools such as edgeR and DESeq.

See the original project wiki at: [wiki](https://github.com/Oshlack/Corset/wiki) for downloads and instructions.


## Changes

This branch is intended to modify the corset library to:

- use cmake to build the corset library and htslib and deflate dependencies. Update corset to use htslib 1.24
- parameterise htslib with number of available hardware cores to allow htslib to use threading when reading BAM files
- read input files in parallel and run compaction in a worker pool. 
- save recovery files for each input file after compaction completes, this allows very fast startup if the process is interrupted and needs to be restarted.
- initialise the clusters in parallel
- perform distance calculations in parallel - only effective for large cluster populations.


After the changes a dramatic speedup is evident, especially in clustering.
Loading from the saved recovery files also makes the ability to repeat the process by incrementally adding other files if required for different experiments.


## Why incrementally build clusters?

Running a process on the HPC on a large dataset can take many hours in excess of multiple weeks if processing 15 billion sequences.
However, maintenance windows or other events will require the process to be terminated.

The results end up being discarded, and the job itself cannot be restarted from the point where it left off.

The aim of this change is to allow incremental processing so that the job may be terminated and restarted from a known point.

# Building

The toolchains for cmake as well as the automake and autoconf dependencies are required (htslib requires the autoconf and automake dependencies).

Libraries needed are:

- CURL::libcurl
- OpenSSL::Crypto 
- ZLIB::ZLIB
- BZip2::BZip2
- LibLZMA::LibLZMA
- pthread

libdeflate and htslib are automatically cloned and compiled statically during the build process and statically linked into the corset binary.

Building is achieved with the build script:

```

rmdir -rf build

mkdir build

cd build
cmake ..
cmake --build .
```

Resulting binaries are located in the build directory:
- build/corset_fasta_ID_changer
    - the original version of the fasta ID changer.

- build/corset_par
    - Parallel reading of input files with serialised binary written after compaction for fast recovery.
    - This version applies a parallel cluster initialisation method which differs from the original, it uses a parallel [disjoint set union](https://en.wikipedia.org/wiki/Disjoint-set_data_structure) algorithm to form the initial clusters.
    - It parallelises the distance calculation within clusters after initialisation for a small speedup (only useful if the cluster population is large).

- build/corset_sync
    - Parallel reading of input files with serialised binary written after compaction for fast recovery.
    - This version uses the synchronous cluster initialisation, but changes the data structures for faster data structure internally for O(1) access in reads.
    - It parallelises the distance calculation within clusters after initialisation for a small speedup (only useful if the cluster population is large).

Both versions read input files in parallel and will write the processed input file after compaction to binary recovery files for fast startup if the process is interrupted if the recovery switch "-R true" is specified.
 
For example the following includes a -R true switch, after reading the data it will write a file into its current working directory with the extension "".corset-recovery" if the process has read and compacted the entire file, and the "-R true" switch is set, it will look for the files having the same original name with the additional extension "corset-recovery" from its current working directory which in the example below is ~/test_corset_par. This means the process can be restarted and if the corset-recovery files already exist it will read the transcript and reads from the binary recovery file, instead of processing the entire bam file again. Reading from the recovery file is very fast (the compaction method was the bottleneck in this case). 
 
I am still debugging the corset_par version in comparison to corset_sync, the corset_par version it produces cluster distributions similar to the corset_sync although clusters are in a different order because of the parallel Disjoint Set Union applied during cluster initialisation.

However, corset_sync should still be a little faster than the default corset as it has parallelisation in the distance calculations within the clusters, but the cluster initialisation remains the bottleneck as it uses the same initialisation procedure. 

Both have the ability to recover after a restart and load data more quickly after having saved the compacted data structures after reading the input files.
 
# Example usage:

Note on our HPC we use a modules environment, this program depends on libbz2 so it needs to be loaded into the environment using the command:

```
module load bzip2/1.0.8-gcc-p2k
```


```
#!/bin/bash
mkdir -p examples/Cx_sitiens_par
cp -f *.corset-recovery examples/Cx_sitiens_par
cp -f build/corset_par examples/Cx_sitiens_par
pushd examples/Cx_sitiens_par
echo "Running corset_par in examples/Cx_sitiens_par $(pwd)"
./corset_par -f true -R true -p Cx_sitiens_par ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam &>> corset_par.log
popd
 
```

and for the synchronous version:

```
#!/bin/bash
mkdir -p examples/Cx_sitiens_sync
cp -f *.corset-recovery examples/Cx_sitiens_sync
cp -f build/corset_sync examples/Cx_sitiens_sync
pushd examples/Cx_sitiens_sync
echo "Running corset_sync in examples/Cx_sitiens_sync $(pwd)"
./corset_sync -f true -R true -p Cx_sitiens_sync ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam &>> corset_sync.log
popd
```

An example pbs job script can be built as follows:

```
#!/bin/bash
#
# Run corset_par
#

#PBS -P MyProjectName
#PBS -l ncpus=128
#PBS -l mem=900gb
#PBS -l walltime=72:00:00
#PBS -l host=targetHost

cd ~/corset/
./test_clust_par.sh

```

The host parameter is optional, but having the ability to see which nodes are busy and which are not, (```pbsnodes -aSj```) its possible 
to target a node with plenty of free resources.

# Testing Parallel Speedup

The resources used on the UniSQ HPC cluster by the parallel processing method for processing 2 BAM files of about 4Gb each with more than 100 million transcripts:

Note that data was loaded from the binary recovery files therefore isolating the comparison to the differences in cluster initialisation.

```
Job Id: 875885
    Job_Name = test_clust_par.pbs
    resources_used.cpupercent = 387
    resources_used.cput = 00:20:18
    resources_used.mem = 3193816kb
    resources_used.ncpus = 128
    resources_used.vmem = 11213100kb
    resources_used.walltime = 00:06:41
```

The resources used for the same data for the synchronous initialisation method was:

```
Job Id: 875887
    Job_Name = test_clust_sync.pbs
    resources_used.cpupercent = 297
    resources_used.cput = 03:18:18
    resources_used.mem = 2694036kb
    resources_used.ncpus = 128
    resources_used.vmem = 11074436kb
    resources_used.walltime = 03:09:55

```

