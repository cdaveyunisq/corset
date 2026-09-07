Corset is a command-line program to go from a de novo transcriptome assembly to gene-level counts. Our software takes a set of reads that have been multi-mapped to the transcriptome (where multiple alignments per read were reported) and hierarchically clusters the transcripts based on the proportion of shared reads and expression patterns. It will report the clusters and gene-level counts for each sample, which are easily tested for differential expression with count based tools such as edgeR and DESeq.

See our [wiki](https://github.com/Oshlack/Corset/wiki) for downloads and instructions.


## Changes

This branch is intended to modify the corset library to:

- use cmake to build the corset library and htslib dependency. Update corset to use htslib 1.24
- incrementally build clusters to enable processes to be stopped and continue from where they previously left off.

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
    - This version applies a parallel cluster initialisation method which differs from the original, it uses a "distributed set union" algorithm to form the initial clusters. I am not sure if its results are consistent with the original algorithm or not and it needs testing.
    - It parallelises the distance calculation in each individual cluster for a small speedup.

- build/corset_sync
    - This version uses the synchronous cluster initialisation, but changes the data structures for faster data structure internally for O(1) access in reads.



Both versions read input files in parallel and will write the processed input file after compaction to binary recovery files for fast startup if the process is interrupted if the recovery switch "-R true" is specified.
 
For example the following includes a -R true switch, after reading the data it will write a file into its current working directory with the extension "".corset-recovery" if the process has read and compacted the entire file, and the "-R true" switch is set, it will look for the files having the same original name with the additional extension "corset-recovery" from its current working directory which in the example below is ~/test_corset_par. This means the process can be restarted and if the corset-recovery files already exist it will read the transcript and reads from the binary recovery file, instead of processing the entire bam file again. Reading from the recovery file is very fast (the compaction method was the bottleneck in this case). 
 
I am still debugging the corset_par version in comparison to corset_sync, the corset_par version I think is not quite right, but if I can get the DSU algorithm right, it will be a magnitude faster than corset_sync. However, corset_sync should still be a little faster than the default corset, and it has the ability to recover after a restart and load data more quickly if it has already read and compacted the BAM input files.
 
Example usage below:


```
#!/bin/bash
mkdir -p examples/Cx_sitiens_par
cp -f *.corset-recovery examples/Cx_sitiens_par
cp -f build/corset_par examples/Cx_sitiens_par
pushd examples/Cx_sitiens_par
echo "Running corset_par in examples/Cx_sitiens_par $(pwd)"
./corset_par -f true -R true -p Cx_sitiens_par ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam
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
./corset_sync -f true -R true -p Cx_sitiens_sync ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam
popd
```


# Testing Parallel execution difference

The resources used on our cluster by the parallel processing method for processing 2 BAM files of about 4Gb each with more than 100 million transcripts:

Note that data was loaded from the binary recovery files.

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

TODO: get the historical walltime as this is will execute overnight.
```
Job Id: 875887.hpc-clm-prd-t1
    Job_Name = test_clust_sync.pbs

```

