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

htslib is automatically cloned and compiled statically during the build process and statically linked into the corset binary.

Building is achieved with the build script:

```

rmdir -rf build

mkdir build

cd build
cmake ..
cmake --build .
```

Resulting binaries are located in the build directory:

- build/corset
- build/corset_fasta_ID_changer