# Input Files

Input files are formats containing sequence alignment maps.

The corset program supports several formats:

- Binary Alignment Map (BAM)
    - https://samtools.github.io/hts-specs/SAMv1.pdf
    - https://en.wikipedia.org/wiki/BAM_(file_format)

- Salmon Equivalence Classes
    - https://salmon.readthedocs.io/en/latest/file_formats.html

- FASTA file format
    - https://en.wikipedia.org/wiki/FASTA_format

# How inputs are normalised.

Inputs are read in the main corset.cc file.

The read signature is common between the different read methods:

```
ReadList *read_bam_file(string all_file_names, TranscriptList *trans, int sample);
```
```
ReadList *read_corset_file(string all_file_names, TranscriptList *trans, int sample)
```
```
ReadList *read_salmon_eq_classes_file(string all_file_names, TranscriptList *trans, int sample);
```

