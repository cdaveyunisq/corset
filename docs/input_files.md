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

# Original example program:

```
corset -f true -p Cx_sitiens U10.k10.sorted.bam U1.k10.sorted.bam U2.k10.sorted.bam U3.k10.sorted.bam U4.k10.sorted.bam U5.k10.sorted.bam U8.k10.sorted.bam U9.k10.sorted.bam W10.k10.sorted.bam W1.k10.sorted.bam W2.k10.sorted.bam W3.k10.sorted.bam W4.k10.sorted.bam W5.k10.sorted.bam W6.k10.sorted.bam W7.k10.sorted.bam W8.k10.sorted.bam W9.k10.sorted.bam
```


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

