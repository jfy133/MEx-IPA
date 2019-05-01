# Test Data Generation

Data is a subset of high-throughput sequencing of ancient calculus samples taken from Weyrich et al. (2017) Nature.

The raw sequencing files were downloaded from OAGR (https://www.oagr.org.au/). Adapter clipping, barcode clipping and merging was performed using AdapterRemoval2. Human reads were removed by mapping to the HG19 reference genome with `bwa` and `samtools`, and unmapped reads exported and ran through MALT v040 to the NCBI nucleotide 'nt' database of December 2017 with 85% precent identity. The resulting RMA6 files were then processed in MALT-Extract as implemented in the 'HOPS' pipeline v0.1, with the 'dietary' taxa supposedly detected in Weyrich et al. as the taxa to extract.
