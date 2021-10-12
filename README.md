![INSeption](/.assets/INSeption.png)
polishing your mapping's structural variants with style.
---

When analyzing for structural variants to a reference genome, often a gap in the query genome will appear due to insufficient information about the query sequence. This gap leads to a known location of a structural variant, but an unknown sequence. INSeption, a bioinformatics workflow responds to this problem by utilizing the unaligned portions of reads (i.e. hanging reads) which correspond to the query genome’s missing sequence to the reference genome by building a consensus sequence from the unaligned reads of the query genome.

__What is A Structural variant?__
Structural variants (SVs) are large genomic alterations, where large is typically (and somewhat arbitrarily) defined as encompassing at least 50 bp. These genomic variants are typically classified as deletions, duplications, insertions, inversions, and translocations describing different combinations of DNA gains, losses, or rearrangements. Copy number variations (CNVs) are a particular subtype of SVs mainly represented by deletions and duplications (reviewed in Carvalho and Lupski). SVs are typically described as single events, although more complex scenarios involving combinations of SV types exist. Chromothripsis, which is a large and complex combination of rearrangements reported in cancer, is an example. While the average genomic variation between two humans is 0.1% in terms of single nucleotide variants (SNVs), when taking SVs into account, this increases to 1.5% . In particular, telomeric regions are affected by a higher rate of SVs.

___
