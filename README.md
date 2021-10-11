![INSeption](/.assets/INSeption.png)
polishing your mapping's structural variants with style.
---

When analyzing for structural variants to a reference genome, often a gap in the query genome will appear due to insufficient information about the query sequence. This gap leads to a known location of a structural variant, but an unknown sequence. INSeption, a bioinformatics workflow responds to this problem by utilizing the unaligned portions of reads (i.e. hanging reads) which correspond to the query genome’s missing sequence to the reference genome by building a consensus sequence from the unaligned reads of the query genome.
