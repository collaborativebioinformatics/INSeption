![INSeption](/.assets/INSeption.png)
polishing your mapping's structural variants with style.
---

When analyzing for structural variants to a reference genome, often a gap in the query genome will appear due to insufficient information about the query sequence. This gap leads to a known location of a structural variant, but an unknown sequence. INSeption, a bioinformatics workflow responds to this problem by utilizing the unaligned portions of reads (i.e. hanging reads) which correspond to the query genome’s missing sequence to the reference genome by building a consensus sequence from the unaligned reads of the query genome.

**What is A Structural variant?**
=================================  

Structural variants (SVs) are large genomic alterations, where large is typically (and somewhat arbitrarily) defined as encompassing at least 50 bp. These genomic variants are typically classified as deletions, duplications, insertions, inversions, and translocations describing different combinations of DNA gains, losses, or rearrangements. Copy number variations (CNVs) are a particular subtype of SVs mainly represented by deletions and duplications. SVs are typically described as single events, although more complex scenarios involving combinations of SV types exist. Chromothripsis, which is a large and complex combination of rearrangements reported in cancer, is an example. While the average genomic variation between two humans is 0.1% in terms of single nucleotide variants (SNVs), when taking SVs into account, this increases to 1.5% . In particular, telomeric regions are affected by a higher rate of SVs [[1]].

___

**Why we shall care?**
======================  

SVs are DNA rearrangements that involve at least 50 nucleotides. By virtue of their size and abundance, SVs represent an important mutational force that shape genome evolution and function, and contribute to germline and somatic diseases. The profound effect of SVs is also attributable to the numerous mechanisms by which they can disrupt protein-coding genes and cis-regulatory architecture. SVs can be grouped into mutational classes that include ‘unbalanced’ gains or losses of DNA (for example, copy-number variants, CNVs), and ‘balanced’ rearrangements that occur without corresponding dosage alterations (such as inversions and translocations). Other common forms of SVs include mobile elements that insert themselves throughout the genome, and multiallelic CNVs (MCNVs) that can exist at high copy numbers. More recently, exotic species of complex SVs have been discovered that involve two or more distinct SV signatures in a single mutational event interleaved on the same allele, and can range from CNV-flanked inversions to rare instances of localized chromosome shattering, such as chromothripsis. The diversity of SVs in humans is therefore far greater than has been widely appreciated, as is their influence on genome structure and function [[2]].





References:
==========
[1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1828-7
[2]: https://www.nature.com/articles/s41586-020-2287-8#citeas
