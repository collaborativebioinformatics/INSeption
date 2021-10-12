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

SVs are DNA rearrangements that involve at least 50 nucleotides. By virtue of their size and abundance, SVs represent an important mutational force that shape genome evolution and function, and contribute to germline and somatic diseases. The profound effect of SVs is also attributable to the numerous mechanisms by which they can disrupt protein-coding genes and cis-regulatory architecture. SVs can be grouped into mutational classes that include ‘unbalanced’ gains or losses of DNA (for example, copy-number variants, CNVs), and ‘balanced’ rearrangements that occur without corresponding dosage alterations (such as inversions and translocations). Other common forms of SVs include mobile elements that insert themselves throughout the genome, and multiallelic CNVs (MCNVs) that can exist at high copy numbers. More recently, exotic species of complex SVs have been discovered that involve two or more distinct SV signatures in a single mutational event interleaved on the same allele, and can range from CNV-flanked inversions to rare instances of localized chromosome shattering, such as chromothripsis. The diversity of SVs in humans is therefore far greater than has been widely appreciated, as is their influence on genome structure and function [[2]]. Insertions plays an important role in human diseases, Almost half of mammalian genomes are comprised of a class of repeat DNA sequences known as transposable elements (TEs). Over evolutionary time, the dynamic nature of a genome is driven, in part, by the activity of transposable elements (TE) such as retrotransposons. On a shorter time scale it has been established that new TE insertions can result in single-gene disease in an individual. In humans, the non-LTR retrotransposon Long INterspersed Element-1 (LINE-1 or L1) is the only active autonomous TE. In addition to mobilizing its own RNA to new genomic locations via a "copy-and-paste" mechanism, LINE-1 is able to retrotranspose other RNAs including Alu, SVA, and occasionally cellular RNAs. To date in humans, 124 LINE-1-mediated insertions which result in genetic diseases have been reported. Disease causing LINE-1 insertions have provided a wealth of insight and the foundation for valuable tools to study these genomic parasites [[4]].


___

**How to detect?**
==================

Algorithms detect SVs from long-reads data by leveraging intra-read and inter-read signatures. Intra-read signatures enable the direct detection of SVs and are derived from reads spanning entire SV events, resulting in a missing/inserted sequence (deletion) or a soft-clip (insertion) within properly aligned flanking sequences. After signature detection, callers typically cluster and merge similar signatures from multiple reads, delineate proximal but different signatures and choose the highest quality reads that support the putative SV [[3]] (**Figure 1**).  


![INSeption](/.assets/ins.png)  

*Figure 1*: A digram shows the signals of an insertion that SV callers uses, reference shown in green, the sample in blue, and the insertion in orange color. we show three reads in red one covers the whole SV and as a result there is part of the read was not aligned (dotted), other two reads are flanking the SV the unaligned part is shown as dotted line.

___

**The problem**
===============

For an insertion (INS), if the read length used to detect the INS is larger than the insertion; then the SV caller will infer the inserted sequence from it, but this does not happen all the time, some insertions are large that there is no read can span the whole insertions. Consecutively, the SV caller will identify the insertion location on reference but will not inform the sequence of that specific insertions (**Figure 2**).  

![INSeption](/.assets/ins2.png)

*Figure 2*: A digram shows reads that signals an insetion but no reads are long enough to span the entire event. we show four reads in red two dotted reads full in the insertion, other two reads are flanking the SV the unaligned part is shown as dotted line.


___

**Our approach to solve the problem**
=====================================

## Data summary

We started by aligning reads to GRCh37 using minimap2, identified SVs using Sniffles 1.0.12, we filter out SVs that are supported by less than 10 reads using bcftools 1.12 where we detected 21,248 SVs in (**Figure 3**) we show the distribution of SVs that we identified. Likewise, we plotted the distribution of reads number supporting each type of SVs we detected (**Figure 4**), as well as the allele frequency (**Figure 5**).

![SV types](/plots/SV_distribution.png)

*Figure 3*: Frequency of identified SVs, on the X axis the SV type, on the Y axis the count of identified SVs.  



![SV types](/plots/re_support.png)

*Figure 4*: distribution of reads count supporting each SV type.  



![SV types](/plots/af.png)

*Figure 5*:Allele frequency distribution across different types of SVs.  



## Workflow

![INSeption](/.assets/workflow.png)
![INSeption](/.assets/workflow_ex.png)




___

**Installation**
=================



___

**How to run**
==============


[1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1828-7
[2]: https://www.nature.com/articles/s41586-020-2287-8#citeas
[3]: https://www.nature.com/articles/s41576-019-0180-9#citeas
[4]: https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-016-0065-9

References:
==========
1. Mahmoud, M., Gobet, N., Cruz-Dávalos, D.I. et al. Structural variant calling: the long and the short of it. Genome Biol 20, 246 (2019). https://doi.org/10.1186/s13059-019-1828-7
2. Collins, R.L., Brand, H., Karczewski, K.J. et al. A structural variation reference for medical and population genetics. Nature 581, 444–451 (2020). https://doi.org/10.1038/s41586-020-2287-8
3. Ho, S.S., Urban, A.E. & Mills, R.E. Structural variation in the sequencing era. Nat Rev Genet 21, 171–189 (2020). https://doi.org/10.1038/s41576-019-0180-9
4. Hancks, D.C., Kazazian, H.H. Roles for retrotransposon insertions in human disease. Mobile DNA 7, 9 (2016). https://doi.org/10.1186/s13100-016-0065-9
