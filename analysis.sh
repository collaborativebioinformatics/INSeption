"""
Extract unaligned reads from bam file.
"""
samtools view -f 4 my.bam | samtools fastq -@ 5  | bgzip > HG002.GRCh37.unmapped.fastq.gz
