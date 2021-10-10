"""
Extract unaligned reads from bam file.
"""
samtools view -f 4 my.bam | samtools fastq -@ 5  | bgzip > HG002.GRCh37.unmapped.fastq.gz


"""
Extract INS that can't be resolved by Sniffles 1.0.12.
"""
bcftools view -i 'INFO/SVTYPE="INS" & INFO/SVLEN=999' HG002.HiFi.GRCh37.SVLEN50.RE10.vcf > HG002.HiFi.GRCh37.SVLEN50.RE10.largeINS.vcf
