"""
Extract unaligned reads from bam file.
"""
samtools view -f 4 my.bam | samtools fastq -@ 5  | bgzip > HG002.GRCh37.unmapped.fastq.gz


"""
Extract INS that can't be resolved by Sniffles 1.0.12.
"""
bcftools view -i 'INFO/SVTYPE="INS" & INFO/SVLEN=999' HG002.HiFi.GRCh37.SVLEN50.RE10.vcf > HG002.HiFi.GRCh37.SVLEN50.RE10.largeINS.vcf


"""
Extract RE and SVTYPE, make plots for RE (task 2), SV type (task 4) and INS-RE (task 5)
"""
bcftools query -f "%INFO/RE\t%INFO/SVTYPE\n" HG002.HiFi.GRCh37.SVLEN50.RE10.vcf > RE_SVTYPE.tsv
Rscript scripts/make_plots_RE_and_SVTYPE.R -i RE_SVTYPE.tsv -p plots/RE_dist_task2.pdf -q plots/SV_type_task4.pdf -r plots/INS_RE_dist_task5.pdf
