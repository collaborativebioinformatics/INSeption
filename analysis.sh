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
bcftools query -f "%INFO/RE\t%INFO/SVTYPE\n" HG002.HiFi.GRCh37.SVLEN50.RE10.vcf > HG002.Hifi.GRCh37.table_RE_SVTYPE.tsv
Rscript scripts/make_plots_RE_and_SVTYPE.R -i HG002.Hifi.GRCh37.table_RE_SVTYPE.tsv -p plots/RE_dist_task2.pdf -q plots/SV_type_task4.pdf -r plots/INS_RE_dist_task5.pdf

"""
Extract RE of long INS, make plots for RE (task 7)
"""
bcftool query -f "%INFO/RE\n" HG002.HiFi.GRCh37.SVLEN50.RE10.largeINS.vcf > HG002.Hifi.GRCh37.table_RE_largeINS.tsv
Rscript scripts/make_plot_RE_dist.R -i ../HG002.Hifi.GRCh37.table_RE_largeINS.tsv -o plots/longIns_RE.pdf

"""
Make a plot showing the number of reads belonging to a cluster, as obtained by CARNAC-LR
"""
wc -l *.fasta > n_seqs_per_cluster.txt
sed -E 's/^[ ]*//g' n_seqs_per_cluster.txt > n_seqs_per_cluster_2.txt
script scripts/clustering_stats.R -i n_fields.txt -o plots/carnac_cluster_stats.pdf

