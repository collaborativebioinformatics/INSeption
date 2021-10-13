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

"""
Clustering unmapped reads with CARNAC-LR
"""
minimap2 HG002.GRCh37.unmapped.fastq HG002.GRCh37.unmapped.fastq -X > minimap_output.paf
python3 CARNAC-LR/scripts/paf_to_CARNAC.py minimap_output.paf HG002.GRCh37.unmapped.fastq input_carnac.txt
ulimit -s unlimited
./CARNAC-LR -f input_carnac.txt -o output_file -t 7
sed 's/@/>@/g' <(awk 'NR%3!=0' <(awk 'NR%4!=0' HG002.GRCh37.unmapped.fastq)) > HG002.GRCh37.unmapped.fasta
./scripts/CARNAC_to_fasta HG002.GRCh37.unmapped.fasta output_file 2 # NOTE: please use fasta, do not ever use fastq

"""
Running assembly using flye flye/2.8.1
"""
flye --pacbio-hifi HG002.GRCh37.unmapped.fastq.gz  --out-dir assembly --threads 32

"""
Running assembly using SPAdes/3.14.0
"""
spades.py -o spades_assembly  -s ../../HG002.GRCh37.unmapped.fastq.gz  --only-assembler


"""
Assemble cluster using flye
"""
flye --pacbio-hifi  cluster62.fasta  --out-dir 62


"""
Get reade name for each SV (INS) and write it out in a file carrying the INS id
"""
bcftools query -f '%ID\t%INFO/RNAMES\n' HG002.HiFi.GRCh37.SVLEN50.RE10.largeINS.vcf |  while read a b ; do c=$(echo $b | tr ',' '\n'); echo $c > ./read_name_per_sv/"${a}".txt ; done

"""
Extract reads from BAM file
"""
for i in $(ls $PWD/read_name_per_sv/*.txt); do a=$(basename $i); c=$(echo $a | cut -d'.' -f 1); samtools view my.bam | grep -f $i | cut -f 1,10 | awk '!/*/ {print ">"$1"\n"$2}' > "${c}".fasta


"""
Run spades (much faster than flye) on every cluster.fasta
"""
python3 scripts/clusterAssembler_general.py multi spades spades.py fastas_to_assemble/ assemblies/

