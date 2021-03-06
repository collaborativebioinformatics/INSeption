#!/bin/bash

error () { printf "%b\n" "[$(date)]: \e[91m$*\e[0m" >&2; exit 1; }
warn () { printf "%b\n" "[$(date)]: \e[93m$*\e[0m" >&2; }
state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
describe () { printf "%b\n" "$*" 2>&1; } 
pizzaz () { printf "%b\n" "\e[92m$*\e[0m" 2>&1; }
collect () { printf "%b\n" "\e[96m$*\e[0m" 2>&1; read var; return $var; }

stop_check() {
  if [ -z "$force" ]; then
    printf "%b\n" "\e[96m$* (y/n)?:\e[0m" 2>&1;
    read cont
    case "$cont" in
      Y|y|yes|Yes|YES|YEs|yES|yEs|yeS )
      ;;
      N|n|no|No|NO|nO )
        error "\"NO\" entered, exiting";
      ;;
      sure )
        warn "sure, I mean, whatever..."
        stop_check
      ;;
      *) 
        warn "Invalid input You entered: $cont";
	      stop_check
      ;;
    esac
  fi
}

display_header(){
  printf "%b\n" "\e[7;40m⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉    ⨉⨉  ⨉⨉⨉⨉⨉⨉⨉  ⨉⨉⨉      ⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉   ⨉⨉⨉⨉⨉⨉  ⨉⨉  ⨉⨉⨉⨉  ⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉    ⨉⨉⨉⨉⨉  ⨉⨉  ⨉⨉⨉⨉  ⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉  ⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉  ⨉⨉  ⨉⨉⨉  ⨉⨉⨉  ⨉⨉⨉⨉⨉⨉⨉⨉   ⨉⨉⨉    ⨉⨉⨉    ⨉⨉  ⨉⨉⨉   ⨉⨉⨉  ⨉ ⨉⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉  ⨉⨉⨉  ⨉⨉  ⨉⨉⨉⨉⨉  ⨉⨉⨉⨉⨉  ⨉  ⨉⨉  ⨉  ⨉⨉⨉  ⨉⨉⨉⨉⨉⨉⨉     ⨉⨉     ⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉  ⨉⨉⨉⨉  ⨉  ⨉⨉⨉⨉⨉⨉⨉  ⨉⨉⨉     ⨉⨉  ⨉  ⨉⨉⨉  ⨉⨉⨉  ⨉⨉  ⨉  ⨉⨉  ⨉  ⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉  ⨉⨉⨉⨉⨉    ⨉⨉  ⨉⨉⨉⨉  ⨉⨉  ⨉⨉⨉⨉⨉    ⨉⨉⨉⨉  ⨉⨉⨉  ⨉⨉  ⨉  ⨉⨉  ⨉  ⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉  ⨉⨉⨉  ⨉⨉⨉⨉⨉⨉   ⨉⨉  ⨉⨉⨉⨉  ⨉⨉  ⨉  ⨉⨉  ⨉⨉⨉⨉⨉⨉  ⨉⨉⨉  ⨉⨉  ⨉  ⨉⨉  ⨉  ⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉    ⨉⨉  ⨉⨉⨉⨉⨉⨉⨉  ⨉⨉⨉      ⨉⨉⨉⨉   ⨉⨉⨉  ⨉⨉⨉⨉⨉⨉   ⨉⨉  ⨉⨉⨉   ⨉⨉⨉  ⨉  ⨉⨉\e[0m"
  printf "%b\n" "\e[7;40m⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉⨉\e[0m"
}

help(){
  display_header
  describe "INSeption Version 0.1.0

  utilities 
    -h or --help
       display this message and exit

    -v or --version
       display otb version and exit

   inputs
    -b or --bam
       bam file which contains mapped reads with unaligned portions over insertions

    -v or --vcf 
       vcf file with insertions
  "
  exit 0;
}

version(){
  describe "INSeption Version 0.1.0
  polishing your mapping's structural variants with style.
  created by:
    Philippe Sanio, Wolfram Höps, Jedrzej Kubica, David Molik, Najeeb syed, and Medhat Mahmoud
  "
  exit 0;
}

while [ $# -gt 0 ] ; do
  case $1 in
    -H | --help) help ;;
    -V | --version) version ;;
    -b | --bam) BAM="$2";;
    -v | --vcf) VCF="$2";;
  esac
  shift
done


state "Extracting unaligned reads from bam file."
unaligned_reads_file=$( echo $BAM | sed 's/\.bam$/\.fastq/' )
samtools view -f 4 $BAM | samtools fastq -@ 5  > $unaligned_reads_file

state "Extracting INS that can't be resolved by Sniffles"
unresolved_ins=$( echo $VCF | sed 's/\.vcf$/\.unresolved\.vcf/' )
bcftools view -i 'INFO/SVTYPE="INS" & INFO/SVLEN=999' $VCF > $unresolved_ins

state "Extracting RE and SVTYPE, makeing plots for RE, SV type and INS-RE"
re_svtype_ins=$( echo $VCF | sed 's/\.vcf$/\.re\.sv\.tsv/' )
bcftools query -f "%INFO/RE\t%INFO/SVTYPE\n" $VCF > $re_svtype_ins
Rscript scripts/make_plots_RE_and_SVTYPE.R -i $re_svtype_ins -p plots/RE_dist.pdf -q plots/SV_type.pdf -r plots/INS_RE_dist.pdf

state "Extracting RE of long INS, make plots for RE"
unresolved_re_ins=$( echo $VCF | sed 's/\.vcf$/\.re\.sv\.unresolved\.tsv/' )
bcftool query -f "%INFO/RE\n" $unresolved_ins > $unresolved_re_ins
Rscript scripts/make_plot_RE_dist.R -i $unresolved_re_ins -o plots/longIns_RE.pdf


# TODO ##########################################
state "Clustering unmapped reads with CARNAC-LR"
unaligned_reads_file_fasta=$( echo $BAM | sed 's/\.bam$/\.fasta/' )
minimap2 $unaligned_reads_file $unaligned_reads_file -X > minimap_output.paf
python3 CARNAC-LR/scripts/paf_to_CARNAC.py minimap_output.paf $unaligned_reads_file input_carnac.txt
ulimit -s unlimited
./CARNAC-LR -f input_carnac.txt -o output_file -t 7
sed 's/@/>@/g' <(awk 'NR%3!=0' <(awk 'NR%4!=0' $unaligned_reads_file)) > $unaligned_reads_file_fasta
./scripts/CARNAC_to_fasta $unaligned_reads_file_fasta output_file 2 # NOTE: please use fasta, do not ever use fastq

state "Making a plot showing the number of reads belonging to a cluster, as obtained by CARNAC-LR"
wc -l *.fasta > n_seqs_per_cluster.txt
sed -E 's/^[ ]*//g' n_seqs_per_cluster.txt > n_seqs_per_cluster_2.txt
script scripts/clustering_stats.R -i n_fields.txt -o plots/carnac_cluster_stats.pdf
#################################################

state "Running assembly using SPAdes"
spades.py -o spades_assembly  -s $unaligned_reads_file  --only-assembler

# TODO ##########################################

# This doesn't look right
state "Getting read names for each SV (INS) and write it out in a file carrying the INS id"
bcftools query -f '%ID\t%INFO/RNAMES\n' $unresolved_ins |  while read a b ; do c=$(echo $b | tr ',' '\n'); echo $c > ./read_name_per_sv/"${a}".txt; done

state "Extracting reads from BAM file"
for file in $(ls $PWD/read_name_per_sv/*.txt); do cluter_name=$(basename $file; cluster=$(echo $cluster_name | cut -d'.' -f 1); samtools view $BAM | grep -f $file | cut -f 1,10 | awk '!/*/ {print ">"$1"\n"$2}' > "${cluster}".fasta; done

state "Running spades on every cluster.fasta"
python3 scripts/clusterAssembler_general.py multi spades spades.py fastas_to_assemble/ assemblies/

state "Getting the contigs.fasta from all successful runs, collect them in one directory." 
bash scripts/collect_assembly_fastas.sh assemblies spades_contigs_fasta

state "Concatenating contigs to use as a pseudoreference for minimap2 in a later step."
cat spades_contigs_fasta/*.fasta > spades_contigs_fasta/spades_pseudoreference.fasta
samtools faidx spades_contigs_fasta/spades_pseudoreference.fasta

state "Merging ins supporting reads, input to minimap"
cat fasta_support_ins/*.fasta > fasta_support_ins/all.fasta

state "Running minimap2: aligning all reads supporting ins against the contigs. using '-P' parameter to retain all alignments, even secondary ones." 
minimap2 -x map-hifi spades_contigs_fasta/spades_pseudoreference.fasta fasta_support_ins/all.fasta -a -P > ins_reads_vs_contigs.sam
samtools sort -O sam -o ins_reads_vs_contigs.sort.sam ins_reads_vs_contigs.sam
#################################################
