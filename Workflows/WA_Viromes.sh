#!/bin/sh
source activate vaplinev2

DATA="/projects/mjolnir1/people/bfg522/04-Wild_Avians/04_PHAGE/01_BACREMOVED/MAGs"
VIR="/projects/mjolnir1/people/bfg522/04-Wild_Avians/04_PHAGE/VIROMES"
cd $DATA
samples=$(ls *_noBac_1.fastq.gz | sed 's/_noBac_1.fastq.gz//')

for sample in A12
  do
    # start adapter removal, using NexteraPE
    trimmomatic PE -threads 10 -phred33 "$sample"_noBac_1.fastq.gz "$sample"_noBac_1.fastq.gz "$sample"_PF1.fq "$sample"_UF1.fq "$sample"_PF2.fq "$sample"_UF2.fq ILLUMINACLIP:$CONDA_PREFIX/bin/Abin/bin/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:50
    # Match paired and un-paired
    cat "$sample"_PF1.fq "$sample"_UF1.fq > "$sample"_forward.fq
    cat "$sample"_PF2.fq "$sample"_UF2.fq > "$sample"_reverse.fq
    rm "$sample"_PF1.fq
    rm "$sample"_UF1.fq
    rm "$sample"_PF2.fq
    rm "$sample"_UF2.fq
    # remove duplicates
    seqkit rmdup "$sample"_forward.fq -s -o "$sample"_DPF1.fq -j 4
    seqkit rmdup "$sample"_reverse.fq -s -o "$sample"_DPF2.fq -j 4
    rm "$sample"_forward.fq
    rm "$sample"_reverse.fq
    # remove contaminant, using kmer, to remove phiX
    bbduk.sh in="$sample"_DPF1.fq out=$VIR/"$sample"_forward.fq  ref=$CONDA_PREFIX/bin/Abin/bin/phi_X174_phage.fa k=31 hdist=1
    bbduk.sh in="$sample"_DPF2.fq out=$VIR/"$sample"_reverse.fq  ref=$CONDA_PREFIX/bin/Abin/bin/phi_X174_phage.fa k=31 hdist=1
    rm "$sample"_DPF1.fq
    rm "$sample"_DPF2.fq
done

# Repairing
for sample in $samples
  do
    echo "$sample"
    # enter repo to process sample-specific reads (forward and reverse)
    fastq_pair "$sample"_forward.fq "$sample"_reverse.fq
    #rm "$sample"_forward.fq
    #rm "$sample"_reverse.fq
    cat "$sample"_single.fq > "$sample"_unpaired.fq
    #rm "$sample"_single.fq
done
# Single Assembly - do as job
#!/bin/sh
source activate vaplinev2
DATA="/projects/mjolnir1/people/bfg522/04-Wild_Avians/04_PHAGE/01_BACREMOVED/MAGs"
VIR="/projects/mjolnir1/people/bfg522/04-Wild_Avians/04_PHAGE/VIROMES"
cd $VIR
for sample in ${1}
  do
    #spades.py --pe1-1 "$sample"_forward.fq.paired.fq --pe1-2 "$sample"_reverse.fq.paired.fq --pe1-s "$sample"_unpaired.fq -o "$sample"_spades_folder -t 10 -m 7 --only-assembler
    megahit -1 "$sample"_forward.fq.paired.fq -2 "$sample"_reverse.fq.paired.fq --min-contig-len 2200 -t 20 --presets meta-sensitive -o "$sample"_spades_folder
    seqkit seq "$sample"_spades_folder/scaffolds.fasta -m 2200 -g > "$sample"_contigs.fasta
    rm -r "$sample"_spades_folder
    seqtk rename "$sample"_contigs.fasta contig-"$sample"_ > "$sample"_contigs_"$sample".fa
    rm "$sample"_contigs.fasta
    cat "$sample"_forward.fq.paired.fq "$sample"_reverse.fq.paired.fq  "$sample"_unpaired.fq | seqtk seq -a  | cut -d ' ' -f 1  >  tmp ; seqtk rename tmp "$sample"_ > "$sample"_reads_"$sample".fa
    #rm "$sample"*fq
    rm tmp
done

##
cp *_reads_* ./00_qual_reads
cp /*/contigs* ./00_contigs_processing
cd 00_contigs_processing/
cat *.fa > All_cross_contigs.fasta
rm *.fa
dedupe.sh in=All_cross_contigs.fasta out=temp.fa minidentity=90
seqkit sort -l -r temp.fa -o Deduplicated.fasta
rm temp.fa
seqtk rename Deduplicated.fasta vOTU_ > h_qual_contigs.fa

#VIRSORTER
#!/bin/sh
###Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00012 -A ku_00012
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N Virsorter_WA
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=20gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=2:00:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load virsorter2/2.2.3
virsorter run -w QC_VIRSORTER2_5K -i 00_contigs_processing/h_qual_contigs.fa -j 40 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.66 --min-length 5000

TAB=$'\t'
cat QC_VIRSORTER2_5K/final-viral-score.tsv | cut -f 1,8 | grep "full"  | grep "dsDNAphage"  | sed 's/|/'"${TAB}"'/g' | cut -f 1 | sed 's/_fragment/*/g' | sed 's/*/'"$TAB"'/g' | cut -f 1 > sQC_VIRSORTER01.txt
cat QC_VIRSORTER2_5K/final-viral-score.tsv | cut -f 1,8 | grep "full"  | grep "RNA"  | sed 's/|/'"${TAB}"'/g' | cut -f 1        | sed 's/_fragment/*/g' | sed 's/*/'"$TAB"'/g' | cut -f 1 > sQC_VIRSORTER02.txt
cat QC_VIRSORTER2_5K/final-viral-score.tsv | cut -f 1,8 | grep "full"  | grep "ssDNA"  | sed 's/|/'"${TAB}"'/g' | cut -f 1      | sed 's/_fragment/*/g' | sed 's/*/'"$TAB"'/g' | cut -f 1 > sQC_VIRSORTER03.txt
cat QC_VIRSORTER2_5K/final-viral-score.tsv | cut -f 1,8 | grep "full"  | grep "lavidaviridae"  | sed 's/|/'"${TAB}"'/g' | cut -f 1      | sed 's/_fragment/*/g' | sed 's/*/'"$TAB"'/g' | cut -f 1 > sQC_VIRSORTER04.txt
cat QC_VIRSORTER2_5K/final-viral-score.tsv | cut -f 1,8 | grep "full"  | grep "NCLDV"  | sed 's/|/'"${TAB}"'/g' | cut -f 1      | sed 's/_fragment/*/g' | sed 's/*/'"$TAB"'/g' | cut -f 1 > sQC_VIRSORTER05.txt
cat sQC* | cut -f1 | sort -u > h_qual_clean_contigs.txt

seqtk subseq h_qual_contigs.fa h_qual_clean_contigs.txt | seqkit sort -l -r -o h_qual_clean_contigs.fa

module load usearch/10.0.240
usearch -ublast h_qual_clean_contigs.fa -db ../bin/vog.lca_206_anello.fa -evalue 1e-3  -blast6out blastout.txt -cover_query  -dbmask fastamino -threads 40
cut -f1 < blastout.txt >  tmp1.txt
TAB=$'\t'
cat  blastout.txt | cut -f2  | sed 's/_/'"$TAB"'/g' | cut -f1  >  tmp2.txt
paste tmp1.txt tmp2.txt > blastout_selected.txt
rm tmp*
