#!/bin/bash
conda activate /projects/mjolnir1/apps/conda/anvio-7
REF_DIR='/projects/mjolnir1/people/bfg522/00_DATA/C2/4-Wild_Avians/refs'
FASTA_DIR='/projects/mjolnir1/people/bfg522/00_DATA/C2/4-Wild_Avians/00_RawData/GreylagGoose'
WORK_DIR='/projects/mjolnir1/people/bfg522/00_DATA/C2/4-Wild_Avians/01_Filtered'
cd $FASTA_DIR/
find *.fastq.gz > temp
sed 's/_[1-2].fastq.gz//g' temp > temp2
uniq temp2 > temp3
echo 'sample' > header
cat header temp3 > samples.txt
rm -f temp*
# A simple loop to serially map all samples.
# referenced from within http://merenlab.org/tutorials/assembly_and_mapping/
# modified by @Jacob Agerbo Rasmussen
# how many threads should each mapping task use?
NUM_THREADS=10 # Number of CPU for processing
REF='Anser_indicus' #name of Bowtie2 DB
for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    echo "$sample"
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls $FASTA_DIR/"$sample"_1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls $FASTA_DIR/"$sample"_2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    #map with bowtie2
    bowtie2 --threads $NUM_THREADS -x $REF_DIR/$REF -1 $R1s -2 $R2s -S $WORK_DIR/"$sample".sam
    #convert sam to bam
    samtools view -@ $NUM_THREADS -bS $WORK_DIR/"$sample".sam > $WORK_DIR/"$sample".bam
    rm $WORK_DIR/"$sample".sam
    # keep only unmapped reads
    samtools view -@ $NUM_THREADS  -b -f 12 -F 256 $WORK_DIR/"$sample".bam > $WORK_DIR/"$sample"_noHost.bam
    # split to paired reads
    samtools sort -n -m 5G -@ 10 $WORK_DIR/"$sample"_noHost.bam -o $WORK_DIR/"$sample"_noHost_sorted.bam
    samtools fastq -@ 10 $WORK_DIR/"$sample"_noHost_sorted.bam \
    -1 $WORK_DIR/"$sample"_noHost_1.fastq.gz \
    -2 $WORK_DIR/"$sample"_noHost_2.fastq.gz \
    -0 /dev/null -s /dev/null -n
    # remove bam files
    rm $WORK_DIR/*.bam
done
