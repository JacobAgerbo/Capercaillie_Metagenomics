conda activate WORK
for host in Chicken RuddyShelduck BH_Goose Capercaillie ALL
  do
    anvi-script-reformat-fasta "$host"/final.contigs.fa -o "$host"/contigs-fixed.fa -l 2000 --simplify-names
    anvi-gen-contigs-database -f "$host"/contigs-fixed.fa -o "$host"/CONTIGS.db -n 'An contigs database of "$sample"'
    anvi-run-hmms -c "$host"/CONTIGS.db
    bowtie2-build "$host"/contigs-fixed.fa "$host"/contigs
done

FASTA_DIR='/projects/mjolnir1/people/bfg522/04-Wild_Avians/01_Filtered'
WORK_DIR='/projects/mjolnir1/people/bfg522/04-Wild_Avians/03_Anvio'
NUM_THREADS=10
for host in RuddyShelduck BH_Goose
  do
    cd $FASTA_DIR/"$host"
    find *.fastq.gz > temp
    sed 's/_noHost_[1-2].fastq.gz//g' temp > temp2
    uniq temp2 > temp3
    echo 'sample' > header
    cat header temp3 > samples.txt
    rm -f temp*
    cd $FASTA_DIR
    for sample in `awk '{print $1}' $FASTA_DIR/"$host"/samples.txt`
      do
    if [ "$sample" == "sample" ]; then continue; fi
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls $FASTA_DIR/"$host"/"$sample"_noHost_1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls $FASTA_DIR/"$host"/"$sample"_noHost_2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    echo "$host"
    echo "$sample"
    bowtie2 --threads $NUM_THREADS -x $WORK_DIR/ALL/contigs -1 $R1s -2 $R2s --no-unal -S $WORK_DIR/"$sample".sam
    samtools view -@ $NUM_THREADS -F 4 -bS $WORK_DIR/"$sample".sam > $WORK_DIR/"$sample"-RAW.bam
    anvi-init-bam $WORK_DIR/"$sample"-RAW.bam -o $WORK_DIR/"$sample".bam
    rm $WORK_DIR/"$sample".sam $WORK_DIR/"$sample"-RAW.bam
    done
done
# Since Chicken data was single read I redo the loop with other parameters
NUM_THREADS=4
FASTA_DIR='/projects/mjolnir1/people/bfg522/04-Wild_Avians/01_Filtered'
WORK_DIR='/projects/mjolnir1/people/bfg522/04-Wild_Avians/03_Anvio'
for host in Chicken
  do
    cd $FASTA_DIR/"$host"
    find *.fastq.gz > temp
    sed 's/.fastq.gz//g' temp > temp2
    uniq temp2 > temp3
    echo 'sample' > header
    cat header temp3 > samples.txt
    rm -f temp*
    cd $FASTA_DIR
    for sample in `awk '{print $1}' $FASTA_DIR/"$host"/samples.txt`
      do
    if [ "$sample" == "sample" ]; then continue; fi
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    READs=`ls $FASTA_DIR/"$host"/"$sample".fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    echo "$host"
    echo "$sample"
    bowtie2 --threads $NUM_THREADS -x $WORK_DIR/ALL/contigs -U $READs --no-unal -S $WORK_DIR/"$sample".sam
    samtools view -@ $NUM_THREADS -F 4 -bS $WORK_DIR/"$sample".sam > $WORK_DIR/"$sample"-RAW.bam
    anvi-init-bam $WORK_DIR/"$sample"-RAW.bam -o $WORK_DIR/"$sample".bam
    rm $WORK_DIR/"$sample".sam #$WORK_DIR/"$sample"-RAW.bam
    done
done
##
for host in Capercaillie
  do
    cd $FASTA_DIR/"$host"
    find *.fq.gz > temp
    sed 's/_metaG_[1-2].fq.gz//g' temp > temp2
    uniq temp2 > temp3
    echo 'sample' > header
    cat header temp3 > samples.txt
    rm -f temp*
    cd $FASTA_DIR
    for sample in `awk '{print $1}' $FASTA_DIR/"$host"/samples.txt`
      do
    if [ "$sample" == "sample" ]; then continue; fi
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls $FASTA_DIR/"$host"/"$sample"_metaG_1.fq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls $FASTA_DIR/"$host"/"$sample"_metaG_2.fq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    echo "$host"
    echo "$sample"
    bowtie2 --threads $NUM_THREADS -x $WORK_DIR/ALL/contigs -1 $R1s -2 $R2s --no-unal -S $WORK_DIR/"$sample".sam
    samtools view -@ $NUM_THREADS -F 4 -bS $WORK_DIR/"$sample".sam > $WORK_DIR/"$sample"-RAW.bam
    anvi-init-bam $WORK_DIR/"$sample"-RAW.bam -o $WORK_DIR/"$sample".bam
    rm $WORK_DIR/"$sample".sam $WORK_DIR/"$sample"-RAW.bam
    done
done
