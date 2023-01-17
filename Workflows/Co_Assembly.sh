#!/bin/bash
### Assembly
conda activate WORK
FASTA_DIR='/projects/mjolnir1/people/bfg522/00_DATA/C2/4-Wild_Avians/01_Filtered/GreylagGoose'
ASSEM_DIR='/projects/mjolnir1/people/bfg522/00_DATA/C2/4-Wild_Avians/02_Assembly'
FASTA=`ls $FASTA_DIR/*.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
MINLENGTH=2000 #minimum length of contigs for downstream analysis
NUM_THREADS=10 # Number of CPU for processing
megahit -r $FASTA --min-contig-len $MINLENGTH -t $NUM_THREADS --presets meta-sensitive -o $ASSEM_DIR/GreylagGoose
