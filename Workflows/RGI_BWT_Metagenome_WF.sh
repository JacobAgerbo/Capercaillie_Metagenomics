#!/bin/bash
module load miniconda
conda activate rgi
DATA="/projects/mjolnir1/people/bfg522/04-Wild_Avians/01_Filtered/Capercaillie"
ARGs="/projects/mjolnir1/people/bfg522/04-Wild_Avians/05_ARGs"
DB="/projects/mjolnir1/people/bfg522/04-Wild_Avians/05_ARGs/.bin"
cd $DATA
samples=$(ls *_metaG_1.fq.gz | sed 's/_metaG_1.fq.gz//')
cd $ARGs
rgi load \
  --card_json .bin/card.json \
  --debug --local \
  --card_annotation .bin/card_database_v3.2.5.fasta \
  --wildcard_annotation .bin/wildcard_database_v4.0.0.fasta \
  --wildcard_index .bin/wildcard/index-for-model-sequences.txt \
  --wildcard_version 4.0.0 \
  --amr_kmers .bin/wildcard/all_amr_61mers.txt \
  --kmer_database .bin/wildcard/61_kmer_db.json \
  --kmer_size 61
for sample in N10
  do
    rgi bwt --read_one $DATA/"$samples"_metaG_1.fq.gz --read_two $DATA/"$samples"_metaG_2.fq.gz --output_file $ARGs/"$samples" --local -n 20 --include_wildcard
done  
