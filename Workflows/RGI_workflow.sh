rgi card_annotation -i ./card.json > card_annotation.log 2>&1
rgi load -i ./card.json --card_annotation card_database_v3.0.1.fasta --local


rgi wildcard_annotation -i wildcard --card_json ./card.json \
  -v v3.2.3 > wildcard_annotation.log 2>&1
rgi load --wildcard_annotation wildcard_database_v3.2.3.fasta \
  --wildcard_index ./wildcard/index-for-model-sequences.txt \
  --card_annotation card_database_v3.2.3.fasta


#
rgi load --kmer_database ./wildcard/61_kmer_db.json \
  --amr_kmers ./wildcard/all_amr_61mers.txt --kmer_size 61 \
  --debug > kmer_load.61.log 2>&1


#
#!/bin/bash
DIR=`find ./DATA -mindepth 0 -type d`
for D in $DIR; do
      NAME=$(basename $D);
      parallel --no-notice --progress -j+0 'rgi main -i {} -o {.} --num_threads 6 --split_prodigal_jobs --clean --debug > {.}.log 2>&1' ::: $NAME/*.{fa,fasta};
done

for sample in *.fa
  do
    echo "$sample"
    rgi main --input_sequence "$sample" --output_file RESULTS/"$sample" --clean --num_threads 6 --split_prodigal_jobs
    mv "$sample" ./done/
done
