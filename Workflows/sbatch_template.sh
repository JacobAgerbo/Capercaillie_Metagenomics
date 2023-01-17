sbatch --job-name remove_host_HI \
       --mail-type ALL \
       --mail-user jacob.rasmussen@bio.ku.dk \
       --time 5-00:00:00 \
       --cpus-per-task 20 \
       --mem-per-cpu 30 \
       -- Co_Assembly_MD.sh
