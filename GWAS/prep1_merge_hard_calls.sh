#!/bin/bash
#SBATCH --time 12:0:0
#SBATCH --account=def-morri
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=julien.st-pierre@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/slurm-prep1-%A_%a.out

module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

echo ~/scratch/hard_calls/hard_call_chr1 > ~/scratch/hard_calls_merge_list.txt
for i in {2..22} 
do
    echo ~/scratch/hard_calls/hard_call_chr${i} >> ~/scratch/hard_calls_merge_list.txt
done

plink \
--threads 8 \
--memory 64000 \
--merge-list ~/scratch/hard_calls_merge_list.txt \
--make-bed \
--out ~/scratch/hard_calls/merged_calls

rm ~/scratch/hard_calls_merge_list.txt
