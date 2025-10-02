#!/bin/bash
#SBATCH --time 4:0:0
#SBATCH --account=def-morri
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs/slurm-prep2-%A_%a.out

# Adjusted from https://saigegit.github.io/SAIGE-doc/docs/set_step1.html#step-1-fitting-the-null-logisticlinear-mixed-model

module load plink/2.00a3.6

#(Optional) get ids for 1000 random markers for each MAC category
## calcuate allele counts for each marker in the large plink file

plink2 \
--threads 4 \
--memory 16000 \
--bfile ~/scratch/hard_calls/merged_calls \
--freq counts \
--out ~/scratch/hard_calls/prep2_freqs

## randomly extract IDs for markers falling in the two MAC categories
cat <(tail -n +2 ~/scratch/hard_calls/prep2_freqs.acount | awk '((2*$6-$5) < 20 && (2*$6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' | shuf -n 1000) \
<(tail -n +2 ~/scratch/hard_calls/prep2_freqs.acount | awk ' $5 >= 20 && (2*$6-$5)>= 20 {print $2}' | shuf -n 1000) > ~/scratch/hard_calls/samples.markerid.list


## extract markers from the large plink file
plink2 \
--threads 4 \
--memory 16000 \
--bfile ~/scratch/hard_calls/merged_calls \
--extract ~/scratch/hard_calls/samples.markerid.list \
--make-bed \
--out ~/scratch/hard_calls/sparsekin_markers


