#!/bin/bash
#SBATCH --array=1-22
#SBATCH --time 2:0:0
#SBATCH --account=def-morri
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=julien.st-pierre@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/slurm-prep0-%A_%a.out

module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
caldir=/project/rpp-aevans-ab/neurohub/ukbb/genetics/cal
outdir=~/projects/def-morri/PRS

## DO NOT UNCOMMENT
# for SLURM_ARRAY_TASK_ID in {1..22}
# do
    plink \
    --memory 16000 \
    --bed $caldir/ukb_cal_chr${SLURM_ARRAY_TASK_ID}_v2.bed \
    --fam $caldir/ukb45551_cal_chr${SLURM_ARRAY_TASK_ID}_v2_s488264.fam \
    --bim $caldir/ukb_snp_chr${SLURM_ARRAY_TASK_ID}_v2.bim \
    --keep <(cat <(awk '(NR>1){print $1,$1}' $outdir/phenotypes/covar_case1.txt) <(awk '(NR>1){print $1,$1}' $outdir/phenotypes/covar_case2.txt)) \
    --remove <(awk '{print $1,$1}' $outdir/QC/exclude_missing.txt) \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --make-bed \
    --out ~/scratch/hard_calls/hard_call_chr${SLURM_ARRAY_TASK_ID}
# done
