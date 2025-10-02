#!/bin/bash
#SBATCH --time=12:0:0
#SBATCH --account=def-morri
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mail-user=julien.st-pierre@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/slurm-saige0-%A_%a.out

module load singularity

singularity exec -B /home -B /scratch -B /project/rpp-aevans-ab/neurohub/ukbb/genetics:/data ~/projects/def-morri/saige.sif \
createSparseGRM.R \
--plinkFile /scratch/gf591137/hard_calls/merged_calls \
--nThreads=32 \
--numRandomMarkerforSparseKin=5000 \
--relatednessCutoff=0.125 \
--outputPrefix /scratch/gf591137/sparse_grm/sparseGRM
