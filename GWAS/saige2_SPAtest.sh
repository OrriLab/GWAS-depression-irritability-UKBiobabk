#!/bin/bash
#SBATCH --array=1-22
#SBATCH --time 24:0:0
#SBATCH --account=def-morri
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=julien.st-pierre@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/slurm-saige2_c1-%A_%a.out

module load singularity

singularity exec -B /home -B /scratch -B /project/rpp-aevans-ab/neurohub/ukbb/genetics:/data ~/projects/def-morri/saige.sif \
step2_SPAtests.R \
--bedFile=/scratch/gf591137/bims/ukbb_chr${SLURM_ARRAY_TASK_ID}.bed \
--bimFile=/scratch/gf591137/bims/ukbb_chr${SLURM_ARRAY_TASK_ID}.bim \
--famFile=/scratch/gf591137/bims/ukbb_chr${SLURM_ARRAY_TASK_ID}.fam \
--AlleleOrder alt-first \
--LOCO FALSE \
--chrom ${SLURM_ARRAY_TASK_ID} \
--is_output_moreDetails TRUE \
--is_fastTest TRUE \
--sparseGRMFile /scratch/gf591137/sparse_grm/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
--sparseGRMSampleIDFile /scratch/gf591137/sparse_grm/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--GMMATmodelFile /scratch/gf591137/sparse_grm/${case}_step1.rda \
--varianceRatioFile /scratch/gf591137/sparse_grm/${case}_step1.varianceRatio.txt \
--is_overwrite_output FALSE \
--SAIGEOutputFile /scratch/gf591137/saige2_outputs/saige2_${case}_chr${SLURM_ARRAY_TASK_ID}

## DO NOT UNCOMMENT
## parallel sbatch ::: --export=case={case1,case2} ::: saige2_SPAtest.sh