#!/bin/bash
#SBATCH --array=1,2
#SBATCH --time=6:0:0
#SBATCH --account=def-morri
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-user=julien.st-pierre@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/slurm-saige1-%A_%a.out

basedir=/home/gf591137/projects/def-morri/PRS

# Remove excluded subjects due to missing or duplicated
grep -v -f <(awk '{print $1}' $basedir/QC/exclude_missing.txt $basedir/QC/exclude_kinship_case${SLURM_ARRAY_TASK_ID}.txt)  $basedir/phenotypes/covar_case${SLURM_ARRAY_TASK_ID}.txt | awk 'NR>1{print $1}' > /home/gf591137/case${SLURM_ARRAY_TASK_ID}_ids

# Copy a temporary phenotype file because SAIGE can't find pheno file
cp $basedir/phenotypes/covar_case${SLURM_ARRAY_TASK_ID}.txt /home/gf591137/covar_case${SLURM_ARRAY_TASK_ID}.txt

module load singularity

singularity exec -B /home -B /scratch -B /project/rpp-aevans-ab/neurohub/ukbb/genetics:/data ~/projects/def-morri/saige.sif \
step1_fitNULLGLMM.R \
--sparseGRMFile /scratch/gf591137/sparse_grm/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
--sparseGRMSampleIDFile /scratch/gf591137/sparse_grm/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--SampleIDIncludeFile /home/gf591137/case${SLURM_ARRAY_TASK_ID}_ids \
--useSparseGRMtoFitNULL TRUE \
--skipVarianceRatioEstimation FALSE \
--plinkFile /scratch/gf591137/hard_calls/sparsekin_markers \
--phenoFile /home/gf591137/covar_case${SLURM_ARRAY_TASK_ID}.txt \
--phenoCol=case${SLURM_ARRAY_TASK_ID} \
--traitType binary \
--covarColList sex,age,batch,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--sampleIDColinphenoFile IID \
--IsOverwriteVarianceRatioFile TRUE \
--outputPrefix /scratch/gf591137/sparse_grm/case${SLURM_ARRAY_TASK_ID}_step1

rm /home/gf591137/case${SLURM_ARRAY_TASK_ID}_ids /home/gf591137/covar_case${SLURM_ARRAY_TASK_ID}.txt
