for casenum in {1..2}
do
echo case$casenum
~/projects/def-morri/combine_without_header.sh \
	/scratch/gf591137/saige2_outputs/saige2_case${casenum}_chr* \
	>> /scratch/gf591137/saige2_outputs/saige2_case${casenum}_full
done
