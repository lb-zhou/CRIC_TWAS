#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=240:00:00
#SBATCH --mem=24GB
#SBATCH --array=1-22


echo "My SLURM_ARRAY_TASK_ID: " ${SLURM_ARRAY_TASK_ID}
module add samtools


subset=1
if [ ${subset} == 1 ]
then
    wkdir="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/intersection/vcf"
    cric="/proj/yunligrp/users/jwen/CRIC/impute_MIS/AFR_r2_hat_v2"
    pos="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/pos"
    bcftools view -R ${pos}/chr${SLURM_ARRAY_TASK_ID}.pos.txt ${cric}/chr${SLURM_ARRAY_TASK_ID}.dose.AFR.snpid.v2.vcf.gz -Oz -o ${wkdir}/cric.chr${SLURM_ARRAY_TASK_ID}.sub.vcf.gz
    bcftools index -t -f ${wkdir}/cric.chr${SLURM_ARRAY_TASK_ID}.sub.vcf.gz
fi

module add vcftools
filter_snp=1
if [ ${filter_snp} == 1 ]
then
    snp="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/pos"
    wkdir="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/intersection"
    vcftools --gzvcf ${wkdir}/vcf/cric.chr${SLURM_ARRAY_TASK_ID}.sub.vcf.gz \
             --snps ${snp}/chr${SLURM_ARRAY_TASK_ID}.intersection.txt \
             --recode --recode-INFO-all --stdout | bgzip -c > ${wkdir}/chr${SLURM_ARRAY_TASK_ID}.subset1.vcf.gz
    bcftools index -t -f ${wkdir}/chr${SLURM_ARRAY_TASK_ID}.subset1.vcf.gz
fi