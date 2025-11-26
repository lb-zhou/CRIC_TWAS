#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=16GB
#SBATCH --array=1-22

echo "My SLURM_ARRAY_TASK_ID: " ${SLURM_ARRAY_TASK_ID}
chr=${SLURM_ARRAY_TASK_ID}

module add samtools
module add datamash

mkdir /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_mesa/prediction/chunking/chunking.chr${chr}/


tabix -p vcf /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_mesa/intersection/chr${chr}.subset1.vcf.gz

chunking=1
if [ $chunking == 1 ]
then
    WKDIR="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_mesa/prediction"
    input="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_mesa/intersection"

    for i in {1..42}; do
        bcftools view -r chr${chr}:$(( ($i-1)*6000000+1 ))-$(( $i*6000000+4000000 )) ${input}/chr${chr}.subset1.vcf.gz -Oz -o ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz

        tabix -p vcf ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz

        /proj/yunligrp/users/jdrosen/bin/DosageConvertor/bin/DosageConvertor --vcfDose ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz \
                                                                             --prefix ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dc --type mach --format 1

        paste <(zcat ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dc.mach.dose.gz | cut -f 1 | awk 'BEGIN{FS=">"}{print $2}') <(zcat ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dc.mach.dose.gz | cut -f 3-) | gzip > ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dosage.gz

        #Remove files from empty chunks
        n_lines="$(zcat ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz | grep -v "^##" | wc -l)"
        echo chunk ${i}: $n_lines
        if [ "$n_lines" -le 1 ]
        then
            rm ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz
            rm ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.vcf.gz.tbi
            rm ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dc.mach.dose.gz
            rm ${WKDIR}/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${i}.dosage.gz
        fi
    done
fi
