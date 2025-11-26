#!/bin/bash

module add r
module add samtools

# Weights path #
path_to_main="/proj/yunligrp/users/jwen/minority_TWAS/update_UKBb/MESA_AA"
out_res="/proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_mesa"

#input variable
chr=$1

mkdir ${out_res}/prediction/temp
mkdir ${out_res}/prediction/temp/chunking
mkdir ${out_res}/prediction/temp/chunking/chunking.chr${chr}
mkdir ${out_res}/prediction/temp/DC
mkdir ${out_res}/prediction/temp/DC/DC.chr${chr}
mkdir ${out_res}/prediction/prediction.chr${chr}


#Loop through genes on list
while read -r chrom start stop commonID; do
    #finding the right chunk
    for i in {1..42}; do
	if [[ $(( $start + 1000000 )) -gt $(( ($i-1)*6000000 + 1 )) && $(( $stop + 1000000 )) -lt $(( ($i-1)*6000000 + 10000000 ))  ]]; then
	                                    chunk=$i
					                    break
							                    fi
	done
    echo ${commonID} ${chunk}

    #extract SNPs in the training model
    sed 1d ${path_to_main}/training/training.chr${chr}/MESA_AA.${chr}.${commonID}.betas.EN.txt | cut -f 1,2 > \
        ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.snppos.temp

    bcftools view -R ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.snppos.temp \
	${out_res}/prediction/chunking/chunking.chr${chr}/cric.chr${chr}.chunk${chunk}.vcf.gz \
	-Oz -o ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.vcf.gz

#Dosage convertor and reformat DCoutput
    /proj/yunligrp/users/jdrosen/bin/DosageConvertor/bin/DosageConvertor --vcfDose ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.vcf.gz \
	--prefix ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC \
	--type mach --format 1

    paste <(zcat ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC.mach.dose.gz | cut -f 1 | awk 'BEGIN{FS=">"}{print $2}') <(zcat ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC.mach.dose.gz | cut -f3-) | gzip \
	> ${out_res}/prediction/temp/DC/DC.chr${chr}/reformat.${commonID}.chr${chr}.dosage.gz

    sed 1d  ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC.mach.info | cut -f1 | awk 'BEGIN {FS=":"; OFS="\t"} {print $1, $2, $3, $4}' \
	      > ${out_res}/prediction/temp/DC/DC.chr${chr}/snplist.chr${chr}.${commonID}.temp

#prediction

    Rscript /proj/yunligrp/users/munan/scripts/pred_0.6.R -d ${out_res}/prediction/temp/DC/DC.chr${chr}/reformat.${commonID}.chr${chr}.dosage.gz \
	-s ${out_res}/prediction/temp/DC/DC.chr${chr}/snplist.chr${chr}.${commonID}.temp \
	-b ${path_to_main}/training/training.chr${chr}/MESA_AA.${chr}.${commonID}.betas.EN.txt \
	-l ${path_to_main}/training/training.chr${chr}/MESA_AA.${chr}.${commonID}.log \
	-o ${out_res}/prediction/prediction.chr${chr}/${commonID}.predicted

    rm ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.snppos.temp
    rm ${out_res}/prediction/temp/chunking/chunking.chr${chr}/chr${chr}.${commonID}.vcf.gz
    rm ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC.mach.dose.gz
    rm ${out_res}/prediction/temp/DC/DC.chr${chr}/chr${chr}.${commonID}.DC.mach.info
    rm ${out_res}/prediction/temp/DC/DC.chr${chr}/snplist.chr${chr}.${commonID}.temp
    rm ${out_res}/prediction/temp/DC/DC.chr${chr}/reformat.${commonID}.chr${chr}.dosage.gz

done  <  gencode.chr${chr}.probeIDs_edit.temp
