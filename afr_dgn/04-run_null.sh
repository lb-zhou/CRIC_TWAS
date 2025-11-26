#!/bin/bash


phe=$1 #number of column in the phenotype file
covar=$2
phename=$3


# Add module for R
module add r

rename .ENSG _ENSG /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/prediction/prediction.chr${chr}/*

Rscript /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/scripts/null_model_0.1.R \
                        -k /proj/yunligrp/users/jwen/CRIC/assoc/CRIC_kinf_plain.kindump \
                        -p /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/pheno/${phename}_AFR.ped \
                        -c ${phe} -x ${covar} \
                        -o /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/afr_dgn/assoc/${phename}
