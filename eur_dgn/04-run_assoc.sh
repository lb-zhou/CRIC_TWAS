#!/bin/bash


chr=$1
phe=$2 #number of column in the phenotype file
phename=$3
out=$4

# Add module for R
module add r/3.4.1

Rscript /proj/yunligrp/users/munan/scripts/assoc_0.3.R  -e /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/eur_dgn/prediction/prediction.chr${chr} \
                        -k /proj/yunligrp/users/jwen/CRIC/assoc/CRIC_kinf_plain.kindump \
                        -p /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/pheno/EUR/${phename}_pheno_eur.txt \
                        -c ${phe} -o /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/eur_dgn/assoc/${out}

