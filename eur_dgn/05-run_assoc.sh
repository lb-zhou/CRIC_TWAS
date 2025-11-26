#!/bin/bash

chr=$1
phename=$2


# Add module for R
module add r

rename .ENSG _ENSG /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/eur_dgn/prediction/prediction.chr${chr}/*

Rscript /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/scripts/assoc_0.2.R \
        -e /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/eur_dgn/prediction/prediction.chr${chr} \
	-t ${phename} \
	-c ${chr} \
                        -n ${phename} \
                        -o /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline/eur_dgn/assoc/${phename}
