CRIC pipeline directory: /proj/yunligrp/users/bmlin/CRIC_DGN_pipeline

Sub-directories (for each race and trained model): 
afr_dgn: CRIC AFR samples, DGN trained models
afr_mesa: CRIC AFR samples, MESA trained models
eur_dgn: CRIC EUR samples, DGN trained models
eur_mesa: CRIC EUR samples, MESA trained models
pheno: phenotype files for each ancestry and outcome

In each of the race/model directories:
code:
	00inter.sh: filter genotypes for the variants included in the trained model
	01-chunking.sh: split each chromosome into chunks
	03-prediction.sh: predict GREX 
	04-run_assoc.sh: script for association