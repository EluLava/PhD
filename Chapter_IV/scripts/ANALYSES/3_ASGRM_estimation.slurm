#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --job-name=GRMAS_est
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

module load gcc r

cd /users/elavanc1/BARNOWL/ID_in3K/data

mkdir -p ./3_ASGRM

R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(SNPRelate)

    #source the get.matchingT function
    source("../scripts/functions/functions.R")

    #read VCF as bed matrix for matching function
    bed = read.VCF("./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz", convert.chr = FALSE)

    #Matching first because we can parallelize !
    matchingMat = get.matchingT(bed, nb.cores = 20)

    #from matching to pairwise kinships
    AsMat = beta.dosage(matchingMat, inb = FALSE, MATCHING = TRUE)

    #from kinship to GRM
    grmas = kinship2grm(AsMat)

    #write GRMas
    saveRDS(grmas, "./3_ASGRM/All3085_AUTOSAUMES_RP502SNPs.RDS")

EOF
