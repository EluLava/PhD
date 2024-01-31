#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=256G
#SBATCH --job-name=FuniWE_est
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

module load gcc r

cd /users/elavanc1/BARNOWL/ID_in3K/data

mkdir -p ./2_Fcoeff/

R --vanilla << EOF

    library(gaston)

    #source the get.funiwT function
    source("../scripts/functions/functions.R")

    #read VCF to BED matrix
    bed = read.vcf("./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz", convert.chr = FALSE)

    #estimate funiWE
    FuniWE = get.funiwT(bed = bed, nb.cores = 30)

    #Add INDIVIDUALS names
    FuniWEdata = as.data.frame(cbind(names(FuniWE),FuniWE))
    colnames(FuniWEdata) = c("INDVs","FuniWE")
    FuniWEdata[,2] = as.numeric(FuniWEdata[,2])

    #save output
    write.table(FuniWEdata, "./2_Fcoeff/FuniWE.txt", quote = F, col.names = T, row.names = F)

EOF
