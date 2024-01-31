#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --job-name=VCFs_3KRP346_merging
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

#Load modules
module load gcc bcftools

#change directory
cd /users/elavanc1/BARNOWL/ID_in3K/data

#3K AUTOSAUMES VCF
VCF3K=/users/elavanc1/3KOWLS/2_GLIMPSE/data/8_GLIMPSES_VCFs_diploidHomozygFemales/3K2768_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz

#create list of SS from this VCF
bcftools query -f'%CHROM\n' ${VCF3K} | sort | uniq > ./inputFILES/SS_AUTOSAUMES_3K.list

#RP502 VCF: several problems: LIBNAMES, all 502 individuals and ZCHR included
RP502VCF=~/3KOWLS/4_REF_PANEL/data/Out_ShapeIt/RPall502/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_allSuperScaffolds.vcf.gz

#first step is to substract only swiss + AUTOSAUMES from vcf
bcftools view -r $(awk '{printf $0","}' ./inputFILES/SS_AUTOSAUMES_3K.list) -S ./inputFILES/RP346CH_LibNames.list -O z -o ./1_VCFs/RPCH346_Libnames_AUTOSAUMES_RP502SNPs.vcf.gz ${RP502VCF}
bcftools index ./1_VCFs/RPCH346_Libnames_AUTOSAUMES_RP502SNPs.vcf.gz

#From LibNames to NewNames
bcftools reheader -s ./inputFILES/RP346CH_LibNames_NewNames.list -o ./1_VCFs/RPCH346_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz ./1_VCFs/RPCH346_Libnames_AUTOSAUMES_RP502SNPs.vcf.gz
bcftools index ./1_VCFs/RPCH346_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz

#merge both VCFs
bcftools merge -O z -o ./1_VCFs/All3114_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz ./1_VCFs/RPCH346_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz ${VCF3K}
bcftools index ./1_VCFs/All3114_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz

#remove the two weird individuals from low cov: 889830 and M038227; plus twins, plus indviduals in^dentified as weird in RP and 3K SEPARATELY
bcftools view -s ^889830,M038227,M026346,M026347,M031634,M040522,M040721,M043108,M043255,M043259,M043262 -O z -o ./1_VCFs/All3102_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz ./1_VCFs/All3114_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz
bcftools index ./1_VCFs/All3102_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz

#rm previous VCF with weird INDVs
rm ./1_VCFs/All3114_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz*

### Then we ran estimated genomic based kinship (with allele sharing method) and compared it to the expected pedigree kinship to remove weird individuals from the dataset (twins, mislabeled individuals, etc.)
# The next lines are for removing these weird individuals:

#count number of indvs to remove
NBLINES=$(cat ./inputFILES/WeirdKinshipINDVs.list | wc -l)

#remove individuals identified as weird
bcftools view -s ^$(awk -v nr=${NBLINES} '{if(NR < nr){printf $0","}else{printf $0}}' ./inputFILES/WeirdKinshipINDVs.list) -O z -o ./1_VCFs/All3085_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz ./1_VCFs/All3102_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz
#index
bcftools index ./1_VCFs/All3085_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz

#rename the three way sample swap individuals
bcftools reheader -s ./inputFILES/swapedSamplesBasedonKinship.list -o ./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz ./1_VCFs/All3085_NewNames_AUTOSAUMES_RP502SNPs.vcf.gz
bcftools index ./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz

#create list of samples
bcftools view -h ./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz | tail -n1 | awk '{for(i=10;i<=NF;i++){print $i}}' > ./inputFILES/All3085_NewNamesCORRECTED.list
