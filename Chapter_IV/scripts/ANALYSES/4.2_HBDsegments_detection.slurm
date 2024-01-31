#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 03-00:00:00
#SBATCH --array [1-53]
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 64G
#SBATCH --job-name RzooRoH3K
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --export=NONE
#SBATCH --account jgoudet_barn_owl

set -x
set -e

cd /users/elavanc1/BARNOWL/ID_in3K/data

#Load the modules
module load gcc r plink-ng/1.9.20200712 bcftools htslib

#Extract ss from SLURM TASK ID
ss=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./inputFILES/SS_AUTOSAUMES_3K.list)

#First step is to create a file with correct format for changing POS field
awk 'NR > 1 {print $1"\t"$2"\t"$3}' ./1_VCFs/SNPsLists/PhysicalANDGeneticPOSplus10_ss_${ss}.txt | bgzip > ./1_VCFs/perSuperScaffolds/TMPfromPHYStoGENpos_All3085_AUTOSAUMES_RP502SNPs_${ss}.tab.gz

#index
tabix -p vcf ./1_VCFs/perSuperScaffolds/TMPfromPHYStoGENpos_All3085_AUTOSAUMES_RP502SNPs_${ss}.tab.gz

#replace POS in VCF
bcftools annotate -a ./1_VCFs/perSuperScaffolds/TMPfromPHYStoGENpos_All3085_AUTOSAUMES_RP502SNPs_${ss}.tab.gz -c CHROM,POS,~POS -O z -o ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10.vcf.gz \
./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_${ss}.vcf.gz

bcftools index ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10.vcf.gz

#rm the TMP file
rm ./1_VCFs/perSuperScaffolds/TMPfromPHYStoGENpos_All3085_AUTOSAUMES_RP502SNPs_${ss}.tab.gz

#1. Convert from VCF to OXford format
plink --vcf ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10.vcf.gz \
--recode oxford --allow-extra-chr --out ./1_VCFs/perSuperScaffolds/TMP_OX_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10

#Open R
R --vanilla <<EOF

    library(RZooRoH)
    library(doParallel)
    library(foreach)
    cl <- 10
    registerDoParallel(cl)

    # read the OX file
    data_Rohs = zoodata("./1_VCFs/perSuperScaffolds/TMP_OX_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10.gen", samplefile = "./inputFILES/All3085_NewNamesCORRECTED.list", zformat = "gp")

    #### RUN THE MODEL

    # Create Mdodel: 13 HBD and 1 non-HBD class.es
    Mod <- zoomodel(K=14, krates=c(2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,8192), err = 3.76e-5)

    # Run the model
    loc_mod <- zoorun(Mod, data_Rohs, localhbd = TRUE, nT = cl)

    # Extract table with all detected segments
    loc_table <- loc_mod@hbdseg

    #We need to update the information about INDVs NAMES
    loc_table[,1] = data_Rohs@sample_ids[loc_table[,1]]

    #Loop through CHR to fill the columns
    for(chr in 1:length(data_Rohs@chrnames)){
	#Fill the CHROM column
	loc_table[,2][loc_table[,2] == chr] = data_Rohs@chrnames[chr]
    }

    #Save this file
    write.table(loc_table, "./4_ROHs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}.hom", quote = F, col.names = T, row.names = F)

    #save the entire SESSION
    save.image(file = "./4_ROHs/EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_${ss}.RData")

EOF

#rm TMP OX file
rm ./1_VCFs/perSuperScaffolds/TMP_OX_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}_GenPOSplus10.*
