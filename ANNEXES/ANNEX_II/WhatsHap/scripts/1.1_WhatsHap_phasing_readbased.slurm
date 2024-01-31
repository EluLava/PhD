#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=10:00:00
#SBATCH --array 1-35280 ##### HOWEVER THIS IS TOO MUCH FOR THE CLUSTER so I SPLITTED THIS SCRIPT INTO 8 SUB-SCRIPTS named 2.2.1, 2.2.2, 2.2.3, etc. in 2.2_filteringonwhatsHap folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --job-name=ReadBasedPhasing
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --account jgoudet_barn_owl

#Load the modules
module load gcc whatshap parallel htslib bcftools samtools

#Shortcut for PATH
myPATH=/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data

#create input file
#parallel --dryrun ::: $(cat ${myPATH}/FamilyInfo/RP504_Libnames.list) ::: $(cat ${myPATH}/filteredVCFs/SuperScaffolds.list) > ${myPATH}/inputFILES/input_whatsHap_ReadBasedPhasing.list

#input argumetnd according to SLURM TASK ID
runLINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/inputFILES/input_whatsHap_ReadBasedPhasing.list)

#define function
WhatsHap_readbased_Phasing(){

    #the first argument is the focal individual
    focus=$1

    #the second argument is the scaffold
    ss=$2

    #The first step is to subset only the individuals we want to phase + ss
    bcftools view -r ${ss} -s ${focus} -o ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz -O z ${myPATH}/filteredVCFs/RP504_Libnames_TF1_Mask_indDP_COMPLETE.vcf.gz
    bcftools index ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz

    echo "main individual is ${focus}, BAM is in ${myPATH}/BAM_BQSRed/${focus}_bqsr.bam"

    #phase indiviual ss with BAM reads
    whatshap phase \
    -o ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf \
    --reference /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta \
    --indels \
    --recombrate 2 \
    ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz \
    ${myPATH}/BAM_BQSRed/${focus}_bqsr.bam

    #rm unphased vcf
    rm ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz*

    # Zip VCF file
    bcftools view ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf -O z -o ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf.gz

    #rm the UNzip file
    rm ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf

    #Index zipped VCF file
    bcftools index ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf.gz

}

export -f WhatsHap_readbased_Phasing

WhatsHap_readbased_Phasing ${runLINE}
