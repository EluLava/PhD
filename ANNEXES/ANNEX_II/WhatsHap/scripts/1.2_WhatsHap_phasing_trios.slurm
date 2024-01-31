#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:00:00
#SBATCH --array 1-13370 #### HOWEVER TOO MANY ARRAYS in this so I just divided it into 3 ARRAYS 2.2.2.1 2.2.2.2 and 2.2.2.3 in 2.2_filteringonwhatsHap folder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=TrioPhasing
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

#Load the modules
module load gcc whatshap parallel htslib bcftools samtools

#Shortcut for PATH
myPATH=/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data

#create runLINE input file: ${myPATH}/inputFILES/input_whatsHap_TriosPhasing.list
#parallel --dryrun ::: $(cut -d' ' -f1 ${myPATH}/FamilyInfo/FullTrios_Libnames.ped) :::+ $(cut -d' ' -f2 ${myPATH}/FamilyInfo/FullTrios_Libnames.ped) :::+ $(cut -d' ' -f3 ${myPATH}/FamilyInfo/FullTrios_Libnames.ped) :::+ $(cut -d' ' -f4 ${myPATH}/FamilyInfo/FullTrios_Libnames.ped) ::: $(cat ${myPATH}/filteredVCFs/SuperScaffolds.list) > ${myPATH}/inputFILES/input_whatsHap_TriosPhasing.list

#runLine has input arguments
runLINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/inputFILES/input_whatsHap_TriosPhasing.list)

#define function
WhatsHap_Pedigree_Phasing(){

    #the first argument is the focal individual
    focus=$1

    #the second argument is the paternal ID
    patID=$2

    #the third argument is the maternal ID
    matID=$3

    #indvs to phase is the forth argument
    indvsToPhase=$4

    #the fifth argument is the scaffold
    ss=$5

    #The first step is to subset the three individuals we want to phase from the trio + corresponding ss
    bcftools view -r ${ss} -s ${focus},${patID},${matID} -o ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz -O z ${myPATH}/filteredVCFs/RP504_Libnames_TF1_Mask_indDP_COMPLETE.vcf.gz
    bcftools index ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz

    echo "main individual is ${focus}, paternal ID is ${patID}, maternal ID is ${matID} and PED file is ${myPATH}/FamilyInfo/tmp_${focus}.ped"
    echo "BAM --> ${myPATH}/1_FamilyBAMs/FAM_${focus}.bam"

    #phase each family and interva+subintervall individually
    whatshap phase \
    --ped ${myPATH}/FamilyInfo/tmp_${focus}.ped \
    -o ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf \
    --reference /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta \
    --indels \
    --recombrate 2 \
    ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz \
    ${myPATH}/1_FamilyBAMs/FAM_${focus}.bam

    #rm unphased vcf
    rm ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz*

    # Subsample INDIVIDUALS we want to phase + zip VCF
    bcftools view -s ${indvsToPhase} ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf -O z -o ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf.gz

    #rm the UNzip and 3 indvs file
    rm ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf

    #Index zipped and final VCF file
    bcftools index ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf.gz
}

export -f WhatsHap_Pedigree_Phasing

WhatsHap_Pedigree_Phasing ${runLINE}
