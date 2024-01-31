#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:00:00
#SBATCH --array 1-70
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --job-name=merging_phased_IND_VCFs
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

#Load the modules
module load gcc parallel bcftools

#save shortcut for path
myPATH=/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data

#command line to create the list of individuals phased within trios (they do not each have on VCF but will be used later to create the list of individuals to get from read-based only phasing)
#cut -d' ' -f4 /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/FullTrios_Libnames.ped | awk -F "," '{for(i=1;i<=NF;i++){print $i}}' > ./FamilyInfo/INdividuals_phased_within_trios.list

#command line to generate the list of individuals phased with reads only (basically lines whcih are in RP504_Libnames.list but not in INdividuals_phased_within_trios.list)
#comm -13 <(sort /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/INdividuals_phased_within_trios.list) <(sort /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/RP504_Libnames.list) > /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/INdividuals_phased_with_reads.list

#Extract input from SLURM ID TASK, each array is one super-scaffold
runLINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p <(cat ${myPATH}/filteredVCFs/SuperScaffolds.list))

#define function which will be launched per super-scaffold
mergingPhasedVCFs(){

    #First argument passed to the function is the super-scaffold
    ss=$1

    #If sexual CHR, just take read-based phasing, else take trio phased for indvs phased in trio !
    if [[ ${ss} == "Super-Scaffold_13" || ${ss} == "Super-Scaffold_42" ]]; then

    #The list of VCFs to extract from read-based phasing is all INDVs
    find ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs -name "Indv_*_ss_${ss}_readbased.vcf.gz"  > ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list

    #sort the file so that all individuals are always in the same order
    sort ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list > ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list

    rm ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list

    #We then simply used the file we created to merge the VCFs
    bcftools merge -l ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list --threads 5 -O z -o ${myPATH}/mergedWhatsHap_phased/RP504_Libnames_TF1_Mask_indDP_whatsHap_phased_ss_${ss}.vcf.gz

    #rm list of files
    rm ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list

    #We then index the file
    bcftools index --threads 5 ${myPATH}/mergedWhatsHap_phased/RP504_Libnames_TF1_Mask_indDP_whatsHap_phased_ss_${ss}.vcf.gz

    else

    #the list of vcf to extract from trios phased VCF is all the VCFs in the folder (easy)
    find ${myPATH}/WhatsHap_phased_TRIOS/trios_VCFs/ -name "Indv_*_ss_${ss}_TrioPhased.vcf.gz" > ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list

    #The list of VCFs to extract from read-based phasing is all INDVs present in INdividuals_phased_with_reads.list
    parallel -j1 find ${myPATH}/WhatsHap_phased_READBASE/individual_VCFs -name "Indv_{}_ss_${ss}_readbased.vcf.gz" ::: $(cat ${myPATH}/FamilyInfo/INdividuals_phased_with_reads.list) >> ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list

    #sort the file so that all individuals are always in the same order
    sort ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list > ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list

    rm ${myPATH}/inputFILES/files_to_merge_after_whatsHap_ss_${ss}.list

    #We then simply used the file we created to merge the VCFs
    bcftools merge -l ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list --threads 5 -O z -o ${myPATH}/mergedWhatsHap_phased/RP504_Libnames_TF1_Mask_indDP_whatsHap_phased_ss_${ss}.vcf.gz

    #rm list of files
    rm ${myPATH}/inputFILES/SORTED_files_to_merge_after_whatsHap_ss_${ss}.list

    #We then index the file
    bcftools index --threads 5 ${myPATH}/mergedWhatsHap_phased/RP504_Libnames_TF1_Mask_indDP_whatsHap_phased_ss_${ss}.vcf.gz

   fi
}

export -f mergingPhasedVCFs

mergingPhasedVCFs ${runLINE}
