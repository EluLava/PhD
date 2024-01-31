#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=MergeSC_RP502
#SBATCH -o STD/%j.stdout
#SBTACH -e STD/%j.stderr
#SBATCH --account jgoudet_barn_owl

module load gcc bcftools shapeit4/4.1.3

### To Run This Script ###
#
# awk '{print "sbatch ./3.1_MergeRefPanels.slurm "$1}' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta.fai| bash
#
### ###

# Uses scaffold name as a parameter of this script
scaffold=$1
echo ${scaffold}

#Shortcut for PATH
myPATH="/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/"
inputVCF_Unrel="${myPATH}/Out_ShapeIt/RPunrelated187/RPunrelated187_Libnames_TF1_Mask_indDP_whatsHap_phased_biall_mac5_missing95_ss_"${scaffold}"_ShapeitPhased.vcf.gz"
inputVCF_Fam="${myPATH}/Out_ShapeIt/RPfamilymembers315/RP315_Libnames_ShapeItPhased_"${scaffold}".vcf.gz"

OutputVCF="${myPATH}/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_"${scaffold}".vcf.gz"

#Define the function
ref='/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta' # reference path

bcftools merge ${inputVCF_Unrel} ${inputVCF_Fam} -O z -o ${OutputVCF}
bcftools index ${OutputVCF}

#bcftools view ${outputVCF}_${scaffold}_ShapeitPhased.vcf -O z -o ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
#bcftools index ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
