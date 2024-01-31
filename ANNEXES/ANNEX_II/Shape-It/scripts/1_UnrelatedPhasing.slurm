#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=12:00:00
# #SBATCH --array 1-70
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --job-name=ShapeIt_Unrelated
#SBATCH -o STD/%j.stdout
#SBTACH -e STD/%j.stderr
#SBATCH --account jgoudet_barn_owl

module load gcc bcftools shapeit4/4.1.3

### To Run This Script ###
#
# awk '{print "sbatch ./2.5.1_ShapeIt_RefPanelUnrelated.slurm "$1}' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta.fai| bash
#
### ###

# Uses scaffold name as a parameter of this script
scaffold=$1
echo ${scaffold}

#Shortcut for PATH
myPATH=/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/
inputVCF="${myPATH}/mergedWhatsHap_phased_filtered/RPunrelated187_Libnames_TF1_Mask_indDP_whatsHap_phased_biall_mac5_missing95_ss_"${scaffold}".vcf.gz"
outputVCF="${myPATH}/Out_ShapeIt/RPunrelated187/RPunrelated187_Libnames_TF1_Mask_indDP_whatsHap_phased_biall_mac5_missing95_ss"

#Define the function
ref='/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta' # reference path

map="${myPATH}/RecombinationMap/"${scaffold}".map.txt"

if [ -f "$map" ]; then
	echo "I found a Map Captain"
	shapeit4 --input ${inputVCF} \
	         -T 24 \
	         -R ${scaffold} --sequencing \
		 --map ${myPATH}/RecombinationMap/${scaffold}.map.txt \
	         --use-PS 0.0001 \
	         --pbwt-depth 8 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
	         --output ${outputVCF}_${scaffold}_ShapeitPhased.vcf
else
	echo "I dont have a Map"
	shapeit4 --input ${inputVCF} \
	         -T 24 \
	         -R ${scaffold} --sequencing \
	         --use-PS 0.0001 \
	         --pbwt-depth 8 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
	         --output ${outputVCF}_${scaffold}_ShapeitPhased.vcf
fi

bcftools view ${outputVCF}_${scaffold}_ShapeitPhased.vcf -O z -o ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
bcftools index ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
