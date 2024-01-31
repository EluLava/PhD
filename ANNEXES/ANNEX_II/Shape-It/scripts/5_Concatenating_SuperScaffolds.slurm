#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=Concat_PR502
#SBATCH -o STD/%j.stdout
#SBTACH -e STD/%j.stderr
#SBATCH --account jgoudet_barn_owl

module load gcc bcftools shapeit4/4.1.3

### To Run This Script ###
#
# sbatch ./3.2_ConcatSuperScaffolds.slurm
#
### ###

# Prepare the file with the list of files to concatenate, in order

rm STD/ListSuperScafoldsFilesAll502.txt
touch STD/ListSuperScafoldsFilesAll502.txt
awk '{print "/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_"$1".vcf.gz"}' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta.fai | awk '{print "ls -l "$1}' | bash | awk '{print $9}' > STD/ListSuperScafoldsFilesAll502.txt

myPATH="/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/"
OutputVCF="${myPATH}/Out_ShapeIt/RPall502/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_allSuperScaffolds.vcf.gz"

#Define the function
ref='/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta' # reference path

bcftools concat -f STD/ListSuperScafoldsFilesAll502.txt -O z -o ${OutputVCF}
bcftools index ${OutputVCF}

#bcftools view ${outputVCF}_${scaffold}_ShapeitPhased.vcf -O z -o ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
#bcftools index ${outputVCF}_${scaffold}_ShapeitPhased.vcf.gz
