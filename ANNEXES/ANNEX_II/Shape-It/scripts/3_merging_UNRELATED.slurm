#! /bin/bash

#SBATCH --partition cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --job-name=MergeIndVCF
#SBATCH -o STD/%j.stdout
#SBTACH -e STD/%j.stderr
#SBATCH --account jgoudet_barn_owl

### Run this script with ###
#
# awk '{print "sbatch ./2.6.1_MergeIndividualVCFs_RP315.slurm "$1}' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta.fai | bash
#
### ###

set -e
set -x

module load gcc vcftools bcftools shapeit4/4.1.3 whatshap htslib/1.15.1

#SCAFFOLD=Super-Scaffold_1
SCAFFOLD=$1

myPATH="/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/"

# Produce the list of imput files
awk -v ss=$SCAFFOLD '{print "/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPfamilymembers315/IndFiles/"$1"_"ss"_ShapeItPhased.vcf.gz"}' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/FamilyMembers_Libnames.list > STD/$SCAFFOLD.list

# use the file to merge the VCFs
bcftools merge -l STD/$SCAFFOLD.list --threads 5 -O z -o ${myPATH}//Out_ShapeIt/RPfamilymembers315/RP315_Libnames_ShapeItPhased_${SCAFFOLD}.vcf.gz

# Index the vcf
bcftools index --threads 5 ${myPATH}/Out_ShapeIt/RPfamilymembers315/RP315_Libnames_ShapeItPhased_${SCAFFOLD}.vcf.gz

#rm list of files
rm STD/$SCAFFOLD.list
