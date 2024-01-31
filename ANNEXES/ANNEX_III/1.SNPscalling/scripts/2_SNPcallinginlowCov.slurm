#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=4G
#SBATCH --job-name=SNPsCalling3Kowls
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

cd /users/elavanc1/3KOWLS/2_GLIMPSE/data

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source /users/elavanc1/3KOWLS/2_GLIMPSE/scripts/functions_GLIMPSE_1.0_3Kowls.sh

#BAMpath
BAMpath=~/3KOWLS/1_TrimmAlign/data/5_renamed_BAMs/

#Ref genome fasta
fastaREF=/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta

REF_SNPs_Calling_perINDV ${BAMpath} ./inputFILES/REFPANELSuperScaffolds.list ${fastaREF} ./1_SNPcalling_3Kowls 30
