#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=GLIMPSEchunking
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

#First step is to re-create the list of ss for GLIMPSE (becasue we splitted SS42 into PAR and unPAR regions)
find ./2_merged3KowlsVCFs -name "*.vcf.gz" | xargs -I {} bash -c 'basename -s ".vcf.gz" {}' > ./inputFILES/GLIMPSESuperScaffolds.list

#Then simply launch the function
CHR_chunking ./inputFILES/GLIMPSESuperScaffolds.list /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds ./3_GLIMPSE_chunks 2
