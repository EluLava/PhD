#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=500G
#SBATCH --job-name=VCFsmerging3Kowls
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

#Launch the function
INDVs_VCF_merging ./inputFILES/REFPANELSuperScaffolds.list ./2_merged3KowlsVCFs 30

### For GLIMPSE we need to PAR and UNPAR region of SS42 separately !!! So we'll just split it there

#PAR region
bcftools view -t Super-Scaffold_42:1-4415000,Super-Scaffold_42:7364001-7389000 -O z -o ./2_merged3KowlsVCFs/Super-Scaffold_42PAR.vcf.gz ./2_merged3KowlsVCFs/Super-Scaffold_42.vcf.gz
bcftools index ./2_merged3KowlsVCFs/Super-Scaffold_42PAR.vcf.gz
#Non PAR region
bcftools view -t Super-Scaffold_42:4415001-7364000,Super-Scaffold_42:7389001-41426968 -O z -o ./2_merged3KowlsVCFs/Super-Scaffold_42nonPAR.vcf.gz ./2_merged3KowlsVCFs/Super-Scaffold_42.vcf.gz
bcftools index ./2_merged3KowlsVCFs/Super-Scaffold_42nonPAR.vcf.gz

#rm the UNSPLIT SS42
rm ./2_merged3KowlsVCFs/Super-Scaffold_42.vcf.gz*

#I'll split the REF PANEL SS42 as well !
#Par region
bcftools view -t Super-Scaffold_42:1-4415000,Super-Scaffold_42:7364001-7389000 -O z -o \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42PAR.vcf.gz \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42.vcf.gz
bcftools index /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42PAR.vcf.gz
#NonPAR region
bcftools view -t Super-Scaffold_42:4415001-7364000,Super-Scaffold_42:7389001-41426968 -O z -o \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42nonPAR.vcf.gz \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42.vcf.gz
bcftools index /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_42nonPAR.vcf.gz
