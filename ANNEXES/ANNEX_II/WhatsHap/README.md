# WhatsHap PHASING

- scripts contains the raw codes

- inputFILES contains important information such as list of trios and whihc individual was phased with which trio

## READ-BASED PHASING

All individulas were phased with reads (only) with WhatsHap (v.1.4). This was done per individual (focus) and super-scaffold (ss)

``
#${focus} bash variable is the individual LIB name
#${ss} bash variable is the super-scaffold name
``

``
#The first step was to subsample the individual and super-scaffold
bcftools view -r ${ss} -s ${focus} -o /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz -O z /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/filteredVCFs/RP504_Libnames_TF1_Mask_indDP_COMPLETE.vcf.gz
#Index the VCF
bcftools index /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz
``

``whatshap phase \
-o /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_READBASE/individual_VCFs/Indv_${focus}_ss_${ss}_readbased.vcf \
--reference /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta \
--indels \
--recombrate 2 \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data//WhatsHap_phased_READBASE/individual_VCFs/${focus}_ss_${ss}.vcf.gz \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data//BAM_BQSRed/${focus}_bqsr.bam
``


## TRIOS-BASED PHASING

All individuals belonging to one (or more) trio·s were phased with the pedigree option of WhatsHap. Parents which we re present in more than one trios (which have more than one kid) were phased along with the kid who had the higher coverage.

``
#${focus} bash variable is the individual LIB name
#${patID} bash variable is the paternal ID
#${matID} bash variable is the maternal ID
#${indvsToPhase} bash is the list of individual·s to phase (the kid always + for some trios one or two parents) 
#${ss} bash variable is the super-scaffold name
``

``
#The first step was to subsample the individual·s and super-scaffold
bcftools view -r ${ss} -s ${indvsToPhase} -o /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz -O z /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/filteredVCFs/RP504_Libnames_TF1_Mask_indDP_COMPLETE.vcf.gz
bcftools index /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz
``

``
whatshap phase \
--ped /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/FamilyInfo/tmp_${focus}.ped \
-o /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_TRIOS/trios_VCFs/Indv_${focus}_ss_${ss}_TrioPhased.vcf \
--reference /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta \
--indels \
--recombrate 2 \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/WhatsHap_phased_TRIOS/trios_VCFs/${focus}_ss_${ss}.vcf.gz \
/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/1_FamilyBAMs/FAM_${focus}.bam
``
