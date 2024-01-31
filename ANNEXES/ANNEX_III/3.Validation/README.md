# GLIMPSE PIEPLINE VALIDATION

We used two different steps to validate the GLIMPSE pipeline:

## The DUPLICATES

The duplicates consist of 22 individuals sequenced once at high coverage and once at low coverage. The aim is to compare the same individuals' genotypes between hogh coverage classical sequencing and the GLIMPSE pipeline. In the main run of GLIMPSE they were excluded from the low coverage individuals and inclduded in the reference panel !

For the DUPLICATES run, these individuals were included in the low coverage set (resulting in 2,790 low coverage individuals) and excluded from the reference panel (resulting in a reference panel containing 480 individuals).

The scripts are the same as the main GLIMPSE run except for the list of individuals and the names of the final VCFs (RPnoDUP480_Libnames_Filtered_WhatshapPhased_ShapeitPhased.vcf.gz for the reference panel and 3K2790withDUP_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz for the low coverage). The list of duplicates can be found in the 3K metadata and in the 3Kduplicates.list file (new names) and in the RefPanelduplicates.list (library names) in the inputFILES folder.

## The TRIPLICATES

The triplicates consist of 10 individuals sequenced once at high coverage and three times at low coverage (in three different plates sent for sequencing). The aim is to increase the sample size for the comparison of individuals' genotypes between high coverage classical seuqencing and low coverage GLIMPSE piepline but also to account for an eventual effect of sequencing between the different plates. In the main GLIMPSE pipeline, these individuals were excluded from the low coverage set and included in the reference panel.

The triplicates GLIMPSE pipeline was run three times (a,b and c) because we cannot have the same individuals several times in one run. In the a run all the triplicates names are NAMEa, in the b run all the triplicates names are NAMEb, in the c run all the triplicates names are NAMEc.

The scripts are the same as the main GLIMPSE run except for the list of individuals and the names of the final VCFs (RPnoTRIP492_Libnames_Filtered_WhatshapPhased_ShapeitPhased.vcf.gz for the reference panel and 3K2778withTRIPa_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz for the low coverage). The list of triplicates can be found in the 3K metadata and in the 3Ktriplicates.[a,b,c].list (new names) and in the RefPaneltriplicates.list (library names) in the inputFILES folder.
