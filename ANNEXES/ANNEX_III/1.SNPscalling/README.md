# SNPs Calling

SNP calling was performed with BCFTools (v.1.15.1) on the SNPs present in the reference panel.

Scripts used for this step are in the script folder but I'll show a bunch of mock commands here for one Super-Scaffold: Super-Scaffold 2

We started by extracting the list of SNPs present in the reference panel for this super scaffold

`bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.vcf.gz \
| bgzip -c > 1_SNPcalling_3Kowls/SNPs_RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.tsv.gz`

Then we index this list of SNPs

`tabix -s1 -b2 -e2 SNPs_RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.tsv.gz`

Then we use BCFTools to call the SNPs at these positions in the 2,768 low coverage individuals. This function calls all individuals in parallel:

`parallel --rpl '{/..} s:^.*/::;s:_markdup\.bam::;' -j30 "bcftools mpileup -f /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta -I -E -a 'FORMAT/DP' -T \
        /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.vcf.gz -r Super-Scaffold_2 {} -O u | \
	bcftools call --ploidy-file ./inputFILES/ploidyFile.txt -S ./inputFILES/TMPlist_{/..}.tsv -Aim -C alleles -T \
	./1_SNPcalling_3Kowls/SNPs_RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.tsv.gz -O z -o 1_SNPcalling_3Kowls/{/..}_Super-Scaffold_2.vcf.gz \
	::: $(find "~/3KOWLS/1_TrimmAlign/data/5_renamed_BAMs/" -maxdepth 1 -name "*.bam")`

The ploidy file is simply diploid for all individuals in this case. See scripts folder for sexual chromosome

Then we merge the individual VCFs into one VCF per super scaffold

`bcftools merge --threads 5 -m none -r Super-Scaffold_2 -l <(ls ./1_SNPcalling_3Kowls/*_Super-Scaffold_2.vcf.gz) -O z -o /Super-Scaffold_2.vcf.gz`
