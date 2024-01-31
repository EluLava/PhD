# IMPUTATION

Imputation and phasing were performed with GLIMPSE (v.1.1.1). We used GLIMPSE v.1 as recommended with small reference panels.

Scripts used for this step are in the script folder but I'll show a bunch of mock commands here for one Super-Scaffold: Super-Scaffold 2

The first step is to divide the super scaffold into chunks

`CHR_chunking --input /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.vcf.gz \
        --region Super-Scaffold_2 --window-size 2000000 --buffer-size 200000 --output ./3_GLIMPSE_chunks/Chunks_Super-Scaffold_2.txt`

The second step is the actual imputation

`while IFS="" read -r LINE || [ -n "$LINE" ];
do
  	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT=./4_GLIMPSE_phasing/Super-Scaffold_2_${ID}.bcf
        GLIMPSE_phase --input ./2_merged3KowlsVCFs/Super-Scaffold_2.vcf.gz --reference /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_Super-Scaffold_2.vcf.gz \
	 --thread 30 --ne 10000 --burnin 100 --main 15 --input-region ${IRG} --map /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/RecombinationMap/Super-Scaffold_2.map.txt \
        --samples-file /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/2_GLIMPSE/data/inputFILES/GLIMPSEploidyfiles/Super-Scaffold_2.ploidy --output-region ${ORG} --output ${OUT}
        bcftools index -f ${OUT}
done < ${CHUNKfile}`

Here ploidy file is also just 2 for all individuals, it would change only for the sexual chromosome (see inputFILES).

The third step is to ligate the imputed chunks

`GLIMPSE_ligate --input <(find ./4_GLIMPSE_phasing/ -name "Super-Scaffold_2_*.bcf") --output ./5_GLIMPSE_ligating/Super-Scaffold_2_phased_imputed.bcf`

The forth step is to extract the most likely haplotype

`GLIMPSE_sample --input ./5_GLIMPSE_ligating/Super-Scaffold_2_phased_imputed.bcf --solve --output ./6_GLIMPSE_haplotypes/Super-Scaffold_2_phased_imputed_haplotypes.bcf`
