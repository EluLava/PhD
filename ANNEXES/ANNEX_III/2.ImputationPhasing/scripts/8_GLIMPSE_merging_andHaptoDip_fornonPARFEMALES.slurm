#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --job-name=Hap2DipVcfmerging
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

cd /users/elavanc1/3KOWLS/2_GLIMPSE/data

#Load the modules
module load gcc bcftools parallel htslib curl boost

#Pass all the BCFs into zipped VCF
parallel --rpl '{/..} s:^.*/::;s:\.bcf::;' -j 20 'bcftools view -O z -o ./7_GLIMPSED_VCFs_HaploidFemales/perSuper-Scaffold/{/..}.vcf.gz {}' ::: $(find ./6_GLIMPSE_haplotypes -name "*.bcf")

#index the VCFs
parallel -j 20 'bcftools index {}' ::: $(find ./7_GLIMPSED_VCFs_HaploidFemales/perSuper-Scaffold/ -name "*.vcf.gz")

## merge all VCFs into one big VCF and index

#List VCFs
find ./7_GLIMPSED_VCFs_HaploidFemales/perSuper-Scaffold/ -name "*.vcf.gz" | sort -V > ./inputFILES/TMP_OutputGLIMPSEvcfs.list
#ulimit to avoid problem of opening mroe files than allowed at the time
ulimit -n 10000
#Merge
bcftools concat -f ./inputFILES/TMP_OutputGLIMPSEvcfs.list -O z --threads 20 -a | bcftools sort -m 126G -O z -o ./7_GLIMPSED_VCFs_HaploidFemales/3K2768_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz
#Index
bcftools index ./7_GLIMPSED_VCFs_HaploidFemales/3K2768_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz
#rm list file
rm ./inputFILES/TMP_OutputGLIMPSEvcfs.list

## Pass females haploid to homozygous diploid for ss13 and ss42nonPAR

mkdir -p 8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold

#get all VCFs into non-zipped format in the new folder
parallel --rpl '{/..} s:^.*/::;s:\.vcf.gz::;' -j 20 'bcftools view -O v -o ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/{/..}.vcf {}' ::: $(find ./7_GLIMPSED_VCFs_HaploidFemales/perSuper-Scaffold/ -name "*.vcf.gz")

## sed the SS13 and SS42nonPAR

#ss13
#sed -r -e 's:/([[:space:]])([0,1])(^[|])/\1\2|\2\3/g' ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes.vcf > ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes_allDIP.vcf
#FOR SOME REASON DOES NOT WORK --> try awk
#Get header
bcftools view -h ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes.vcf > ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes_allDIP.vcf

#Grep the non header, pipe that into awk, loop trhough columns
bcftools view -H ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes.vcf | awk '{for(i=1;i<=NF;i++){

	#if the first 9 columns (INFO FORMAT) etc. just print the column with a tab
	if(i < 10){

	    #print the column with a tab
	    printf $i "\t"

	#else if the GT is 3 digits (diploid)
	} else if(length($i) == 3){

	    #print GT
	    printf $i "\t"

	#else if the GT field is only one digit (haploid) and not the last column
	} else if((length($i) == 1) && (i < NF)){

	    #print GT | GT and a tab
	    printf $i "|" $i "\t"

	#else if the GT is one digit (haploid) AND the last column
	} else if((length($i) == 1) && (i == NF)){

            #print GT | GT
            print $i "|" $i

	}
    }
}' >> ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes_allDIP.vcf

#change name because for some reason I cannot run sed -ire and I am tired (cannot directly modify the file with awk so all same thing here)
mv ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes_allDIP.vcf ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes.vcf
#change FPLOIDY in header
sed -i 's/##FPLOIDY=-2/##FPLOIDY=2/' ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_13_phased_imputed_haplotypes.vcf

#ss42nonPAR
bcftools view -h ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes.vcf > ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes_allDIP.vcf

#Grep the non header, pipe that into awk, loop trhough columns
bcftools view -H ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes.vcf | awk '{for(i=1;i<=NF;i++){

        #if the first 9 columns (INFO FORMAT) etc. just print the column with a tab
        if(i < 10){

            #print the column with a tab
            printf $i "\t"

        #else if the GT is 3 digits (diploid)
        } else if(length($i) == 3){

            #print GT
            printf $i "\t"

        #else if the GT field is only one digit (haploid) and not the last column
        } else if((length($i) == 1) && (i < NF)){

            #print GT | GT and a tab
            printf $i "|" $i "\t"

        #else if the GT is one digit (haploid) AND the last column
        } else if((length($i) == 1) && (i == NF)){

            #print GT | GT
            print $i "|" $i

        }
    }
}' >> ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes_allDIP.vcf

#change name because for some reason I cannot run sed -ire and I am tired (cannot directly modify the file with awk so all same thing here)
mv ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes_allDIP.vcf ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes.vcf
#change FPLOIDY in header
sed -i 's/##FPLOIDY=-2/##FPLOIDY=2/' ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/Super-Scaffold_42nonPAR_phased_imputed_haplotypes.vcf

#zip all VCFs
parallel --rpl '{/..} s:^.*/::;s:\.vcf::;' -j 30 'bcftools view -O z -o ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/{/..}.vcf.gz {}' ::: $(find ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/ -name "*.vcf")

#index
parallel  -j 30 'bcftools index {}' ::: $(find ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/ -name "*.vcf.gz")

#rm unzipped VCFs
rm ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/*.vcf

#merge all into one big VCF.gz
bcftools concat -f <(find ./8_GLIMPSES_VCFs_diploidHomozygFemales/perSuper-Scaffold/ -name "*.vcf.gz" | sort -V) -O z --threads 20 -a | bcftools sort -m 126G -O z -o ./8_GLIMPSES_VCFs_diploidHomozygFemales/3K2768_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz
#indexing
bcftools index ./8_GLIMPSES_VCFs_diploidHomozygFemales/3K2768_NewNames_ALL_SCAFFOLDS_RP502SNPs.vcf.gz
