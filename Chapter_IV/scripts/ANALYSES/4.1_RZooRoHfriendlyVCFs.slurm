#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 10:00:00
#SBATCH --array 1-53
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --job-name RecMapIntervalplus10
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --export=NONE

set -e
set -x

cd /users/elavanc1/BARNOWL/ID_in3K/data

mkdir -p 4_ROHs

#Load the modules
module load gcc r bcftools

#Extract SuperScaffold
ss=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./inputFILES/SS_AUTOSAUMES_3K.list)

#Extract ss  number only
ssnb=$(echo ${ss} | cut -d'_' -f2)

#Subsample the SNPs for this ss
bcftools view -r ${ss} -O z -o ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_${ss}.vcf.gz ./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz
#index
bcftools index ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_${ss}.vcf.gz

#Extract the list of SNPs we have for this ss
bcftools query -f '%POS\n' ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_${ss}.vcf.gz > ./1_VCFs/SNPsLists/SNPsLISTPHYSICALPOS_ss_${ss}.list

#check if the recombination map exists
if [[ -f /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/genmaps/ss${ssnb}_genetic_map.txt ]]; then

    #If the file exists we'll just extrapolate with Tristan script (thanks <3)
    R --vanilla << EOF

	#Read the rec map
	RecMap = read.table("/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/genmaps/ss${ssnb}_genetic_map.txt", header = T)[,c(1,3)]
	#Update colnames so it matches what the functions expects
	colnames(RecMap) = c("BP", "CM")

	#Read the list of SNPs we have in the VCF
	SNPslist = read.table("./1_VCFs/SNPsLists/SNPsLISTPHYSICALPOS_ss_${ss}.list", header = F)
	colnames(SNPslist) = "BP"

	#Load the function
	source("../scripts/functions/functions.R")

	#Create the new map
	NewMap = interpolate_cm(RecMap, SNPslist)
	#Multiply cM position by 500,000 (because this is what RZoo takes as input, conversion of cM to (roughly) bp and in our case (barn owl) 1cM = roughly 500,000 bp)
	NewMap[,2] = NewMap[,2]*500000
	colnames(NewMap)[2] = "cMx500000"
	#round the positions in cM
        NewMap[,2] = round(NewMap[,2], digits = 0)
        NewMap[,2] = as.integer(NewMap[,2])

	#Add the ss as column
	NewMap.2 = as.data.frame(cbind(CHROM=rep("${ss}", nrow(NewMap)), NewMap))

	## Then create the same map but with increase of one (which roughly correpsond to one bp) in cMx.6e6 if we have many positions with same values !
	#Loop through rows of the dataframe
	for(row in 2:nrow(NewMap.2)){

	    #extract previous position
	    previouspos = NewMap.2[(row-1),3]
	    #extract current pos
	    pos = NewMap.2[row,3]

	    #extract diff
	    diff = abs(pos - previouspos)

	    #if same
	    if(diff < 10){

		#all the next positions will be added a plus one !
		NewMap.2[row:nrow(NewMap.2),3] = NewMap.2[row:nrow(NewMap.2),3] + 10

	    }
	}

	NewMap.2[,3] = round(NewMap.2[,3], digits = 0)
	NewMap.2[,3] = as.integer(NewMap.2[,3])

	#write the new MAP
	write.table(NewMap.2, "./1_VCFs/SNPsLists/PhysicalANDGeneticPOSplus10_ss_${ss}.txt", quote = F, col.names = T, row.names = F)


EOF

else

    R --vanilla << EOF

	#Read list of SNPs
	SNPslist = try(read.table("./1_VCFs/SNPsLists/SNPsLISTPHYSICALPOS_ss_${ss}.list", header = F))

	if(inherits(SNPslist, "try-error")){

	    print(paste0("No SNPs to proceeed for ", ${ss}))

	} else {

	    colnames(SNPslist) = "BP"

	    #Set value to go from physic positions to genetic distances because not the same for autosaumes and sexual chromosomes
	    if("${ss}" %in% c("Super-Scaffold_13","Super-Scaffold_42")){
		val = 750000
	    } else {
		val = 500000
	    }

	    #Create CONSTANT map
	    NewMap = as.data.frame(cbind(CHROM=rep("${ss}", nrow(SNPslist)), BP=as.numeric(SNPslist[,1]), CM=as.numeric(SNPslist[,1]/val)))
	    #Pass second and third column to numbers
	    NewMap[,2] = as.integer(NewMap[,2])
	    NewMap[,3] = as.numeric(NewMap[,3])
	    #Multiply cM position by 500,000 (because this is what RZoo takes as input, conversion of cM to (roughly) bp and in our case (barn owl) 1cM = roughly 500,000 bp, EXCEPT sexual CHR)
            NewMap[,3] = NewMap[,3]*val
	    colnames(NewMap)[3] = "cMx500000"
	    #round the positions in cM
	    NewMap[,3] = round(NewMap[,3], digits = 0)
	    NewMap[,3] = as.integer(NewMap[,3])

	    ## Then create the same map but with increase of one (which roughly correpsond to one bp) in cMx.5e6 if we have many positions with same values !
	    #Loop through rows of the dataframe
	    for(row in 2:nrow(NewMap)){

		#extract previous position
    		previouspos = NewMap[(row-1),3]
    		#extract current pos
    		pos = NewMap[row,3]

		#extract diff
		diff = abs(pos - previouspos)
		#if same
		if(diff < 10){

		    #all the next positions will be added a plus one !
		    NewMap[row:nrow(NewMap),3] = NewMap[row:nrow(NewMap),3] + 10
		}
	    }

	    NewMap[,3] = round(NewMap[,3], digits = 0)
	    NewMap[,3] = as.integer(NewMap[,3])

	    #write the new MAP
	    write.table(NewMap, "./1_VCFs/SNPsLists/PhysicalANDGeneticPOSplus10_ss_${ss}.txt", quote = F, col.names = T, row.names = F)

	}

EOF

fi
