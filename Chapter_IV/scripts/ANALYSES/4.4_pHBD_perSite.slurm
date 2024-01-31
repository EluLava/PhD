#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 11:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 500G
#SBATCH --job-name HBDprobperSNPandHBDclass3K
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --export=NONE
#SBATCH --account jgoudet_barn_owl

#Load the modules
module load gcc r

cd /users/elavanc1/BARNOWL/ID_in3K/data

R --vanilla <<EOF

    library("RZooRoH")
    library("data.table")

    #read sample IDs
    sampleIDs = as.vector(read.table("./inputFILES/All3085_NewNamesCORRECTED.list")[,1])

    #read ss list
    SuperScaffolds = as.vector(read.table("./inputFILES/SS_AUTOSAUMES_3K.list")[,1])

    #Read NB SNPs per ss
    SNPsNB = as.vector(read.table("./inputFILES/SNPsnb_per_ss.list", h = F)[,2])

    #Create empty dataframe that we'll fill with SuperScaffolds, ncol = 16: nb classes (13 HBD 1 non-HBD) + 2 (CHR + POS)
    HBDprob = as.data.frame(matrix(ncol = 16, nrow = 0))
    colnames(HBDprob) = c("CHROM","POS",paste0("class",seq(1,14)))

    #keep a number for SNPsNB index
    i = 0

    #Read the super scaffolds R ENV files
    for(ss in SuperScaffolds){

        print(paste0("Starting ", ss))

	i = i + 1

        #Load the R file
        load(list.files(path="./4_ROHs", pattern=paste0("EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_",ss,".RData"), full.names = T))

        #Create empty dataframe that we'll fill with eah individuals
	HBDprobint = as.data.frame(matrix(ncol = 16, nrow = SNPsNB[i]))
	colnames(HBDprobint) = colnames(HBDprob)

	#Fill chrom column
	HBDprobint[,1] = ss
	#Fill position column
	HBDprobint[,2] = data_Rohs@bp

        #Loop through HBD classes
        for(class in 1:14){

	    #create empty dataframe that we'll fill with all INDIVIDUALS so we can average among individuals later on
	    HBDprobclass = as.data.frame(matrix(ncol = 3085, nrow = SNPsNB[i]))

	    #Loop throgh individuals to extract their prob for this specific class
	    for(indv in 1:3085){

		#Extract values from model
		HBDprobclass[,indv] = loc_mod@hbdp[[indv]][class,]

	    }

	    #Get the mean per per SNP and fill the interval df
	    HBDprobint[,class+2] = apply(HBDprobclass, 1, mean)

        }

        #merge the interval HBD prob with the overall HBD prob
        HBDprob = rbind(HBDprob, HBDprobint)

    }

    #Write the final df
    write.table(HBDprob, "./4_ROHs/HBD_prob_per_SNP_perHBDclass.txt", sep = "\t", row.names = F, col.names = T, quote = F)

    #get the sum of pHBD from classes 1 to 10
    pHBD = fread("./4_ROHs/HBD_prob_per_SNP_perHBDclass.txt", header = T)
    pHBD = as.data.frame(pHBD)

    #rm ss with less than 1000 SNPs
    largeSS = names(table(pHBD[,1])[table(pHBD[,1]) > 10000])
    pHBD = pHBD[pHBD[,1] %in% largeSS,]

    #get new dataframe with CHR, POS, pHBD classes 1 to 10
    pHBDuptoClass10 = as.data.frame(matrix(nrow = nrow(pHBD), ncol = 3))
    colnames(pHBDuptoClass10) = c("CHROM", "GenPOS", "pHBD")

    #fill CHROM aand POS columns
    pHBDuptoClass10[,1] = pHBD[,1]
    pHBDuptoClass10[,2] = pHBD[,2]

    #get the pHBD
    pHBDuptoClass10[,3] = rowSums(pHBD[,3:12])

    #merge with physical position
    dtaGenPhyPos = read.table("./1_VCFs/SNPsLists/PhysicalANDGeneticPOSplus10_ALL_Super-Scaffolds.list", h = T)
    colnames(dtaGenPhyPos)[3] = "GenPOS"

    pHBDuptoClass10.2 = merge(pHBDuptoClass10, dtaGenPhyPos, by = c("CHROM","GenPOS"))

    #sort the output
    pHBDuptoClass10.3 = pHBDuptoClass10.2[order(pHBDuptoClass10.2[,1], pHBDuptoClass10.2[,4]), ]

    #Add a new column with PLOTTING POS
    pHBDuptoClass10.3 = as.data.frame(cbind(pHBDuptoClass10.3, PLOTPOS = vector(mode = "numeric", length = nrow(pHBDuptoClass10.3))))
    #Add a new column with plotting colot (factor 1 or 2)
    pHBDuptoClass10.3 = as.data.frame(cbind(pHBDuptoClass10.3, COLPOS = vector(mode = "numeric", length = nrow(pHBDuptoClass10.3))))

    #Fill this new plotting position column
    pHBDuptoClass10.3[1,5] = pHBDuptoClass10.3[1,4]
    pHBDuptoClass10.3[1,6] = 1

    #loop through rows
    for(row in 2:nrow(pHBDuptoClass10.3)){

	#if same ss
	if(pHBDuptoClass10.3[row,1] == pHBDuptoClass10.3[(row - 1),1]){
	    pHBDuptoClass10.3[row,5] = (pHBDuptoClass10.3[row,4] - pHBDuptoClass10.3[(row - 1),4]) + pHBDuptoClass10.3[(row - 1), 5]
	    pHBDuptoClass10.3[row,6] = pHBDuptoClass10.3[(row - 1),6]
	} else {
	    pHBDuptoClass10.3[row,5] = pHBDuptoClass10.3[(row - 1),5] + 100
	    pHBDuptoClass10.3[row,6] = pHBDuptoClass10.3[(row - 1),6] + 1
	}
    }

    #savee the df
    write.table(pHBDuptoClass10.3, "./4_ROHs/HBD_prob_per_SNP_HBDClasses1to10_PLOTPOS.txt", quote = F, col.names = T, row.names = F)

EOF

