#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 11:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 512G
#SBATCH --job-name Fest3K
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

#Load the modules
module load gcc r plink-ng/1.9.20200712 bcftools

cd /users/elavanc1/BARNOWL/ID_in3K/data

#check if the file exists, if it does not, create if

if [[ ! -f ./inputFILES/SNPsnb_per_ss.list ]]; then

#Get nb of SNPs per interval (for F weighted mean among CHRs)
for ss in $(cat ./inputFILES/SS_AUTOSAUMES_3K.list); do echo ${ss} $(bcftools view -H ./1_VCFs/perSuperScaffolds/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_${ss}.vcf.gz | wc -l) >> ./inputFILES/SNPsnb_per_ss.list; done

fi

#Open R
R --vanilla <<EOF

    library("RZooRoH")

    #read sample IDs Libnames
    sampleIDs = as.vector(read.table("./inputFILES/All3085_NewNamesCORRECTED.list")[,1])

    #read ss list
    SuperScaffolds = as.vector(read.table("./inputFILES/SS_AUTOSAUMES_3K.list")[,1])

    #### F ####

    #Create DF with just sample IDs that we'll complete with each F interval
    F_HBD_512genint = as.data.frame(sampleIDs)

    #Read the 53 ss R ENV files and extract F
    for(ss in SuperScaffolds){

        print(paste0("Starting SS ", ss))

        #Load the R file
        load(list.files(path="./4_ROHs", pattern=paste0("EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_",ss,".RData"), full.names = T))

	#Estimate Fgen 512
	Fgen512int = cumhbd(zres=loc_mod, T = 1024)

	#merge with sample IDs
	Fgen512 = as.data.frame(cbind(loc_mod@sampleids,Fgen512int))
	colnames(Fgen512)[1] = "sampleIDs"

	#merge Fgen512int with the full dataframe
	F_HBD_512genint = merge(F_HBD_512genint, Fgen512, by = "sampleIDs")
	#Set new column name
	colnames(F_HBD_512genint)[ncol(F_HBD_512genint)] = paste0("FHBD_ss_",ss)

    }

    #Read the file with nb of SNPs per interval (for F weigthed mean)
    SNPsNB = as.vector(read.table("./inputFILES/SNPsnb_per_ss.list", h = F)[,2])

    #pass FHBD columns of F_HBD_512genint into umeric
    for(col in 2:ncol(F_HBD_512genint)){F_HBD_512genint[,col] = as.numeric(as.character(F_HBD_512genint[,col]))}

    #Create new vectors with weighted means (per intervals)
    F_HBD = as.data.frame(cbind(INDVs=F_HBD_512genint[,1], FHBD512gen=as.numeric(apply(F_HBD_512genint[,2:ncol(F_HBD_512genint)],1 ,function(x){mean(x, weigths = SNPsNB)}))))

    #pass FHBD as numeric
    F_HBD[,2] = as.numeric(as.character(F_HBD[,2]))

    colnames(F_HBD) = c("INDVs","FHBD512gen")

    #read the other Fs estimated previously
    Funiwe = read.table("./2_Fcoeff/FuniWE.txt", header = T)

    #merge
    F_estimates = merge(F_HBD, Funiwe, by = "INDVs")

    #save the FHBD "final" table
    write.table(F_estimates, "./2_Fcoeff/All3085_FuniWE_FHBD512g.txt", sep = "\t", quote = F, col.names = T, row.names = F)

EOF

#For the ROHs distFile, I just need to paste all .hom files one after the other
#header first
head -n1 ./4_ROHs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_Super-Scaffold_1.hom > ./4_ROHs/HBDsegments_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.hom
#then the files
for ss in $(cat ./inputFILES/SS_AUTOSAUMES_3K.list); do awk 'NR > 1 {print $0}' ./4_ROHs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_${ss}.hom >> ./4_ROHs/HBDsegments_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.hom; done
