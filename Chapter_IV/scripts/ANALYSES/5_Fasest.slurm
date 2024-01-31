#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 11:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 256G
#SBATCH --job-name MatchingMatandFas3K
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --export=NONE
#SBATCH --account jgoudet_barn_owl

#Load the modules
module load gcc r

cd /users/elavanc1/BARNOWL/ID_in3K/data

R --vanilla <<EOF

    library(hierfstat)
    library(gaston)

    #source for fas funnction
    source("../scripts/functions/functions.R")

    #read the vcf
    bed = read.VCF("./1_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.vcf.gz", convert.chr = F)

    #getMatching Matrix full dataset
    MatchingMat = get.matchingT(bed, nb.cores = 2)

    #Extract self-Allele matching
    Mii = diag(MatchingMat)

    #get beta kinships (coancestries) from the matching matrix
    BetaKin = beta.dosage(MatchingMat, inb = F, MATCHING = T)

    #get rid of the diag
    diag(BetaKin) = NA

    #only keel the upper triangle of the matrix
    BetaKin[upper.tri(BetaKin, diag = F)] = NA

    #select the list of UNrelated samples
    UNR = find.unrelated(BetaKin, 0.05)

    #get the list of unrelated individuals
    UNR.list = rownames(UNR)

    #save lit of unrelated
    write.table(UNR.list, "./inputFILES/UNRELATED_indvs.list", quote = F, row.names = F, col.names = F)

    #subsample only the unrelated individuals
    bedUNR = select.inds(bed, id %in% UNR.list)
    #remove monomorphic SNPs
    bedUNR.snps = select.snps(bedUNR, maf > 0)

    #get matching matrix UNRELATED set
    MatchingMatUNR = get.matchingT(bedUNR.snps, nb.cores = 2)

    #Mean matching UNR (= MB in beta.dosage function)
    MB = mean(MatchingMatUNR[upper.tri(MatchingMatUNR, diag = F)])

    #get Kinship with Mij full dataset but MB unrelated samples
    KinshipUNRallfrq = (MatchingMat - MB)/(1 - MB)

    #save the BetaKinship matrix
    saveRDS(object = KinshipUNRallfrq, file = "./2_Fcoeff/AS_UNRallfrq_GRM.RDS")

    #Get Fas for all INDVs (Mii from full MatchingMat) but MB ony unrelated
    Fas = ((Mii * 2 - 1) - MB)/(1 - MB)

    #merge with INDVs names
    dtaF = as.data.frame(cbind(INDVs = rownames(MatchingMat), Fas = Fas))
    dtaF[,2] = as.numeric(dtaF[,2])

    #save dtaF dataframe
    write.table(dtaF, "./2_Fcoeff/Fas_UNRallfrq_3085.txt", quote = F, col.names = T, row.names = F)

EOF
