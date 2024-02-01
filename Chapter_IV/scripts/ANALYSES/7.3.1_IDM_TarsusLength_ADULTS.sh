set -x
set -e

#Load modules
module load gcc r

#change directory
cd /users/elavanc1/BARNOWL/ID_in3K/data

R --vanilla << EOF

    #load libraries
    library(MCMCglmm)
    library(pedigree)

    #Read the F data
    dtaF = read.table("./2_Fcoeff/FuniWE.txt", h = T)
    #list of RingIDs
    RingIDs = dtaF[,1]
    #Update colnames
    colnames(dtaF) = c("RingId","FuniWE")
    #Read genetic sex dta (from coverage WGS)
    dtaGenSex = read.table("./5_ID_Models/5.1_data/GeneticSex_3K_RP548.txt", h = T)

    #merge F and sex data
    dtaF = merge(dtaGenSex, dtaF, by = "RingId")

    #Read the morpho data
    dtaMorpho = read.csv("./5_ID_Models/5.1_data/BarnOwls_Legacy_2023.10.23/BarnOwls_Legacy_20231023111917_BirdMeasurement.csv")
    #rm Genetic sex column
    dtaMorpho = dtaMorpho[,-12]

    #Read the pedigree
    ped <- read.table("./5_ID_Models/5.1_data/pedigreeCORRECTED.tab", h = T)

    #Read the bird file
    birds = read.csv("./5_ID_Models/5.1_data/BarnOwls_Legacy_2023.10.23/BarnOwls_Legacy_20231023111917_Bird.csv")
    #add a column merging born and raised clutch IDs
    birds[,33] = paste0(as.character(birds[,10]), "_", as.character(birds[,12]))
    birds[,33][birds[,33] == "NA_NA"] = NA
    names(birds)[33] = "ClutchIDBorn_Raised"
    birds[,33] = as.factor(birds[,33])
    #new column with correct format for HatchYear
    birds[,34] = gsub(' [0-9][0-9]:[0-9][0-9]:[0-9][0-9]','', birds[,6])
    birds[,34] = as.Date(birds[,34], format = "%d/%m/%Y")
    names(birds)[34] = "HatchDate.2"
    #Some RingIDs are encoded as "", pass these to NA
    birds[birds[,2] == "",2] = NA
    #remove lines with NA as RingId
    birds = birds[!is.na(birds[,2]),]

    #Extract year
    birds[,35] = gsub('-[0-9][0-9]-[0-9][0-9]','', birds[,34])
    names(birds)[35] = "YearBorn"

    #Create empty df to fill
    RanKindv = as.data.frame(cbind(RingId = vector(mode = "character", length = 0), Rank = vector(mode = "numeric", length = 0)))

    #For some reasons some clutchID are not included in the clutch dataset so I'll just use INDV ID from birds RaisedClutchID (because the clutchID from Morpho data is shit)

    #loop through individuals
    for(indv in unique(birds[,2])){

        #check if individual has already been ranked with one of his siblings
        if(indv %in% RanKindv[,1]){

            next

        } else {

            #extract clutch ID
            clutchid = unique(birds[birds[,2] == indv,12][!is.na(birds[birds[,2] == indv,12])])

            #if we don't know clutch ID forget it
            if(length(clutchid) == 0){

                linestoadd = as.data.frame(cbind(RingId = indv, Rank = NA))
                #Add this to RanKindv dataframe
                RanKindv = rbind(RanKindv, linestoadd)

            } else {

                #extract individuals within the same clutch ID
                otherINDVs = unique(birds[(birds[,12] == clutchid) & (!is.na(birds[,12])),2])

                #Rank the individuals
                linestoadd = as.data.frame(cbind(RingId = birds[birds[,2] %in% otherINDVs, 2],
                                         Rank = rank(birds[birds[,2] %in% otherINDVs, 34])))

                #Add this to RanKindv dataframe
                RanKindv = rbind(RanKindv, linestoadd)

            }
	}
    }

    #remove duplicated lines
    RanKindv = unique(RanKindv)

    #############################################################################
    ############################### TARSUS LENGTH ###############################
    #############################################################################


    ############################### ANIMAL MODEL ###############################

    #update pedigree format for modeling
    p2 <- as.data.frame(ped[,c(1,2,3)])
    ord <- pedigree::orderPed(p2)
    p2 <- p2[order(ord),] #now you can use p2 to fit your animal model

    #Sub-sample data to remove NA (wing length) and only keep adults!
    Tarsuslengthdta <- dtaMorpho[ (!is.na(dtaMorpho[,16])) &
				(dtaMorpho[,13] == 'Adult'),
                 		c(2,9,8,16)]

    #some oberver are "" --> NA
    Tarsuslengthdta[,2][Tarsuslengthdta[,2] == ""] = NA

    Tarsuslengthdta = unique(Tarsuslengthdta)

    #Add clutch ID RAISED and YearBorn
    Tarsuslengthdta = merge(Tarsuslengthdta, unique(birds[,c(2,33,35)]), by = "RingId")

    #add rank in the dataframe
    Tarsuslengthdta = merge(Tarsuslengthdta, RanKindv, by = "RingId")
    Tarsuslengthdta[,7] = as.numeric(as.character(Tarsuslengthdta[,7]))

    #add sex
    Tarsuslengthdta[,8] <- ped[match(Tarsuslengthdta[,1], ped[,1]),4]
    colnames(Tarsuslengthdta)[8] = "sex"
    #sex 3 -> NA
    Tarsuslengthdta[,8][Tarsuslengthdta[,8] == 3] = NA
    #sex as factor
    Tarsuslengthdta[,8] = factor(Tarsuslengthdta[,8], levels = c(1,2))

    # non - na covariates
    Tarsuslengthdta <- Tarsuslengthdta[!is.na(Tarsuslengthdta[,4]),]
    Tarsuslengthdta <- Tarsuslengthdta[!is.na(Tarsuslengthdta[,7]),]

    #subsample only individuals sequenced
    TarsuslengthdtaSEQ = Tarsuslengthdta[Tarsuslengthdta[,1] %in% RingIDs,]
    #merge with InbreedingCoeff
    TarsuslengthdtaSEQF = merge(TarsuslengthdtaSEQ, dtaF, by = "RingId")

    #GeneticSex as factor
    TarsuslengthdtaSEQF[,9] = factor(TarsuslengthdtaSEQF[,9], levels = c(1,2))

    #change INDV colnames
    colnames(TarsuslengthdtaSEQF)[1] = "id"

    #set factor levels with only observations we have
    TarsuslengthdtaSEQF[,1] = factor(TarsuslengthdtaSEQF[,1], levels = unique(as.character(TarsuslengthdtaSEQF[,1]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,1])))])
    TarsuslengthdtaSEQF[,2] = factor(TarsuslengthdtaSEQF[,2], levels = unique(as.character(TarsuslengthdtaSEQF[,2]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,2])))])
    TarsuslengthdtaSEQF[,5] = factor(TarsuslengthdtaSEQF[,5], levels = unique(as.character(TarsuslengthdtaSEQF[,5]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,5])))])
    TarsuslengthdtaSEQF[,6] = factor(TarsuslengthdtaSEQF[,6], levels = unique(as.character(TarsuslengthdtaSEQF[,6]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,6])))])
    TarsuslengthdtaSEQF[,8] = factor(TarsuslengthdtaSEQF[,8], levels = unique(as.character(TarsuslengthdtaSEQF[,8]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,8])))])
    TarsuslengthdtaSEQF[,9] = factor(TarsuslengthdtaSEQF[,9], levels = unique(as.character(TarsuslengthdtaSEQF[,9]))[!is.na(unique(as.character(TarsuslengthdtaSEQF[,9])))])

    # as data frame
    TarsuslengthdtaSEQF <- as.data.frame(TarsuslengthdtaSEQF)

    #There is one individual with a very low (TO LOW) tarsusLength and is far from the other measures of the same individual so we'll remove it !
    TarsuslengthdtaSEQF = TarsuslengthdtaSEQF[TarsuslengthdtaSEQF[,4] > 600,]


    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, nu = 0.002), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.02, alpha.V = 1000),
			   G2 = list(V = 1, nu = 0.1, alpha.V = 1000),
			   G3 = list(V = 1, nu = 0.5, alpha.V = 1000),
			   G4 = list(V = 1, nu = 0.5, alpha.V = 1000)))

    #Model itself
    TarsuslengthIDM <- MCMCglmm(LeftTarsus ~ GeneticSex + Rank + FuniWE,
                    random = ~ id + Observer + ClutchIDBorn_Raised + YearBorn,
                    family='gaussian',
                    data = TarsuslengthdtaSEQF,
                    prior = prior,
                    nitt = 1e6,
                    burnin = 50000,
                    thin=500)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/IDM_TarsusLength_ADULTSRank.RData")

EOF
