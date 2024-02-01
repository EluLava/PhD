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

    # HatchDate in Date Format
    birds[,34] = gsub(' [0-9][0-9]:[0-9][0-9]:[0-9][0-9]','', birds[,6])
    birds[,34] = as.Date(birds[,34], format = "%d/%m/%Y")
    names(birds)[34] = "HatchDate.2"

    #Add YearBorn
    birds[,35] = gsub('-[0-9][0-9]-[0-9][0-9]','', birds[,34])
    names(birds)[35] = "YearBorn"

    #Some RingIDs are encoded as "", pass these to NA
    birds[birds[,2] == "",2] = NA
    #remove lines with NA as RingId
    birds = birds[!is.na(birds[,2]),]

    # read clutch data for day of last chick which hatched
    clutch = read.csv("./5_ID_Models/5.1_data/BarnOwls_Legacy_2023.10.23/BarnOwls_Legacy_20231023111917_Clutch.csv")

    ## get hatch date per clutch

    # create empty df to fill
    LastHatchDayperClutchID = as.data.frame(matrix(ncol = 2, nrow = 0))
    colnames(LastHatchDayperClutchID) = c("clutchID", "lastHatchday")

    #Get last day of hatching per clutch
    for(clutchID in unique(clutch[,1])) {

        #extract last hatching day
        lastHatchday = max(birds[,34][(birds[,11] == clutchID) & (!is.na(birds[,11]))])
        #tmp dta to merge
        tmpdat = as.data.frame(matrix(c(clutchID, as.character(lastHatchday)), ncol = 2, nrow = 1, byrow = T))
        #return clutch ID and last Hatchday
        LastHatchDayperClutchID = rbind(LastHatchDayperClutchID, tmpdat)
    }

    colnames(LastHatchDayperClutchID) = c("clutchID", "lastHatchday")

    #back to date
    LastHatchDayperClutchID[,2] = as.Date(LastHatchDayperClutchID[,2])

    #New column in dtaMorpho for whether the last chick hatched or not
    dtaMorpho = as.data.frame(cbind(dtaMorpho, LastChickHatched = vector(mode = "logical", length = nrow(dtaMorpho))))
    ## New column in dtaMorpho for date format for observation date
    dtaMorpho = as.data.frame(cbind(dtaMorpho, ObservationDate.2 = vector(mode = "character", length = nrow(dtaMorpho))))
    #remove time
    dtaMorpho[,49] = gsub(' [0-9][0-9]:[0-9][0-9]:[0-9][0-9]','', dtaMorpho[,8])
    #asDate
    dtaMorpho[,49] = as.Date(dtaMorpho[,49], format = "%d/%m/%Y")


    #Loop through observations in the Morpho data and look if the last chick had already hatched !
    for(obs in 1:nrow(dtaMorpho)){

        #Extract clutch ID for this obs
        clutchID = dtaMorpho[obs,6]

        #if we don't know the clutch ID, just NA
        if(is.na(clutchID)){dtaMorpho[obs,48] = NA

        #else if the clutch ID is not included in our LastHatchDayperCliutchID dta, also NA
        } else if(!(clutchID %in% LastHatchDayperClutchID[,1])){dtaMorpho[obs,48] = NA

        #else if the last hatch date is NA, also NA
        } else if(is.na(LastHatchDayperClutchID[LastHatchDayperClutchID[,1] == clutchID,2])){dtaMorpho[obs,48] = NA

        #else, if the observation date is before the last hatching date, FALSE
        } else if(dtaMorpho[obs,49] < LastHatchDayperClutchID[,2][LastHatchDayperClutchID[,1] == clutchID]){
            dtaMorpho[obs,48] = FALSE

        #else if the observation date is after the last chick hatched, TRUE
        } else if(dtaMorpho[obs,49] > LastHatchDayperClutchID[,2][LastHatchDayperClutchID[,1] == clutchID]){
            dtaMorpho[obs,48] = TRUE
        }
    }

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

    #remove duplicated info
    RanKindv = unique(RanKindv)

    ####################################################################
    ############################### MASS ###############################
    ####################################################################


    ############################### ANIMAL MODEL ###############################

    #update pedigree format for modeling
    p2 <- as.data.frame(ped[,c(1,2,3)])
    ord <- pedigree::orderPed(p2)
    p2 <- p2[order(ord),] #now you can use p2 to fit your animal model

    #Sub-sample data to remove NA (wing length) and only keep adults!
    Massdta <- dtaMorpho[ (!is.na(dtaMorpho[,19])) & (!is.na(dtaMorpho[,48])) &
                                (dtaMorpho[,13] == 'Adult'),
                                c(2,9,8,19,48)]

    #Some Observer encoded as "" --> NA
    Massdta[,2][Massdta[,2] == ""] = NA

    #remove duplicated rows
    Massdta = unique(Massdta)

    #Add clutch ID RAISED
    Massdta = merge(Massdta, unique(birds[,c(2,33,35)]), by = "RingId")

    #add rank in the dataframe
    Massdta = merge(Massdta, RanKindv, by = "RingId")
    Massdta[,8] = as.numeric(as.character(Massdta[,8]))

    #add sex
    Massdta[,9] <- ped[match(Massdta[,1], ped[,1]),4]
    colnames(Massdta)[9] = "sex"
    #sex 3 -> NA
    Massdta[,9][Massdta[,9] == 3] = NA
    #sex as factor
    Massdta[,9] = factor(Massdta[,9], levels = c(1,2))

    # non - na covariates
    Massdta <- Massdta[!is.na(Massdta[,4]),]
    Massdta <- Massdta[!is.na(Massdta[,5]),]
    Massdta[,5] = as.factor(Massdta[,5])
    #Massdta <- Massdta[!is.na(Massdta[,8]),]

    #subsample individuals we sequenced
    MassdtaSEQ = Massdta[Massdta[,1] %in% RingIDs,]
    #merge with dtaF
    MassdtaSEQF = merge(MassdtaSEQ, dtaF, by = "RingId")

    #GeneticSex as factor
    MassdtaSEQF[,10] = factor(MassdtaSEQF[,10], levels = c(1,2))

    #change colnames of first covariates
    colnames(MassdtaSEQF)[1] = "id"

    #set factor levels with only observations we have
    MassdtaSEQF[,1] = factor(MassdtaSEQF[,1], levels = unique(as.character(MassdtaSEQF[,1]))[!is.na(unique(as.character(MassdtaSEQF[,1])))])
    MassdtaSEQF[,2] = factor(MassdtaSEQF[,2], levels = unique(as.character(MassdtaSEQF[,2]))[!is.na(unique(as.character(MassdtaSEQF[,2])))])
    MassdtaSEQF[,6] = factor(MassdtaSEQF[,6], levels = unique(as.character(MassdtaSEQF[,6]))[!is.na(unique(as.character(MassdtaSEQF[,6])))])
    MassdtaSEQF[,7] = factor(MassdtaSEQF[,7], levels = unique(as.character(MassdtaSEQF[,7]))[!is.na(unique(as.character(MassdtaSEQF[,7])))])

    # as data frame
    MassdtaSEQF <- as.data.frame(MassdtaSEQF)

    #remove the observation with REALLY big mass 4000 (IMPOSSIBLE); this individuals have other (possible) measurments so it iss fine we'll just remove that observation)
    MassdtaSEQF = MassdtaSEQF[MassdtaSEQF[,4] < 1000,]

    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, nu = 0.002), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.02, alpha.V = 1000),
			   G2 = list(V = 1, nu = 0.002, alpha.V = 1000),
			   G3 = list(V = 1, nu = 0.1, alpha.V = 1000),
			   G4 = list(V = 1, nu = 0.1, alpha.V = 1000)))

    #Model itself
    MassIDMRank <- MCMCglmm(CalculatedMass ~ GeneticSex * LastChickHatched + FuniWE,
                    random = ~ id + Observer + ClutchIDBorn_Raised + YearBorn,
                    family='gaussian',
                    data = MassdtaSEQF,
                    prior = prior,
                    nitt = 1e6,
                    burnin = 50000,
                    thin=500)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/IDM_Mass_ADULTS_noRank.RData")

EOF
