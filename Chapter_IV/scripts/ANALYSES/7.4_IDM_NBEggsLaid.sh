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
    names(birds)[34] = "YearBorn"
    #Extract year
    birds[,34] = gsub('-[0-9][0-9]-[0-9][0-9]','', birds[,34])

    #Read clutch dataset
    clutch = read.csv("./5_ID_Models/5.1_data/BarnOwls_Legacy_2023.10.23/BarnOwls_Legacy_20231023111917_Clutch.csv")
    #Subset FemaleRing and number of eggs laid as well as year
    clutch = clutch[,c(2,3,12,17,20)]
    #Update column names for merging with dtaF
    colnames(clutch)[c(1,4)] = c("YearClutch", "RingId")
    #remove clutch where we don't know the mom
    clutch[clutch[,4] == "",4] = NA
    #remove clutches for which mom = NA
    clutch = clutch[!is.na(clutch[,4]),]
    #remove clutches for which Size is unknown
    clutch = clutch[!is.na(clutch[,5]),]
    #change Year clutch to year format
    clutch[,1] = gsub(' [0-9][0-9]:[0-9][0-9]:[0-9][0-9]','', clutch[,3])
    clutch[,1] = as.Date(clutch[,1], format = "%d/%m/%Y")
    #Extract year
    clutch[,1] = gsub('-[0-9][0-9]-[0-9][0-9]','', clutch[,1])

    #Read color data
    dtaColor = read.csv("./5_ID_Models/5.1_data/BarnOwls_Legacy_2023.10.23/BarnOwls_Legacy_20231023111917_BirdColour.csv")
    #Extract year of measurment
    dtaColor[,29] = gsub(' [0-9][0-9]:[0-9][0-9]:[0-9][0-9]','', dtaColor[,5])
    dtaColor[,29] = as.Date(dtaColor[,29], format = "%d/%m/%Y")
    names(dtaColor)[29] = "YearObs"
    #Extract year
    dtaColor[,29] = gsub('-[0-9][0-9]-[0-9][0-9]','', dtaColor[,29])

    #Subset only the individuals we observed in the clutch data and for the years they had clutches (because we don't want nestlings to infleuce color scores)
    dtaColor = dtaColor[paste0(dtaColor[,2], "_", dtaColor[,29]) %in% paste0(clutch[,4], "_", clutch[,1]),]

    #Get mean color per individual (mean of Breast and Belly per Obs, then mean per INDV)
    dtaColorMEAN = aggregate((dtaColor[,11] + dtaColor[,12])/2, by = list(dtaColor[,2]), FUN = mean)
    #Update col names
    colnames(dtaColorMEAN) = c("RingId", "Color")

    #Get mean spotiness Volume per individual (sum of Breast and Belly volume (nb * (diam/2)^2 * pi) per Obs, then mean per INDV)
    dtaSpotinessVolumeMEAN = aggregate((dtaColor[,25] * (((dtaColor[,26]/2)^2)*pi)) + (dtaColor[,27] * (((dtaColor[,28]/2)^2)*pi)), by = list(dtaColor[,2]), FUN = mean)
    #Update col names
    colnames(dtaSpotinessVolumeMEAN) = c("RingId", "SpotinessVolume")

    ###############################################################################
    ############################### NUMBER OF EGGS  ###############################
    ###############################################################################


    ############################### ANIMAL MODEL ###############################

    #update pedigree format for modeling
    p2 <- as.data.frame(ped[,c(1,2,3)])
    ord <- pedigree::orderPed(p2)
    p2 <- p2[order(ord),] #now you can use p2 to fit your animal model

    #Sub-sample data to remove NA (tarsus length) and only keep adults!
    Tarsuslengthdta <- dtaMorpho[ (!is.na(dtaMorpho[,16])) &
				(dtaMorpho[,13] == 'Adult'),
                 		c(2,16)]

    #Mean TarsusLength value per INDV
    TarsuslengthdtaMEAN = aggregate(Tarsuslengthdta[,2], by = list(Tarsuslengthdta[,1]), FUN = mean)
    #change colnames
    colnames(TarsuslengthdtaMEAN) = c("RingId", "MeanLeftTarsus")

    #Add clutch ID RAISED and Year born
    TarsuslengthdtaMEAN = merge(TarsuslengthdtaMEAN, unique(birds[,c(2,33,34)]), by = "RingId", all.x = T)

    #merge with clutch data
    Clutchdta = merge(clutch, TarsuslengthdtaMEAN, by = "RingId")

    #add color
    Clutchdta = merge(Clutchdta, dtaColorMEAN, by = "RingId")

    #add SpotinessVolume
    Clutchdta = merge(Clutchdta, dtaSpotinessVolumeMEAN, by = "RingId")

    #Laying date as Date
    Clutchdta[,4] = as.Date(Clutchdta[,4], format = "%d/%m/%Y")

    #Extract age of the parent for this year
    Clutchdta[,11] = as.numeric(as.character(Clutchdta[,2])) - as.numeric(as.character(Clutchdta[,8]))
    colnames(Clutchdta)[11] = "FemaleAge"

    #Remove observations for whihc female age is lower than 1 (only 1)
    Clutchdta = Clutchdta[Clutchdta[11] > 0,]

    #add one variable with date in julian Day (nth day of the year)
    Clutchdta[,12] = (Clutchdta[,4] - as.Date(paste0(Clutchdta[,2], "-01-01")))
    colnames(Clutchdta)[12] = "JulDay"

    #merge with dtaF
    Clutchdta = merge(Clutchdta, dtaF, by = "RingId")

    #change RingID to animal
    names(Clutchdta)[1] <- 'id'

    #add sex
    Clutchdta[,15] <- ped[match(Clutchdta[,1], ped[,1]),4]
    colnames(Clutchdta)[15] = "sex"
    #sex 3 -> NA
    Clutchdta[,15][Clutchdta[,15] == 3] = NA
    #sex as factor
    Clutchdta[,15] = factor(Clutchdta[,15], levels = c(1,2))

    #Year of clutch as factor
    Clutchdta[,2] <- as.factor(Clutchdta[,2])
    #SiteId as factor
    Clutchdta[,3] <- as.factor(Clutchdta[,3])

    # non - na covariates
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,5]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,6]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,9]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,10]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,11]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,12]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,14]),]

    #set factor levels with only observations we have
    Clutchdta[,1] = factor(Clutchdta[,1], levels = unique(as.character(Clutchdta[,1]))[!is.na(unique(as.character(Clutchdta[,1])))])
    Clutchdta[,2] = factor(Clutchdta[,2], levels = unique(as.character(Clutchdta[,2]))[!is.na(unique(as.character(Clutchdta[,2])))])
    Clutchdta[,3] = factor(Clutchdta[,3], levels = unique(as.character(Clutchdta[,3]))[!is.na(unique(as.character(Clutchdta[,3])))])
    Clutchdta[,7] = factor(Clutchdta[,7], levels = unique(as.character(Clutchdta[,7]))[!is.na(unique(as.character(Clutchdta[,7])))])
    Clutchdta[,8] = factor(Clutchdta[,8], levels = unique(as.character(Clutchdta[,8]))[!is.na(unique(as.character(Clutchdta[,8])))])
    Clutchdta[,13] = factor(Clutchdta[,13], levels = unique(as.character(Clutchdta[,13]))[!is.na(unique(as.character(Clutchdta[,13])))])

    # as data frame
    Clutchdta <- as.data.frame(Clutchdta)
    #JulDay as numeric
    Clutchdta[,12] = as.numeric(Clutchdta[,12])
    #remove weird observation where JulDay is smaller than 30
    Clutchdta = Clutchdta[Clutchdta[,12] > 30,]

    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, nu = 0.002), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.02, alpha.V = 1000),
                           G2 = list(V = 1, nu = 0.02, alpha.V = 1000),
			   G3 = list(V = 1, nu = 0.02, alpha.V = 1000),
			   G4 = list(V = 1, nu = 0.02, alpha.V = 1000)))

    #Model itself
    NBEggslaidIDM <- MCMCglmm(SizeOriginally ~ FuniWE + MeanLeftTarsus + JulDay + FemaleAge,
                    random = ~ id + YearClutch + SiteId + ClutchIDBorn_Raised,
                    family='gaussian',
                    data = Clutchdta,
                    prior = prior,
                    nitt = 1e6,
                    burnin = 50000,
                    thin=500)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/IDM_NBEggsLaid_ADULTS1.RData")

EOF
