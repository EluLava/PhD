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

    #only keep unique info
    RanKindv = unique(RanKindv)


    ###########################################################################
    ############################### BILL LENGTH ###############################
    ###########################################################################


    ############################### ANIMAL MODEL ###############################

    #update pedigree format for modeling
    p2 <- as.data.frame(ped[,c(1,2,3)])
    ord <- pedigree::orderPed(p2)
    p2 <- p2[order(ord),] #now you can use p2 to fit your animal model

    #Sub-sample data to remove NA (bill length) and only keep adults!
    Billlengthdta <- dtaMorpho[ (!is.na(dtaMorpho[,28])) &
				(dtaMorpho[,13] == 'Adult'),
                 		c(2,9,8,28)]

    #remove potentially duplicated visits
    Billlengthdta = unique(Billlengthdta)

    #Add clutch ID RAISED and YearBorn
    Billlengthdta = merge(Billlengthdta, unique(birds[,c(2,33,35)]), by = "RingId")

    #add rank in the dataframe
    Billlengthdta = merge(Billlengthdta, RanKindv, by = "RingId")
    Billlengthdta[,7] = as.numeric(as.character(Billlengthdta[,7]))

    #add sex
    Billlengthdta[,8] <- ped[match(Billlengthdta[,1], ped[,1]),4]
    colnames(Billlengthdta)[8] = "sex"
    #sex 3 -> NA
    Billlengthdta[,8][Billlengthdta[,8] == 3] = NA
    #sex as factor
    Billlengthdta[,8] = factor(Billlengthdta[,8], levels = c(1,2))

    #chnage RingId name
    colnames(Billlengthdta)[1] = "id"

    #add animal variable
    Billlengthdta[,9] = Billlengthdta[,1]
    colnames(Billlengthdta)[9] = "animal"

    # non - na covariates
    Billlengthdta <- Billlengthdta[!is.na(Billlengthdta[,4]),]
    Billlengthdta <- Billlengthdta[!is.na(Billlengthdta[,8]),]

    #set factor levels with only observations we have
    Billlengthdta[,1] = factor(Billlengthdta[,1], levels = unique(as.character(Billlengthdta[,1]))[!is.na(unique(as.character(Billlengthdta[,1])))])
    Billlengthdta[,2] = factor(Billlengthdta[,2], levels = unique(as.character(Billlengthdta[,2]))[!is.na(unique(as.character(Billlengthdta[,2])))])
    Billlengthdta[,5] = factor(Billlengthdta[,5], levels = unique(as.character(Billlengthdta[,5]))[!is.na(unique(as.character(Billlengthdta[,5])))])
    Billlengthdta[,6] = factor(Billlengthdta[,6], levels = unique(as.character(Billlengthdta[,6]))[!is.na(unique(as.character(Billlengthdta[,6])))])

    # as data frame
    Billlengthdta <- as.data.frame(Billlengthdta)

    #subsample only the individuals present in the pedigree
    Billlengthdta = Billlengthdta[Billlengthdta[,9] %in% p2[,1],]

    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, nu = 0.002), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.002, alpha.V = 1000),
			   G2 = list(V = 1, nu = 0.05, alpha.V = 1000), # we are more sure that INDV is flat !
			   G3 = list(V = 1, nu = 0.02, alpha.V = 1000), # same Obs mostly flat
			   G4 = list(V = 1, nu = 0.02, alpha.V = 1000), # same here mostly flat
			   G5 = list(V = 1, nu = 0.02, alpha.V = 1000))) # same here mostly flat

    #Model itself
    BilllengthIDMnoRank <- MCMCglmm(BillLength ~ sex,
                    random = ~ animal + id + Observer + ClutchIDBorn_Raised + YearBorn,
                    family='gaussian',
		            pedigree = p2,
                    data = Billlengthdta,
                    prior = prior,
                    nitt = 1e6,
                    burnin = 50000,
                    thin=500)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/AM_BillLength_ADULTS_noRank.RData")

EOF