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

    #subset only individuals sequenced
    BilllengthdtaSEQ = Billlengthdta[Billlengthdta[,1] %in% RingIDs,]

    #merge with dtaF
    BilllengthdtaSEQF = merge(BilllengthdtaSEQ, dtaF, by = "RingId")

    #chnage RingId name
    colnames(BilllengthdtaSEQF)[1] = "id"

    # non - na covariates
    BilllengthdtaSEQF <- BilllengthdtaSEQF[!is.na(BilllengthdtaSEQF[,4]),]
    #BilllengthdtaSEQF <- BilllengthdtaSEQF[!is.na(BilllengthdtaSEQF[,7]),]
    BilllengthdtaSEQF <- BilllengthdtaSEQF[!is.na(BilllengthdtaSEQF[,9]),]

    #set factor levels with only observations we have
    BilllengthdtaSEQF[,1] = factor(BilllengthdtaSEQF[,1], levels = unique(as.character(BilllengthdtaSEQF[,1]))[!is.na(unique(as.character(BilllengthdtaSEQF[,1])))])
    BilllengthdtaSEQF[,2] = factor(BilllengthdtaSEQF[,2], levels = unique(as.character(BilllengthdtaSEQF[,2]))[!is.na(unique(as.character(BilllengthdtaSEQF[,2])))])
    BilllengthdtaSEQF[,5] = factor(BilllengthdtaSEQF[,5], levels = unique(as.character(BilllengthdtaSEQF[,5]))[!is.na(unique(as.character(BilllengthdtaSEQF[,5])))])
    BilllengthdtaSEQF[,6] = factor(BilllengthdtaSEQF[,6], levels = unique(as.character(BilllengthdtaSEQF[,6]))[!is.na(unique(as.character(BilllengthdtaSEQF[,6])))])
    BilllengthdtaSEQF[,9] = factor(BilllengthdtaSEQF[,9], levels = unique(as.character(BilllengthdtaSEQF[,9]))[!is.na(unique(as.character(BilllengthdtaSEQF[,9])))])

    # as data frame
    BilllengthdtaSEQF <- as.data.frame(BilllengthdtaSEQF)

    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, nu = 0.002), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.02, alpha.V = 1000), # we are more sure that INDV is flat !
			   G2 = list(V = 1, nu = 0.02, alpha.V = 1000), # same Obs mostly flat
			   G3 = list(V = 1, nu = 0.02, alpha.V = 1000), # same here mostly flat
			   G4 = list(V = 1, nu = 0.02, alpha.V = 1000))) # same here mostly flat

    #Model itself
    BilllengthIDMnoRank <- MCMCglmm(BillLength ~ GeneticSex + FuniWE,
                    random = ~ id + Observer + ClutchIDBorn_Raised + YearBorn,
                    family='gaussian',
                    data = BilllengthdtaSEQF,
                    prior = prior,
                    nitt = 1e6,
                    burnin = 50000,
                    thin=500)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/IDM_BillLength_ADULTS_noRank.RData")

EOF
