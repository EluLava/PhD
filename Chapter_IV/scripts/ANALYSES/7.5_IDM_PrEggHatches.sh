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
    #removed cross-fostered clutch because it sucks BECAUSE THEY DON'T KNOW WHICH EGG THEY SWAPPED
    clutch = clutch[clutch[,36] != 1,]
    clutch = clutch[clutch[,37] != 1,]
    #Pass Mom and Dad IDs "" to NA
    clutch[(clutch[,16] == "") & (!is.na(clutch[,16])), 16] = NA
    clutch[(clutch[,17] == "") & (!is.na(clutch[,17])), 17] = NA
    #remove clutches for which we don't know one of the parents
    clutch = clutch[(!is.na(clutch[,16])) | (!is.na(clutch[,17])),]

    #Extract clutches for which we have the Size
    clutch = clutch[!is.na(clutch[20]),]

    #create new final dta dataframe
    Clutchdta = as.data.frame(matrix(ncol = 7, nrow = 0))
    colnames(Clutchdta) = c("Survival","ClutchID","MomID","DadID","Year","SiteID", "JulDay")

    #Loop through clutch and create a new dataframe with one observation per egg
    for(rw in 1:nrow(clutch)){

	#extract initial number of eggs
	nbeggs = as.numeric(as.character(clutch[rw,20]))
	#extarct number of eggs which actually hatched
	nbeggsHATCHED = as.numeric(as.character(clutch[rw,52]))
	#extract clutch ID
	clutchID = as.character(clutch[rw,1])
	#extract mother ID
	MotherID = as.character(clutch[rw,17])
	#extract father id
	FatherID = as.character(clutch[rw,16])
	#Extract SiteID
	SiteId = as.character(clutch[rw,3])
	#extract year
	Year = as.character(clutch[rw,2])

	#Extract layingDate
	layDate = clutch[rw,12]
	layDate = as.Date(layDate, format = "%d/%m/%Y")
	#Get JulDay
	JulDay = layDate - as.Date(paste0("01/01/", Year), format = "%d/%m/%Y")

	#create new output dataframe
	output = as.data.frame(matrix(ncol = 7, nrow = nbeggs))
	colnames(output) = c("Survival","ClutchID","MomID","DadID","Year","SiteID", "JulDay")

	#Get 0 and 1 vector of survival (or hatched)
	if(nbeggs == nbeggsHATCHED){

	    survival = rep(1, nbeggs)

	} else {

	    survival = c(rep(1, nbeggsHATCHED), rep(0, (nbeggs - nbeggsHATCHED)))

	}

	#Fill the output
	output[,1] = survival
	output[,2] = clutchID
	output[,3] = MotherID
        output[,4] = FatherID
        output[,5] = Year
        output[,6] = SiteId
        output[,7] = JulDay

	#merge the output with other clutches
	Clutchdta = rbind(Clutchdta, output)

    }

    ###############################################################################
    ############################### NUMBER OF EGGS  ###############################
    ###############################################################################


    ############################### ANIMAL MODEL ###############################

    #update pedigree format for modeling
    p2 <- as.data.frame(ped[,c(1,2,3)])
    ord <- pedigree::orderPed(p2)
    p2 <- p2[order(ord),] #now you can use p2 to fit your animal model

    #merge with dtaF for FEMALES
    names(Clutchdta)[3] = "RingId"
    Clutchdta = merge(Clutchdta, dtaF[,c(1,3)], by = "RingId", all.x = T)
    #Reset MomID name
    names(Clutchdta)[1] = "MotherID"
    #Change mom F names
    names(Clutchdta)[8] = "MotherFuniWE"

    #merge with dtaF for MALES
    names(Clutchdta)[4] = "RingId"
    Clutchdta = merge(Clutchdta, dtaF[,c(1,3)], by = "RingId", all.x = T)
    #Reset MomID name
    names(Clutchdta)[1] = "FatherID"
    #Change mom F names
    names(Clutchdta)[9] = "FatherFuniWE"

    # non - na covariates
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,3]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,7]),]
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,8]),] #mom's F
    Clutchdta <- Clutchdta[!is.na(Clutchdta[,9]),] #dad's F

    #set factor levels with only observations we have
    Clutchdta[,1] = factor(Clutchdta[,1], levels = unique(as.character(Clutchdta[,1]))[!is.na(unique(as.character(Clutchdta[,1])))])
    Clutchdta[,2] = factor(Clutchdta[,2], levels = unique(as.character(Clutchdta[,2]))[!is.na(unique(as.character(Clutchdta[,2])))])
    Clutchdta[,3] = factor(Clutchdta[,3], levels = unique(as.character(Clutchdta[,3]))[!is.na(unique(as.character(Clutchdta[,3])))])
    Clutchdta[,4] = factor(Clutchdta[,4], levels = unique(as.character(Clutchdta[,4]))[!is.na(unique(as.character(Clutchdta[,4])))])
    Clutchdta[,5] = factor(Clutchdta[,5], levels = unique(as.character(Clutchdta[,5]))[!is.na(unique(as.character(Clutchdta[,5])))])
    Clutchdta[,6] = factor(Clutchdta[,6], levels = unique(as.character(Clutchdta[,6]))[!is.na(unique(as.character(Clutchdta[,6])))])

    # as data frame
    Clutchdta <- as.data.frame(Clutchdta)

    #rm weird obs where JulDate < 30
    Clutchdta = Clutchdta[Clutchdta[,7] > 30,]

    #WHEN WE RUN WITH BOTH PARENTS F, ONLY OBSERVATIONS IN 2017 so no Year in the model
    #NO SiteID either btw

    # PRIORS
    # R ARE THE FIXED, G ARE THE RANDOM EFFECT PRIORS
    # ADD PRIORS FOR EVEERY EFFECT
    prior <- list(R = list(V = 1, fix = 1), # default values of inverse wishart for var
                  G = list(G1 = list(V = 1, nu = 0.02, alpha.V = 1000)))

    #Model itself with logit link
    NBEggslaidIDMlogitLINK <- MCMCglmm(Survival ~ JulDay + MotherFuniWE * FatherFuniWE,
                    random = ~ ClutchID,
                    family='categorical',
                    data = Clutchdta,
                    prior = prior,
                    nitt = 5e6,
                    burnin = 50000,
                    thin=1000)

    #Model itself with pobit link for comparison
    NBEggslaidIDMprobitLINK <- MCMCglmm(Survival ~ JulDay + MotherFuniWE * FatherFuniWE,
                    random = ~ ClutchID,
                    family='ordinal',
                    data = Clutchdta,
                    prior = prior,
                    nitt = 5e6,
                    burnin = 50000,
                    thin=1000)

    #save the R Env. for model analysis
    save.image("./5_ID_Models/5.2_modelsOutput/IDM_EggsHatch_BothparentsRFIX1_INT.RData")

EOF
