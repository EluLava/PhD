#################################
##### SUPPLEMENTARY FIGURES #####
#################################

library(vioplot)

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/RData/")

###########################################################################
##### FIGURE S1, heritability estimation vioplots adults vs nestlings #####
###########################################################################

#load the models nestlings
load("./Funi_all_model_outputs.RData")
#load the modelsAdults
load("./noF_AnimalModels_ADULTS.RData")

#Create dataframe for h2 estimate
h2est = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(h2est) = c("h2", "trait", "age")

#Get h2 estimates Bill length nestlings
heritBillNest = mod_bill.2.ped[["VCV"]][ , "RingId"] / (mod_bill.2.ped[["VCV"]][ , "RingId"] + mod_bill.2.ped[["VCV"]][ , "RingId_pe"] + mod_bill.2.ped[["VCV"]][ , "Observer"] +
                                                          mod_bill.2.ped[["VCV"]][ , "clutch_merge"] + mod_bill.2.ped[["VCV"]][ , "year"] + mod_bill.2.ped[["VCV"]][ , "units"])
#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritBillNest, trait = rep("BillLength", length(heritBillNest)), age = rep("Nestlings", length(heritBillNest))))
h2est = rbind(h2est, tomerge)

#Get h2 estimates Bill length adults
heritBillAd = BilllengthAM[["VCV"]][ , "animal"] / (BilllengthAM[["VCV"]][ , "animal"] + BilllengthAM[["VCV"]][ , "id"] + BilllengthAM[["VCV"]][ , "Observer"] +
                                                        BilllengthAM[["VCV"]][ , "ClutchIDBorn_Raised"] + BilllengthAM[["VCV"]][ , "YearBorn"] + BilllengthAM[["VCV"]][ , "units"])
#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritBillAd, trait = rep("BillLength", length(heritBillAd)), age = rep("Adults", length(heritBillAd))))
h2est = rbind(h2est, tomerge)

#Get h2 estimates Mass nestlings
heritMassNest = mod_mass.2.ped[["VCV"]][ , "RingId"] / (mod_mass.2.ped[["VCV"]][ , "RingId"] + mod_mass.2.ped[["VCV"]][ , "RingId_pe"] + mod_mass.2.ped[["VCV"]][ , "Observer"] +
                                                          mod_mass.2.ped[["VCV"]][ , "clutch_merge"] + mod_mass.2.ped[["VCV"]][ , "year"] + mod_mass.2.ped[["VCV"]][ , "units"])

#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritMassNest, trait = rep("Mass", length(heritMassNest)), age = rep("Nestlings", length(heritMassNest))))
h2est = rbind(h2est, tomerge)

#Get h2 estimates Bill length adults
heritMassAd = MassAM[["VCV"]][ , "animal"] / (MassAM[["VCV"]][ , "animal"] + MassAM[["VCV"]][ , "id"] + MassAM[["VCV"]][ , "Observer"] +
                                                MassAM[["VCV"]][ , "ClutchIDBorn_Raised"] + MassAM[["VCV"]][ , "YearBorn"] + MassAM[["VCV"]][ , "units"])
#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritMassAd, trait = rep("Mass", length(heritMassAd)), age = rep("Adults", length(heritMassAd))))
h2est = rbind(h2est, tomerge)

#Get h2 estimates Tarsus Length nestlings
heritTarsusNest = mod_tarsus.2.ped[["VCV"]][ , "RingId"] / (mod_tarsus.2.ped[["VCV"]][ , "RingId"] + mod_tarsus.2.ped[["VCV"]][ , "RingId_pe"] + mod_tarsus.2.ped[["VCV"]][ , "Observer"] +
                                                            mod_tarsus.2.ped[["VCV"]][ , "clutch_merge"] + mod_tarsus.2.ped[["VCV"]][ , "year"] + mod_tarsus.2.ped[["VCV"]][ , "units"])

#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritTarsusNest, trait = rep("TarsusLength", length(heritTarsusNest)), age = rep("Nestlings", length(heritTarsusNest))))
h2est = rbind(h2est, tomerge)

#Get h2 estimates Bill length adults
heritTarsusAd = TarsuslengthAM[["VCV"]][ , "animal"] / (TarsuslengthAM[["VCV"]][ , "animal"] + TarsuslengthAM[["VCV"]][ , "id"] + TarsuslengthAM[["VCV"]][ , "Observer"] +
                                                          TarsuslengthAM[["VCV"]][ , "ClutchIDBorn_Raised"] + TarsuslengthAM[["VCV"]][ , "YearBorn"] + TarsuslengthAM[["VCV"]][ , "units"])
#Add it to the dataframe
tomerge = as.data.frame(cbind(h2 = heritTarsusAd, trait = rep("TarsusLength", length(heritTarsusAd)), age = rep("Adults", length(heritTarsusAd))))
h2est = rbind(h2est, tomerge)



#h2 as numeric
h2est$h2 = as.numeric(h2est$h2)
#Traits as factor
h2est$trait = factor(h2est$trait, levels = c("BillLength","Mass","TarsusLength"))
#Age as factor
h2est$age = factor(h2est$age, levels = c("Nestlings","Adults"))

pdf("../../../PLOTS/SUPPLOTS/FIGS1_H2MophoJuvvsAd.pdf", width = 30, height = 10)

#get margins
par(mar = c(15,15,5,5))

vioplot(h2est$h2 ~ h2est$age + h2est$trait, xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL, col = "#9966CC")
#add juv vs ad
axis(1, at = 1:6, labels = rep(c("Juvenile","Adult"),3), padj = 1.75, cex.axis = 2.25, lwd = 2)
#Add traits
mtext(side = 1, text = "Bill Length", at = 1.5, cex = 3, line = 10)
mtext(side = 1, text = "Mass", at = 3.5, cex = 3, line = 10)
mtext(side = 1, text = "Tarsus Length", at = 5.5, cex = 3, line = 10)
#Add y axis title
mtext(side = 2, text =
        expression(V[A] ~ "(heritability)"),
      cex = 2.5, line = 10)
#Add h2 axis
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)

dev.off()

################################################
##### FIGURE S2, beta estimates all models #####
################################################

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/RData/")

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)

#load the models juveniles both F
load("./Funi_all_model_outputs.RData")
load("./FROH_all_model_outputs_Juveniles.RData")
#load the models Adults both F
load("./FuniWE_AnimalModels_ADULTS.RData")
load("./FuniWE_SimpleIDModels_ADULTS.RData")
load("./FHBD_AnimalModels_ADULTS.RData")
load("./FHBD_SimpleIDModels_ADULTS.RData")

F_hat=function(model_input, ibc){
  F_est=as.data.frame(model_input$Sol)%>%select(starts_with(ibc))
  #ibc_eff=F_est*0.25 
  mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
  pred
}


all_f=list(mod_tarsus.1.2,mod_tarsus.2.ped, 
           mod_mass.1.2, mod_mass.2.ped, 
           mod_bill.1.2, mod_bill.2.ped,
           mod_tarsus.1.2.FROH,mod_tarsus.2.ped.FROH, 
           mod_mass.1.2.FROH, mod_mass.2.ped.FROH, 
           mod_bill.1.2.FROH, mod_bill.2.ped.FROH)

all_f_hat_juve=NULL
for (i in all_f){
  
  hat=F_hat(i, "F")
  all_f_hat_juve=rbind(all_f_hat_juve,hat)
}


all_f_hat_juve_df=all_f_hat_juve%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c("Tarsus length.Funi.(simple model)","Tarsus length.Funi.(animal model)", 
                        "Mass.Funi.(simple model)", "Mass.Funi.(animal model)", 
                        "Bill length.Funi.(simple model)", "Bill length.Funi.(animal model)",
                        "Tarsus length.Froh.(simple model)","Tarsus length.Froh.(animal model)", 
                        "Mass.Froh.(simple model)", "Mass.Froh.(animal model)", 
                        "Bill length.Froh.(simple model)", "Bill length.Froh.(animal model)"))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","F", "model"), cols_remove = F)%>%
  mutate(mod_name =as.character(mod_name))

## just getting the models in the right order
all_f_hat_juve_df$mod_name<-ordered(all_f_hat_juve_df$mod_name,
                                    levels = c(
                                      
                                      "Tarsus length.Froh.(animal model)", "Tarsus length.Froh.(simple model)",
                                      "Tarsus length.Funi.(animal model)", "Tarsus length.Funi.(simple model)",
                                      "Mass.Froh.(animal model)", "Mass.Froh.(simple model)", 
                                      "Mass.Funi.(animal model)", "Mass.Funi.(simple model)",
                                      "Bill length.Froh.(animal model)", "Bill length.Froh.(simple model)",
                                      "Bill length.Funi.(animal model)", "Bill length.Funi.(simple model)"
                                      
                                    ))

## just getting the F in the right order
all_f_hat_juve_df$F <- ordered(all_f_hat_juve_df$F, levels = c("Funi","Froh"))

labels = c("Tarsus Length (Froh)\n(animal model)\nn = 6'490", "Tarsus Length (Froh)\n(simple model)\nn = 6'490",
           "Tarsus Length (Funi)\n(animal model)\nn = 6'490", "Tarsus Length (Funi)\n(simple model)\nn = 6'490",
           "Mass (Froh)\n(animal model)\nn = 9'864", "Mass (Froh)\n(simple model)\nn = 9'864",
           "Mass (Funi)\n(animal model)\nn = 9'864", "Mass (Funi)\n(simple model)\nn = 9'864",
           "Bill Length (Froh)\n(animal model)\nn = 6'128", "Bill Length (Froh)\n(simple model)\nn = 6'128",
           "Bill Length (Funi)\n(animal model)\nn = 6'128", "Bill Length (Funi)\n(simple model)\nn = 6'128"
           )

juve_model_com = ggplot(data = all_f_hat_juve_df, aes(x = mod_name, y = est, ymin = upr, ymax = lwr, col = trait, shape = F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2")+ ## actually looks alright just black?? to get cols just change fill to colour
  scale_x_discrete(labels = labels) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() +
  theme(plot.margin=unit(c(1,1,1.25,1), 'cm'),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = -5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 0)) +
  labs(y="Inbreeding depression (\u03B2)", shape="Inbreeding coefficient", col="Trait", 
       x="Model used", title = "Estimates in juveniles")

juve_model_com

### Now Elo's adult models 

all_f_adults=list(TarsuslengthIDM.FuniWE, TarsuslengthAM.FuniWE,
                  TarsuslengthIDM.FHBD, TarsuslengthAM.FHBD,
                  MassIDM.FuniWE, MassAM.FuniWE,
                  MassIDM.FHBD, MassAM.FHBD,
                  BilllengthIDM.FuniWE, BilllengthAM.FuniWE,
                  BilllengthIDM.FHBD, BilllengthAM.FHBD)



all_f_hat_adults = NULL
for (i in all_f_adults){
  
  hat=F_hat(i, "F")
  all_f_hat_adults=rbind(all_f_hat_adults,hat)
}


all_f_hat_adults_df=all_f_hat_adults%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c("Tarsus length.Funi.(simple model)", "Tarsus length.Funi.(animal model)",
                        "Tarsus length.Froh.(simple model)", "Tarsus length.Froh.(animal model)",
                        "Mass.Funi.(simple model)", "Mass.Funi.(animal model)", 
                        "Mass.Froh.(simple model)", "Mass.Froh.(animal model)",
                        "Bill length.Funi.(simple model)", "Bill length.Funi.(animal model)",
                        "Bill length.Froh.(simple model)", "Bill length.Froh.(animal model)"
  ))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","F", "model"), cols_remove = F)%>%
  mutate(mod_name=as.character(mod_name))

## just getting the models in the right order
all_f_hat_adults_df$mod_name<-ordered(all_f_hat_adults_df$mod_name,
                                    levels = c(
                                      
                                      "Tarsus length.Froh.(animal model)", "Tarsus length.Froh.(simple model)",
                                      "Tarsus length.Funi.(animal model)", "Tarsus length.Funi.(simple model)",
                                      "Mass.Froh.(animal model)", "Mass.Froh.(simple model)", 
                                      "Mass.Funi.(animal model)", "Mass.Funi.(simple model)",
                                      "Bill length.Froh.(animal model)", "Bill length.Froh.(simple model)",
                                      "Bill length.Funi.(animal model)", "Bill length.Funi.(simple model)"
                                      
                                    ))

## just getting the F in the right order
all_f_hat_adults_df$F <- ordered(all_f_hat_adults_df$F, levels = c("Funi","Froh"))

labels = c("Tarsus Length (Froh)\n(animal model) n = 773", "Tarsus Length (Froh)\n(simple model) n = 773",
           "Tarsus Length (Funi)\n(animal model) n = 773", "Tarsus Length (Funi)\n(simple model) n = 773",
           "Mass (Froh)\n(animal model) n = 3'677", "Mass (Froh)\n(simple model) n = 3'677",
           "Mass (Funi)\n(animal model) n = 3'677", "Mass (Funi)\n(simple model) n = 3'677",
           "Bill Length (Froh)\n(animal model) n = 2'061", "Bill Length (Froh)\n(simple model) n = 2'061",
           "Bill Length (Funi)\n(animal model) n = 2'061", "Bill Length (Funi)\n(simple model) n = 2'061"
)

adult_model_com = ggplot(data = all_f_hat_adults_df, aes(x = mod_name, y = est, ymin = upr, ymax = lwr, col = trait, shape = F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2")+ ## actually looks alright just black?? to get cols just change fill to colour
  scale_x_discrete(labels = labels) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() +
  theme(plot.margin = unit(c(1,1,1.25,1), 'cm'),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = -5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 0)) +
  labs(y = "Inbreeding depression (\u03B2)", shape = "Inbreeding coefficient", col = "Trait", 
       x = "Model used", title = "Estimates in adults")

adult_model_com

## also eggs 

## have to modify function slightly as interaction means it picks up interaction too
F_hat_2=function(model_input, ibc){
  F_est=as.data.frame(model_input$Sol)%>%select(all_of(ibc))
  #ibc_eff=F_est*0.25 
  mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
  pred
}


hat_eggslaid_Funi = F_hat_2(NBEggslaidIDM.FuniWE, "FuniWE")
hat_eggslaid_Froh = F_hat_2(NBEggslaidIDM.FHBD, "FHBD")
hat_eggslaid_Funi.ped = F_hat_2(NBEggslaidAM.FuniWE, "FuniWE")
hat_eggslaid_Froh.ped = F_hat_2(NBEggslaidAM.FHBD, "FHBD")
hat_eggs_hatching_mum_Funi = F_hat_2(EggsHatchingIDMprobitLINK.FuniWE, "MotherFuniWE")
hat_eggs_hatching_dad_Funi = F_hat_2(EggsHatchingIDMprobitLINK.FuniWE, "FatherFuniWE")
hat_eggs_hatching_mum_Froh = F_hat_2(EggsHatchingIDMprobitLINK.FHBD, "MotherFHBD")
hat_eggs_hatching_dad_Froh = F_hat_2(EggsHatchingIDMprobitLINK.FHBD, "FatherFHBD")


eggies=rbind(hat_eggslaid_Funi,hat_eggslaid_Froh)%>%
  rbind(hat_eggslaid_Funi.ped)%>%
  rbind(hat_eggslaid_Froh.ped)%>%
  rbind(hat_eggs_hatching_mum_Funi)%>%
  rbind(hat_eggs_hatching_dad_Funi)%>%
  rbind(hat_eggs_hatching_mum_Froh)%>%
  rbind(hat_eggs_hatching_dad_Froh)


fhat_eggies=eggies%>%
  rename("est"="V1","upr"="V2","lwr"="V3")%>%
  add_column(mod_name=c(  "# eggs laid.Funi.(simple model)", "# eggs laid.Froh.(simple model)",
                          "# eggs laid.Funi.(animal model)", "# eggs laid.Froh.(animal model)",
                          "Probability of\nEggs hatching \n(Female).Funi.(simple model)", "Probability of\nEggs hatching \n(Male).Funi.(simple model)",
                          "Probability of\nEggs hatching \n(Female).Froh.(simple model)", "Probability of\nEggs hatching \n(Male).Froh.(simple model)"
                          
  ))%>%
  separate_wider_delim(mod_name, delim = ".", names = c("trait","F", "model"), cols_remove = F)%>%
  mutate(mod_name=as.character(mod_name))

## just getting the models in the right order
fhat_eggies$mod_name<-ordered(fhat_eggies$mod_name,
                                      levels = c(
                                        
                                        "# eggs laid.Froh.(animal model)", "# eggs laid.Froh.(simple model)",
                                        "# eggs laid.Funi.(animal model)", "# eggs laid.Funi.(simple model)",
                                        "Probability of\nEggs hatching \n(Male).Froh.(simple model)", "Probability of\nEggs hatching \n(Male).Funi.(simple model)", 
                                        "Probability of\nEggs hatching \n(Female).Froh.(simple model)", "Probability of\nEggs hatching \n(Female).Funi.(simple model)"

                                      ))

## just getting the F in the right order
fhat_eggies$F <- ordered(fhat_eggies$F, levels = c("Funi","Froh"))
## just getting the traits in the right order
fhat_eggies$trait <- ordered(fhat_eggies$trait, levels = c("Probability of\nEggs hatching \n(Female)", "Probability of\nEggs hatching \n(Male)", "# eggs laid"))

fhat_eggies

labels = c("# eggs laid (Froh)\n(animal model)\nn= 489", "# eggs laid (Froh)\n(simple model)\nn= 489",
           "# eggs laid (Funi)\n(animal model)\nn= 489","# eggs laid (Funi)\n(simple model)\nn= 489",
           "Probability of Eggs hatching \n(Male Froh & simple model)\nn = 389", "Probability of Eggs hatching \n(Male Funi & simple model)\nn = 389",
           "Probability of Eggs hatching \n(Female Froh & simple model)\nn = 389", "Probability of Eggs hatching \n(Female Funi & simple model)\nn = 389")

adult_model_com_eggies=ggplot(data=fhat_eggies, aes(x=mod_name, y=est, ymin=upr, ymax=lwr, col=trait, shape=F)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Set2") +
  scale_x_discrete(labels = labels) +
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() +
  theme(plot.margin = unit(c(1,1,1.25,1), 'cm'),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = -5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 0)) +
  labs(y="Inbreeding depression (\u03B2)", shape="Inbreeding coefficient", col="Trait", 
       x="Model used", title = "Estimates in adults (egg traits)")

adult_model_com_eggies

all_comparison=juve_model_com+adult_model_com/adult_model_com_eggies+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')

all_comparison


####################################################################
##### FIGURE S3, Tarsus length in adults with and without rank #####
####################################################################

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/")

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)

#### INCLUDING RANK ####

#load models
load("./IDM_TarsusLength_ADULTSRank.RData")
TarsuslengthIDMrank = TarsuslengthIDM

#### PLOTTING

#Create en empty list of predictors that we'll fill
predictors = list()

# 1. Extract the estimates of the model for each fixed predictor

## Intercept first
Intercept = as.data.frame(matrix(TarsuslengthIDMrank$Sol[,1]))
colnames(Intercept) = "Intercept"

#Add intercept to the predictors list
predictors[[1]] = Intercept
names(predictors)[1] = "Intercept"

## sex, we have two sexes so we'll loop through both

#create dfs
GeneticSex = as.data.frame(cbind(male = vector(mode = "numeric", length = nrow(TarsuslengthIDMrank$Sol)), 
                                 female = vector(mode = "numeric", length = nrow(TarsuslengthIDMrank$Sol))))
#loop through sexes
for(sex in 1:2){
  
  #estimate sex effect size multiplied by sex value (0 for males, 1 for females)
  TMPGeneticSex = as.vector(TarsuslengthIDMrank$Sol[,2])*(sex - 1)
  #Add the intercept effect and fill the dataframe
  GeneticSex[,sex] = TMPGeneticSex + predictors[[1]]
  
}

#Add sex to the predictors list
predictors[[2]] = GeneticSex
names(predictors)[2] = "Intercept_GeneticSex"

## Rank here
# We won't have a regression line per rank individual here (because we don't care (I think)) so we'll get the mean effect for mean rank value

## Rank is the estimate multiplied by mean rank first
TMPRankMEAN = as.vector(TarsuslengthIDMrank$Sol[,3])*mean(TarsuslengthdtaSEQF$Rank)

#Get one estimate including RANK for each sex effect

#create dfs
RankMEAN = as.data.frame(cbind(male = vector(mode = "numeric", length = nrow(TarsuslengthIDMrank$Sol)), 
                               female = vector(mode = "numeric", length = nrow(TarsuslengthIDMrank$Sol))))
#Fill column for males
RankMEAN[,1] = predictors[[2]][,1] + TMPRankMEAN
#Fill column for females
RankMEAN[,2] = predictors[[2]][,2] + TMPRankMEAN

#Add intercept to the predictors list
predictors[[3]] = RankMEAN
names(predictors)[3] = "Intercept_GeneticSex_RankMEAN"

## Inbreeding coefficients, we'll just use the quantiles here as Anna suggested, we have a linear relationship between 
# Would need more (potentially looping through all values of F lol) if the relationship was not linear !

#Create prediction output data frame
predictionsDATA = as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(predictionsDATA) = c("mean", "lwrCI", "uprCI", "FuniWE", "sex")

#estimate quantiles (from the original data not the model)
quant = quantile(TarsuslengthdtaSEQF$FuniWE)

#create dfs
FuniWE = as.data.frame(matrix(nrow = nrow(TarsuslengthIDMrank$Sol), ncol = length(quant)))
colnames(FuniWE) = names(quant)

#loop through quantiles
for(qu in 1:length(quant)){
  
  #multiply the Funi effect size by Funi quantile
  FuniWE[,qu] = TarsuslengthIDMrank$Sol[,4] * quant[qu]
  
  #Loop through levels combination of the previous predictors (here only two since we had two levels of sex and intercept)
  for (lev in 1:ncol(predictors[[length(predictors)]])) {
    
    #combine inbreeding effect with all other effects
    predictionsTMP = unlist(as.vector(predictors[[length(predictors)]][lev] + FuniWE[,qu]))
    
    #get the mean among all iterations
    meanPRED = mean(predictionsTMP)
    #get the lower confidence interval
    lwrCI = quantile(predictionsTMP, probs = c(0.025))
    #get the upper confidence interval
    uprCI = quantile(predictionsTMP, probs = c(0.975))
    
    #Combine three four with other factors levels (sex here) AND F quantile value
    tmpOUTPUT = matrix(c(meanPRED, lwrCI, uprCI, quant[qu], lev), nrow = 1, ncol = 5)
    colnames(tmpOUTPUT) = c("mean", "lwrCI", "uprCI", "FuniWE", "sex")
    
    #Combine with output dataframe
    predictionsDATA = rbind(predictionsDATA, tmpOUTPUT)
    
  }
}

#Pass predictions to numeric
predictionsDATA[,1] = as.numeric(predictionsDATA[,1])
predictionsDATA[,2] = as.numeric(predictionsDATA[,2])
predictionsDATA[,3] = as.numeric(predictionsDATA[,3])
predictionsDATA[,4] = as.numeric(predictionsDATA[,4])
predictionsDATA[,5] = as.factor(predictionsDATA[,5])

#### FINALLY, the plotting ####

PLOT_Tarsus_IDM_Adults_RANK = ggplot(data = predictionsDATA, aes(x = FuniWE, y = mean, col = sex)) +
  scale_color_manual(name = "Sex", labels = c("males", "females"), values = c("#8B1A1A", "#104E8B")) +
  scale_fill_manual(name = "Sex", labels = c("males", "females"), values=c("#8B1A1A", "#104E8B")) +
  geom_point(data = TarsuslengthdtaSEQF, (aes(x = FuniWE, y = LeftTarsus, col = GeneticSex)), inherit.aes = F, alpha = 0.2) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(x = expression("Inbreeding coefficient (" ~ F[UNI] ~ ")"), y = "Tarsus length [mm]") +
  geom_ribbon(aes(ymin = lwrCI, ymax = uprCI, fill = sex), alpha = 0.2, colour = NA) +
  theme(legend.position = "right", text = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 20))) +
  theme(axis.title.y = element_text(margin = margin(r = 20))) +
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dashed"))

PLOT_Tarsus_IDM_Adults_RANK

#### NO RANK ####

load("./IDM_TarsusLength_ADULTSnoRank.RData")
TarsuslengthIDMnorank = TarsuslengthIDM

#### PLOTTING

#Create en empty list of predictors that we'll fill
predictors = list()

# 1. Extract the estimates of the model for each fixed predictor

## Intercept first
Intercept = as.data.frame(matrix(TarsuslengthIDMnorank$Sol[,1]))
colnames(Intercept) = "Intercept"

#Add intercept to the predictors list
predictors[[1]] = Intercept
names(predictors)[1] = "Intercept"

## sex, we have two sexes so we'll loop through both

#create dfs
GeneticSex = as.data.frame(cbind(male = vector(mode = "numeric", length = nrow(TarsuslengthIDMnorank$Sol)), 
                                 female = vector(mode = "numeric", length = nrow(TarsuslengthIDMnorank$Sol))))
#loop through sexes
for(sex in 1:2){
  
  #estimate sex effect size multiplied by sex value (0 for males, 1 for females)
  TMPGeneticSex = as.vector(TarsuslengthIDMnorank$Sol[,2])*(sex - 1)
  #Add the intercept effect and fill the dataframe
  GeneticSex[,sex] = TMPGeneticSex + predictors[[1]]
  
}

#Add sex to the predictors list
predictors[[2]] = GeneticSex
names(predictors)[2] = "Intercept_GeneticSex"

## Inbreeding coefficients, we'll just use the quantiles here as Anna suggested, we have a linear relationship between 
# Would need more (potentially looping through all values of F lol) if the relationship was not linear !

#Create prediction output data frame
predictionsDATA = as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(predictionsDATA) = c("mean", "lwrCI", "uprCI", "FuniWE", "sex")

#estimate quantiles (from the original data not the model)
quant = quantile(TarsuslengthdtaSEQF$FuniWE)

#create dfs
FuniWE = as.data.frame(matrix(nrow = nrow(TarsuslengthIDMnorank$Sol), ncol = length(quant)))
colnames(FuniWE) = names(quant)

#loop through quantiles
for(qu in 1:length(quant)){
  
  #multiply the Funi effect size by Funi quantile
  FuniWE[,qu] = TarsuslengthIDMnorank$Sol[,3] * quant[qu]
  
  #Loop through levels combination of the previous predictors (here only two since we had two levels of sex and intercept)
  for (lev in 1:ncol(predictors[[length(predictors)]])) {
    
    #combine inbreeding effect with all other effects
    predictionsTMP = unlist(as.vector(predictors[[length(predictors)]][lev] + FuniWE[,qu]))
    
    #get the mean among all iterations
    meanPRED = mean(predictionsTMP)
    #get the lower confidence interval
    lwrCI = quantile(predictionsTMP, probs = c(0.025))
    #get the upper confidence interval
    uprCI = quantile(predictionsTMP, probs = c(0.975))
    
    #Combine three four with other factors levels (sex here) AND F quantile value
    tmpOUTPUT = matrix(c(meanPRED, lwrCI, uprCI, quant[qu], lev), nrow = 1, ncol = 5)
    colnames(tmpOUTPUT) = c("mean", "lwrCI", "uprCI", "FuniWE", "sex")
    
    #Combine with output dataframe
    predictionsDATA = rbind(predictionsDATA, tmpOUTPUT)
    
  }
}

#Pass predictions to numeric
predictionsDATA[,1] = as.numeric(predictionsDATA[,1])
predictionsDATA[,2] = as.numeric(predictionsDATA[,2])
predictionsDATA[,3] = as.numeric(predictionsDATA[,3])
predictionsDATA[,4] = as.numeric(predictionsDATA[,4])
predictionsDATA[,5] = as.factor(predictionsDATA[,5])

#### FINALLY, the plotting ####

PLOT_Tarsus_IDM_Adults_noRANK = ggplot(data = predictionsDATA, aes(x = FuniWE, y = mean, col = sex)) +
  scale_color_manual(name = "Sex", labels = c("males", "females"), values = c("#8B1A1A", "#104E8B")) +
  scale_fill_manual(name = "Sex", labels = c("males", "females"), values=c("#8B1A1A", "#104E8B")) +
  geom_point(data = TarsuslengthdtaSEQF, (aes(x = FuniWE, y = LeftTarsus, col = GeneticSex)), inherit.aes = F, alpha = 0.2) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(x = expression("Inbreeding coefficient (" ~ F[UNI] ~ ")"), y = "Tarsus length [mm]") +
  geom_ribbon(aes(ymin = lwrCI, ymax = uprCI, fill = sex), alpha = 0.2, colour = NA) +
  theme(legend.position = "right", text = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 20))) +
  theme(axis.title.y = element_text(margin = margin(r = 20))) +
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dashed"))

PLOT_Tarsus_IDM_Adults_noRANK


bothPLOTS = PLOT_Tarsus_IDM_Adults_RANK + PLOT_Tarsus_IDM_Adults_noRANK +
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')

bothPLOTS
