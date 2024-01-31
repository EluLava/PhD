setwd("~/PhD/Barn_owl/ID_in3Kowls/data/ROHs/")

library(data.table)
library(vioplot)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


###############################################################################################
############################## FIGURE 1: Inbreeding coefficients ##############################

#### READ THE DATA

#FHBD
dtaF = read.table("All3085_FuniWE_FHBD512g.txt", h = T)
#Fas
dtaFas = read.table("../Fas_UNRallfrq_3085.txt", h = T)
#merge both
dtaF = merge(dtaF, dtaFas, by = "INDVs")

#### pHBD along the genome

pHBDerClass = fread("./HBD_prob_per_SNP_largesSS_HBDClasses1to10_PLOTPOS.txt", header = T, verbose = T)

##### HBD segments FIGURE 1 #####

width = c(.3,1,.2,2,.1)
heights = c(.1,1,.25,1.5,.1)

pdf("../PLOTS/FIG1_FHBD.pdf", width = sum(width)*5, height = sum(heights)*5)

layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,3,3,0,rep(0,5)), nrow = 5, ncol = 5, byrow = T),
       widths = width, heights = heights)
#layout.show(n=3)

#### PANEL A; FHBD distrbution vioplot ####

#vioplot
vioplot(dtaF$FHBD512gen, col = "grey", ylim = c(0,0.3), rectCol = "grey20", yaxt = 'n', xaxt = 'n', axes = F)
#add the dots
points(dtaF$FHBD512gen ~ jitter(rep(1,length(dtaF$FHBD512gen)), 3), pch = 20, col = add.alpha("black", .3), cex = 3)
#Add y axis
axis(2, at = seq(0,0.3,0.05), labels = seq(0,0.3,0.05), tick = 3, las = 1, hadj = 1.5, cex.axis = 2.5)
#Add y axis name
mtext(text = expression(F[HBD]), side = 2, at = .15, cex = 3, line = 10)
#Add panel letter
mtext(text = "A", side = 3, at = 0, cex = 4, line = 3)

#### PANEL B; Fas vs FHBD ####

plot(dtaF$FHBD512gen ~ dtaF$Fas, ylab = NA, xlab = NA, pch = 19, col = add.alpha("grey20", .3),
     ylim = c(0,0.35), xlim = c(-0.2,0.5), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 1.5)
abline(0,1)
#x axis
axis(side = 1, at = seq(-20,40,10)/100, labels = seq(-20,40,10)/100, cex.axis = 2.5 , padj = 1.5)
mtext(side = 1, text = expression(F[AS]), cex = 2.5, line = 9, at = 0.1)
#y axis
axis(side = 2, at = seq(0,35,5)/100, labels = seq(0,35,5)/100, cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(F[HBD]), cex = 2.5, line = 9)
#Add panel letter
mtext(side = 3, text = "B", at = -0.4, line = 3, cex = 4)

#### PANEL C; pHBD along the genome ####

plot(pHBDerClass$pHBD ~ pHBDerClass$PLOTPOS, col =  c("grey","black")[((pHBDerClass$COLPOS %% 2) + 1)], type = 'n', axes = F, xlab = NA, ylab = NA, ylim = c(-.25,.5))
#Loop trough ss to add the line
for(ss in unique(pHBDerClass$CHROM)){
  
  #add line
  points(pHBDerClass$pHBD[pHBDerClass$CHROM == ss] ~ pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss], col =  c("grey","black")[((pHBDerClass$COLPOS[pHBDerClass$CHROM == ss] %% 2) + 1)], type = 'l')
  #add line below for ss name
  segments(x0 = min(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss])+(.1*(max(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]) - min(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]))),
           y0 = -.05, x1 = max(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss])-(.1*(max(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]) - min(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]))), y1 = -.05, col = "black", lwd = 2)
  #add s name
  par(xpd = TRUE) #Draw outside plot area
  text(x = (max(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]) - (max(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]) - min(pHBDerClass$PLOTPOS[pHBDerClass$CHROM == ss]))/2), y = -.2, labels = ss, srt = 270, cex = 1.5)
  
}

#add ya xis
axis(2, las = 1, at = seq(0,0.5,0.1), labels = seq(0,0.5,0.1), cex.axis = 2.5)
#Add y axis label
mtext(side = 2, text = expression(p[HBD]), at = 0.25, line = 10, cex = 2.5)
#Add panel letter
mtext(side = 3, text = "C", at = -1.8e8, line = 3, cex = 4)

#close the plot
dev.off()


####################################################################################
############################## FIGURE 2: HBD Segments ##############################

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/ROHs/")

library(data.table)
library(vioplot)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### READ THE DATA

#FHBD
dtaF = read.table("All3085_FuniWE_FHBD512g.txt", h = T)

#HBD segments
dtaDIST = fread("./HBDsegments_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs.hom", header = T, verbose = T)

#Create empty dataframe
dtaSROHNROH = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(dtaSROHNROH) = c("INDVs","SROH","NROH")
#Loop through individuals
for(indv in unique(dtaDIST$id)){
  
  linea=as.data.frame(matrix(nrow = 1, ncol = 3))
  colnames(linea) = colnames(dtaSROHNROH)
  
  linea[1,1] = indv
  linea[1,2] = sum(dtaDIST$length[dtaDIST$id == indv])
  linea[1,3] = nrow(dtaDIST[dtaDIST$id == indv,])    
  
  dtaSROHNROH = rbind(dtaSROHNROH, linea)    
}

dtaSROHNROH$SROH = as.integer(dtaSROHNROH$SROH)
dtaSROHNROH$NROH = as.integer(dtaSROHNROH$NROH)

#### without the less than 10 HBD classes

dtaDIST.4 = dtaDIST[!(dtaDIST$HBDclass %in% c(11,12,13)),]

#Create empty dataframe
dtaSROHNROH.5 = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(dtaSROHNROH.5) = c("INDVs","SROH","NROH")
#Loop through individuals
for(indv in unique(dtaDIST$id)){
  
  linea=as.data.frame(matrix(nrow = 1, ncol = 3))
  colnames(linea) = colnames(dtaSROHNROH.5)
  
  linea[1,1] = indv
  linea[1,2] = sum(dtaDIST.4$length[dtaDIST.4$id == indv])
  linea[1,3] = nrow(dtaDIST.4[dtaDIST.4$id == indv,])
  
  dtaSROHNROH.5 = rbind(dtaSROHNROH.5, linea)    
}

dtaSROHNROH.5$SROH = as.integer(dtaSROHNROH.5$SROH)
dtaSROHNROH.5$NROH = as.integer(dtaSROHNROH.5$NROH)

write.table(dtaSROHNROH.5, "./NROH_SROH_3085owls.txt", quote = F, col.names = T, row.names = F)
dtaSROHNROH.5 = read.table("./NROH_SROH_3085owls.txt", h = T)

#Create empty df to fill
meanSUMofLENgthperClass = as.data.frame(matrix(nrow = 10, ncol = 2))
colnames(meanSUMofLENgthperClass) = c("HBDclass","ALL")
#Fill the HDB class column
meanSUMofLENgthperClass$HBDclass = 1:10

## OVER ALL

# Loop through classes
for(class in 1:10){
  
  #Subset the dataset
  dtasub = dtaDIST[dtaDIST$HBDclass == class,]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST$id)){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$id == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$ALL[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## Extract the 3 less and most inbred individuals
INDVs = dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen)[2])]
INDVs = c(INDVs, dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen)[3])])
INDVs = c(INDVs, dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen)[4])])
INDVs = c(INDVs, dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen, decreasing = T)[3])])
INDVs = c(INDVs, dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen, decreasing = T)[2])])
INDVs = c(INDVs, dtaF$INDVs[which(dtaF$FHBD512gen %in% sort(dtaF$FHBD512gen, decreasing = T)[1])])

HBDseg6INDVs = dtaDIST[dtaDIST$id %in% INDVs,]
#keep only the 10 most recent HBD classes
HBDseg6INDVs = HBDseg6INDVs[HBDseg6INDVs$HBDclass < 11,]
#INDVs as factor to keep order
HBDseg6INDVs$id = factor(HBDseg6INDVs$id, levels = c(INDVs))

##### HBD segments FIGURE 2 #####

width = c(.4,1.5,.4,1.75,.1,1.05,.1)
heights = c(.1,1,.25,1.5,.25)

pdf("../PLOTS/FIG2_HBDsegments.pdf", width = sum(width)*5, height = sum(heights)*5)

layout(mat = matrix(c(rep(0,7),0,1,0,2,0,3,0,rep(0,7),0,4,4,4,4,4,0,rep(0,7)), nrow = 5, ncol = 7, byrow = T),
       widths = width, heights = heights)
#layout.show(n=4)

#### PANEL A; SROH vs NROH ####

plot(dtaSROHNROH.5$NROH ~ dtaSROHNROH.5$SROH, ylab = NA, xlab = NA, pch = 20, xlim = c(0, 3.6e8), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 3,
     ylim = c(0, 750), col = add.alpha("black", .3))
#x axis
axis(side = 1, at = seq(0,3e8,1e8), labels = c(expression("0"),expression("1 x 10"^"8"),expression("2 x 10"^"8"),expression("3 x 10"^"8")), cex.axis = 2.5, padj = 1.5)
mtext(side = 1, text = expression(S["HBD segments [bp]"]), cex = 3, line = 10.5, at = 1.5e8)
#y axis
axis(side = 2, at = seq(0,600,200), labels = seq(0,600,200), cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(N["HBD segments"]), cex = 3, line = 9, at = 300)
#Add panel letter
mtext(side = 3, text = "A", at = -1.3e8, line = 3, cex = 3)

#### PANEL B; HBD segments distribution full pop ####

barplotSistCONTRECENTsmall = meanSUMofLENgthperClass[meanSUMofLENgthperClass$HBDclass <= 7,]
barplotSistCONTRECENTlarge = meanSUMofLENgthperClass[meanSUMofLENgthperClass$HBDclass > 7,]

#Two sets of labels, once with both KB and MB for precise dating but not space, one with only Mb for more space
labelsprecise = c("1g\n[50Mb]","2g\n[25Mb]","4g\n[12.5Mb]","8g\n[6.25Mb]","16g\n[3.13Mb]","32g\n[1.56Mb]","64g\n[781.25Kb]",
                  "128g\n[390.63Kb]","256g\n[195.31Kb]","512g\n[97.66Kb]")
labelsspace = c("1g\n[50]","2g\n[25]","4g\n[12.5]","8g\n[6.25]","16g\n[3.13]","32g\n[1.56]","64g\n[0.78]",
                "128g\n[0.39]","256g\n[0.20]","512g\n[0.10]")

#First four classes
barplot(formula = ALL/1000000 ~ HBDclass, data = barplotSistCONTRECENTsmall, axes = F, col = "grey", xlab = NA, names.arg = NA, ylab = NA)
#add x axis
axis(1, at = seq(0.5,6.5,1) + (0.2*(seq(1,7,1))), labels = labelsspace[1:7], line = 1, padj = 1, cex.axis = 2.5)
#add axis title
mtext(side = 1, text = "# generations back to the coalescence event\n[expected mean HBD segments length in Mb]",
      line = 13, cex = 1.75, at = 7.5)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.25)
#Add axis title
mtext(side = 2, text = "Mean sum of length [Mb]", line = 9, cex = 2)
#Add panel letter
mtext(side = 3, text = "B", at = -2.75, line = 3, cex = 3)

#Last 6 classes
barplot(formula = ALL/1000000 ~ HBDclass, data = barplotSistCONTRECENTlarge, beside = T, axes = F, col = "grey", xlab = NA, names.arg = NA, ylab = NA)
#add x axis
axis(1, at = seq(0.5,2.5,1) + (0.2*(seq(1,3,1))), labels = labelsspace[8:10], line = 1, padj = 1, cex.axis = 2.5)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.5)

#### PANEL C; HBD segments distribution 3 less and most inbred individuals ####

vioplot(HBDseg6INDVs$length ~ HBDseg6INDVs$id, horizontal = T, col = add.alpha("darkgrey", .8), frame.plot = FALSE, yaxt = 'n', xaxt = 'n', axes = F,
        xlab = NULL, ylab = NULL)

#add the dots
#create i for y axiy value for dots
i = 0
#Loop throgh individuals
for (indv in levels(HBDseg6INDVs$id)) {
  
  #Add the dots
  points( i + (jitter(rep(1,length(HBDseg6INDVs$length[HBDseg6INDVs$id == indv])), 3)) ~ HBDseg6INDVs$length[HBDseg6INDVs$id == indv],
          pch = 20, col = add.alpha("black", .5), cex = 3)
  
  #add one value to i
  i = i + 1
  
}

#add y axis labels
axis(side = 2, at = 1:6, labels = INDVs, cex.axis = 3, hadj = 1.25, las = 1)
#add x axis labels
axis(side = 1, at = seq(0, 2.5e7, 5e6), labels = c(expression("0"), expression("0.5 x 10"^"7"), expression("1 x 10"^"7"), expression("1.5 x 10"^"7"),
                                                   expression("2 x 10"^"7"), expression("2.5 x 10"^"7")), cex.axis = 3, padj = 1.5)
#add x axis title
mtext(side = 1, text = "HBD segment length [bp]", cex = 2.75, line = 10, at = 1.25e7)
#Add panel letter
mtext(side = 3, text = "C", at = -0.35e7, line = 3, cex = 3)

dev.off()


######################################################################################
############################## FIGURE 3: Variance plots ##############################

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/RData/")

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)

#load the models
load("./noF_AnimalModels_ADULTS.RData")
load("./Funi_all_model_outputs.RData")

#define function to extract variance comonents
get_var_partitions = function(model_input){
  
  tar_randoms=as.data.frame(posterior.mode(model_input$VCV))
  
  total_V_tar=as.numeric(colSums(tar_randoms))
  va=as.numeric(tar_randoms[1,]/total_V_tar)
  
  pe=as.numeric(tar_randoms[2,]/total_V_tar)
  obs=as.numeric(tar_randoms[3,]/total_V_tar)
  clutch=as.numeric(tar_randoms[4,]/total_V_tar)
  year=as.numeric(tar_randoms[5,]/total_V_tar)
  res=as.numeric(tar_randoms[6,]/total_V_tar)
  
  
  Var_ped=rbind(va,pe)%>%
    rbind(obs)%>%
    rbind(clutch)%>%
    rbind(year)%>%
    rbind(res)%>%
    as.data.frame()%>%
    tibble::rownames_to_column("Variance_element")
  
  Var_ped
  
}

#tarsus length

flede_tarsus=get_var_partitions(mod_tarsus.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Tarsus Length")
adult_tarsus=get_var_partitions(TarsuslengthAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Tarsus Length")

#bill length

flede_bill=get_var_partitions(mod_bill.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Bill Length")

adult_bill=get_var_partitions(BilllengthAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Bill Length")

#Mass

flede_mass=get_var_partitions(mod_mass.2.ped)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Mass")

adult_mass=get_var_partitions(MassAM)%>%
  mutate(age="Adult")%>%
  mutate(trait="Mass")

all_var_comp=rbind(flede_tarsus, adult_tarsus)%>%
  rbind(flede_bill)%>%rbind(adult_bill)%>%
  rbind(flede_mass)%>%rbind(adult_mass)%>%
  mutate()

all_var_comp$Variance_element<-as.factor(all_var_comp$Variance_element)
all_var_comp$Variance_element<-ordered(all_var_comp$Variance_element, 
                                       levels = c("year", "pe", "clutch", 
                                                  "obs","res","va"))



var_plot_lengths=ggplot(all_var_comp, aes(fill=Variance_element, x=factor(age, level = c('Juvenile', 'Adult')), y=V1)) +
  geom_bar(position="stack",stat="identity")+
  theme_bw() +
  facet_wrap(~trait)+
  scale_fill_viridis_d(direction = -1, 
                       labels=c(bquote (V["Birth Year"]), 
                                bquote (V["Individual ID"]),
                                bquote (V["Clutch ID"]),
                                bquote (V["Observer ID"]),
                                bquote (V["Residuals"]),
                                bquote (V["A"])
                       )) +
  theme(text = element_text(size = 15), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15)
  ) +
  labs(fill="Variance \ncomponents")


var_plot_lengths



########################################################################################################
############################## FIGURE 4: Beta estimates for simple models ##############################

setwd("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/RData/")

library(tidyverse)
library(matrixStats)
library(viridis)
library(MCMCglmm)
library(patchwork)

#load the models
load("./FuniWE_SimpleIDModels_ADULTS.RData")
load("./Funi_all_model_outputs.RData")


estimate_F_effect=function(model_input){
  F_est=as.data.frame(model_input$Sol)%>%dplyr::select(matches("FuniWE"))
  
  #ibc_eff=F_est*-0.25 
  mean<-apply(F_est,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(F_est,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(F_est,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  pred=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))
  pred
}


flede_tarsus_est=estimate_F_effect(mod_tarsus.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Tarsus Length")
adult_tarsus_est=estimate_F_effect(TarsuslengthIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Tarsus Length")

flede_mass_est=estimate_F_effect(mod_mass.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Mass")
adult_mass_est=estimate_F_effect(MassIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Mass")

flede_bill_est=estimate_F_effect(mod_bill.1.2)%>%
  mutate(age="Juvenile")%>%
  mutate(trait="Bill Length")
adult_bill_est=estimate_F_effect(BilllengthIDM.FuniWE)%>%
  mutate(age="Adult")%>%
  mutate(trait="Bill Length")

all_ests_F=rbind(flede_tarsus_est, adult_tarsus_est)%>%
  rbind(flede_mass_est)%>%rbind(adult_mass_est)%>%
  rbind(flede_bill_est)%>%rbind(adult_bill_est)

all_ests_F$trait = ordered(all_ests_F$trait, 
                          levels = c("Tarsus Length", "Mass", "Bill Length"))

all_ests_F$age = ordered(all_ests_F$age, 
                           levels = c("Adult", "Juvenile"))

labels = c("Tarsus Length\nn.juv = 6'490\nn.ad = 773",
                    "Mass\nn.juv = 9'864\nn.ad = 3'677",
                      "Bill Length\nn.juv = 6'128\nn.ad = 2'061")

beta_ests_Funi = ggplot(data = all_ests_F, aes(x = trait, y = V1, ymin = V2, ymax = V3, shape = age)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_colour_brewer(palette = "Dark2") + ## actually looks alright just black?? to get cols just change fill to colour
  scale_x_discrete(labels = labels) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() +
  theme(plot.margin=unit(c(1,1,1.25,1), 'cm'),
        title = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = -3),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = -2),
        axis.text.y = element_text(hjust = 0),
        legend.text =  element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "left") +
  labs(y="Inbreeding depression (\u03B2)", shape="Life stage", title = "Adult and juvenile traits")+
  guides(colour = "none")+
  geom_text(aes(x = 3, y = -5, label = "ns"),colour="black", size = 3) +
  geom_text(aes(x = 3.35, y = -17, label = "*** \n(p < 0.01)"), colour = "black", size = 3) +
  geom_text(aes(x = 2, y = 18, label = "ns"),colour="black", size = 3) +
  geom_text(aes(x = 2.35, y = -28, label = ". \n(p = 0.08)"), colour = "black", size = 3) +
  geom_text(aes(x = 1, y = 31, label = "ns"),colour="black", size = 3) +
  geom_text(aes(x = 1.35, y = -28, label = ". \n(p = 0.07)"), colour = "black", size = 3)


  beta_ests_Funi

##and for eggs 

number_eggs = estimate_F_effect(NBEggslaidIDM.FuniWE)%>%
  mutate(age = "# eggs laid")%>%
  mutate(trait = "# eggs laid")


F_est_mother=as.data.frame(EggsHatchingIDMprobitLINK.FuniWE$Sol)%>%dplyr::select(matches("MotherFuniWE"))
#ibc_eff=F_est*-0.25 
mean<-apply(F_est_mother,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(F_est_mother,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(F_est_mother,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
pred_mother=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))%>%
  add_column(age="Probability of \nEggs hatching \n(Female F)")%>%
  add_column(trait="Probability of \nEggs hatching")


F_est_father=as.data.frame(EggsHatchingIDMprobitLINK.FuniWE$Sol)%>%dplyr::select(matches("FatherFuniWE"))
#ibc_eff=F_est*-0.25 
mean<-apply(F_est_father,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(F_est_father,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(F_est_father,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
pred_father=as.data.frame(matrix(c(mean,CI_upper,CI_lower), nrow = 1, ncol = 3))%>%
  add_column(age="Probability of \nEggs hatching \n(Male F)")%>%
  add_column(trait="Probability of \nEggs hatching")

eggs=pred_father%>%
  rbind(pred_mother)%>%
  rbind(number_eggs)
eggs

labels = c("# eggs laid\nn= 489", "Probability of \nEggs hatching \n(Female Funi)\nn = 389", "Probability of \nEggs hatching \n(Male Funi)\nn = 389")

beta_ests_Funi_eggs=ggplot(data=eggs, aes(x=age, y=V1, ymin=V2, ymax=V3)) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  # scale_colour_brewer(palette = "Greens", direction = -1)+ ## actually looks alright just black?? to get cols just change fill to colour
  scale_x_discrete(labels = labels) +
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip()+
  theme(title = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = -3),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = -2),
        axis.text.y = element_text(hjust = 0)) +        
  labs(y="Inbreeding depression (\u03B2)", title = "Adult only traits")+
  guides(color = "none")+
  geom_text(aes(x=1.15,y=0.9 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=2.15,y=1.8 ,label="ns"),colour="black",size=3)+
  geom_text(aes(x=3.15,y=4 ,label="ns"),colour="black",size=3)

#Combine both

all_beta_traits=beta_ests_Funi+beta_ests_Funi_eggs+ plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')

all_beta_traits

#### REPEATABILITY ####

library(ICC)

## Bill LENGTH ADULTS

load("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/AM_BillLength_ADULTS_noRank.RData")

ICC::ICCest(x = id, y = BillLength, data = Billlengthdta)

load("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/AM_Mass_ADULTS_noRank.RData")

ICC::ICCest(x = id, y = CalculatedMass, data = Massdta)

load("~/PhD/Barn_owl/ID_in3Kowls/data/MCMCglmm_MODELS/SendingAnna/AM_Mass_ADULTS_noRank.RData")

ICC::ICCest(x = id, y = CalculatedMass, data = Massdta)
