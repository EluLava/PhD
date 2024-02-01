# scripts

- functions.R: contain different R functions used in the followinf scripts

- 1_Merging_HighandLow_coverages.slurm: merges high-coverage (the refernce panel) and low coverage VCFs. Remove weird individuals we want to exclude from our analyses

- 2_FuniWE_estimation.slurm: Estimated Funi (weighted version) inbreeding coefficient for the 3,085 owls

- 3_ASGRM_estimation.slurm: Estimates the allele-sharing-based genetic relatedness matrix in the 3,085 owls

- 4.1_RZooRoHfriendlyVCFs.slurm: pass VCF from physic to genetic positions and make sure each genetic position has at least 10 units more than the previous position (for model convergence issue)

- 4.2_HBDsegments_detection.slurm: calls HBD segments

- 4.3_HBDsegments_FandDistEstimation.sh: estimates HBD-based inbreeding coefficietn as well as HBD segments distribution

- 4.4_pHBD_perSite.slurm: estimate the probability to be HBD per site and merges it with physical position

- 5_Fasest.slurm: estimates Fas (inbreeding coefficient) for the 3,085 owls

For MCMCglmm models, we call 'animal model (AM)' any script which includes the pedigree relatedness matrix as a random factor and 'inbreeding depression model (IDM)' all modesl which do not include the pedigree relatedness matrix.

- 6.1.1_AM_BillLength_ADULTS.sh: code for running the animal model for bill length in adults

- 6.2.1_AM_Mass_ADULTS.sh: code for running the animal model for mass in adults

- 6.3.1_AM_TarsusLength_ADULTS.sh: code for running the animal model for tarsus length in adults

- 7.1.1_IDM_BillLength_ADULTS.sh : code for running the inbreeding depression model for bill length in adults

- 7.2.1_IDM_Mass_ADULTS.sh: code for running the inbreeding depression model for mass in adults

- 7.3.1_IDM_TarsusLength_ADULTS.sh: code for running the inbreeding depression model for tarsus length in adults

- 7.4_IDM_NBEggsLaid.sh: code for running the inbreeding depression model for number of eggs laid (in adult females only)

- 7.5_IDM_PrEggHatches.sh: code for running the logistic inbreeding depression model for modelling the probabiloty that an egg hatches (in adults only)

All codes for modeling in juveniles are in Anna Hewett GitHub: https://github.com/annamayh/Owls/tree/main/Inbreeding_depression
