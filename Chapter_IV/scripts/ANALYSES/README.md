

- functions.R: contain different R functions used in the followinf scripts

- 1_Merging_HighandLow_coverages.slurm: merges high-coverage (the refernce panel) and low coverage VCFs. Remove weird individuals we want to exclude from our analyses

- 2_FuniWE_estimation.slurm: Estimated Funi (weighted version) inbreeding coefficient for the 3,085 owls

- 3_ASGRM_estimation.slurm: Estimates the allele-sharing-based genetic relatedness matrix in the 3,085 owls

- 4.1_RZooRoHfriendlyVCFs.slurm: pass VCF from physic to genetic positions and make sure each genetic position has at least 10 units more than the previous position (for model convergence issue)

- 4.2_HBDsegments_detection.slurm: calls HBD segments

- 4.3_HBDsegments_FandDistEstimation.sh: estimates HBD-based inbreeding coefficietn as well as HBD segments distribution

- 4.4_pHBD_perSite.slurm: estimate the probability to be HBD per site and merges it with physical position

- 5_Fasest.slurm: estimates Fas (inbreeding coefficient) for the 3,085 owls

- 
