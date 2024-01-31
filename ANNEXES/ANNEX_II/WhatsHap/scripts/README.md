# scripts

the following scripts have been used to phase the reference panel with whatsHap

- 1.1_WhatsHap_phasing_readbased.slurm: phases (all) individuals with read-based only (no pedigree information)

- 1.2_WhatsHap_phasing_trios.slurm: phases individuals in trios with both read-bsed and trio (pedigree) information

- 2_postWhatsHap_individuals_merging.slurm: merge outputs individuals from trios are taken from trio-phasing and parents are phased with their offspring with higher coverage. List of individuals in inputFILES/Phasing_merging.list

