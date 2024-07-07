# Telomere-Length-and-Genomic-Offset
Here we combine models of genomic offsets with measures of telomere shortening and abundance trends in a migratory bird, the yellow warbler (Setophaga petechia), 

# YEWA_TL_data_2023_noHY.csv contains the metadata needed for the all analyses. 
The 'TL.Cq' and 'GAP.Cq' columns are the PCR Cq values for telomere (TL) and single copy gene (GAP). The 'TL.eff' and 'GAP.eff' are the qPCR efficiencies for the plate the sample was run on. The 'Plate' column indicates qPCR plate. 

The 'GO' column is the genomic offset measure. To estimate genomic offset, we ran gradient forest (Fitzpatrick et al. 2021), a machine-learning regression tree-based approach implemented in the R package gradientforest (Ellis et al. 2012) on a subset of 1694 unlinked candidate single-nucleotide polymorphisms (SNPs) that were significantly associated with BIOCLIM variables based on LFMM analysis from Bay et al. (2018). We then used the predict function within gradient forest to weight the environmental response variables for both current (2021-2040) and future (2041-2060) predicted climates at 10,000 random locations across the yellow warbler breeding range. We then interpolated across the entire breeding range to form a continuous map of genomic offset.Estimates of genomic offset were then calculated for each specific sample site. 

'TS' is the TS ratio calculated as a measure of telomere length. Telomere length is calculated as the ratio (T/S) of telomere repeat copy number (T) to a control single gene copy number (S), which is standardized to a reference sample and expressed as relative telomere length. We used the glyceraldehyde-3-phosphate dehydrogenase (GAPDH) as the single control gene.

The 'abundance' column is the measure of yellow warbler abundance trend values for each sample location extracted from Breeding Bird Survey (BBS) data (Ziolkowski et al. 2022).

The '*_Hdiff' columns indicate differences in the respective bioclim variable between historical (1970-2000) and present time (2021-2040). The '*_Fdiff' columns indicate differences in the respective bioclim variables between current (2021-2040) and future time (2041-2060).

# The TL.popdiffs.R script includes code for testing differences in telomere length across populations, AIC model selection for ranking models explaining telomere length, and for testing for differences in qPCR Cq and efficiencies across plates and populations.

# The pr_explore.R script includes code for analyzing the relationship between climate data (past, present, and future) and genomic offset values (past and future).

