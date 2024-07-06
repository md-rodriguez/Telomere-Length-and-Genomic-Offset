# Telomere-Length-and-Genomic-Offset
Analyses of telomere length as a response to genomic offset and population declines in the yellow warbler

YEWA_TL_data_2023_noHY.csv contains the metadata needed for the all analyses. The 'sample' column contains unique identifiers for each yellow warbler sampled. The 'TL.Cq' and 'GAP.Cq' columns are the PCR Cq values for telomere (TL) and single copy gene (GAP). The 'TL.eff' and 'GAP.eff' are the qPCR efficiencies for the plate the sample was run on. The 'Plate' column indicates qPCR plate. The 'GO' column is the genomic offset measure. 'TS' is the TS ratio calculated as a measure of telomere length. The 'abundance' column is the measure of abundance trends calculated from Breeding Bird Survey (BBS) data. The '*_Hdiff' columns indicate differences in the respective bioclim variable between historical (1970-2000) and present time (2021-2040). The '*_Fdiff' columns indicate differences in the respective bioclim variables between current (2021-2040) and future time (2041-2060).

The TL.popdiffs.R script includes code for testing differences in telomere length across populations, AIC model selection for ranking models explaining telomere length, and for testing for differences in qPCR Cq and efficiencies across plates and populations.

The pr_explore.R script includes code for analyzing the relationship between climate data (past, present, and future) and genomic offset values (past and future).

