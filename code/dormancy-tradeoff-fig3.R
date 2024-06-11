####################################################################
## dormancy-tradeoff-fig3.R:  ######################################
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
####################################################################

### Load required packages and project code -----------------------------------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot()) 
source("code/dormancy-tradeoff-functions.R")

### Read individual csv files -------------------------------------------------
sppcombo <- "N1N2N3" # Change for pairwise combinations, e.g. "N1N2")
dat_list <- extract_simfiles(sppcombo)

### Summarize coexistence/exclusion outcomes ----------------------------------
outcomes <- summarise_outcomes(dat_list)

### Make heat map of competitive outcomes -------------------------------------
## Uncomment`load(file = "data/3spp-coexistoutline.Rdata")` to load the outline 
## for the 3 species coexistence region and delete `= NULL` from call to `plot_comp_resource`.
#load(file = "data/3spp-coexistoutline.Rdata") 
simplot <- plot_comp_resource(outcomes, outline = outline) 
simplot

############################

