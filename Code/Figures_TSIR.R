library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("Code/Functions_TSIR.R")

###################
## Make Figure 3 ##
###################

load("Data/data_1sero_1997.RData")

Fig_3A <- Run_TSIR_1sero(data=data_EVA71, start_year_infer=start_year_infer, max_year=max_year, alpha=0.975, color="darkseagreen", which_plot_return="beta")
Fig_3B <- Run_TSIR_1sero(data=data_EVA71, start_year_infer=start_year_infer, max_year=max_year, alpha=0.975, color="darkseagreen", which_plot_return="fwdsim")
Fig_3C <- Run_TSIR_1sero(data=data_CVA16, start_year_infer=start_year_infer, max_year=max_year, alpha=0.975, color="firebrick", which_plot_return="beta")
Fig_3D <- Run_TSIR_1sero(data=data_CVA16, start_year_infer=start_year_infer, max_year=max_year, alpha=0.975, color="firebrick", which_plot_return="fwdsim")
blank <- grid.rect(gp=gpar(col="white"))

par(mar=c(2,2,0.5,0.5), oma=c(0.5,0.5,0.5,0.5))
grid.arrange(Fig_3A, Fig_3B, blank, Fig_3C, Fig_3D, blank, blank, blank, layout_matrix=rbind(c(1,1,2,2,2,2,3,3), c(1,1,2,2,2,2,7,7), c(4,4,5,5,5,5,6,6), c(4,4,5,5,5,5,8,8)))

###################
## Make Figure 5 ##
###################

load("Data/data_2sero_1997.RData")

par(mar=c(2,2,0.5,0.5), oma=c(0.5,0.5,0.5,0.5))
Run_TSIR_2sero(data=data_both, data_train=data_both_train, data_test=data_both_test, start_year_infer=start_year_infer, max_year=max_year, k_EVA71=k_EVA71, k_CVA16=k_CVA16, alpha=0.975)
