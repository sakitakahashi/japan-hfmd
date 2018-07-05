source("Code/Functions_TSIR.R")

data_raw <- read.csv("Data/data_raw.csv")

###########################################
## Get data for 1-serotype TSIR analysis ##
###########################################

## Global parameters
which_proportion <- "WMA"
start_year_infer <- 1997
max_year <- 2015

## Save EV-A71 data
which_serotype <- "EVA71"
data_EVA71 <- Est_UR_1sero_linear(data=data_raw, which_serotype=which_serotype, which_proportion=which_proportion, start_year_infer=start_year_infer, max_year=max_year)

## Save CV-A16 data
which_serotype <- "CVA16"
data_CVA16 <- Est_UR_1sero_linear(data=data_raw, which_serotype=which_serotype, which_proportion=which_proportion, start_year_infer=start_year_infer, max_year=max_year)

save(data_EVA71, data_CVA16, start_year_infer, max_year, file="data_1sero_1997.RData")

###########################################
## Get data for 2-serotype TSIR analysis ##
###########################################

## Global parameters
which_proportion <- "WMA"
min_year <- 1982
start_year_infer <- 1997
last_year_included <- 2006
max_year <- 2015
k_EVA71 <- 8
k_CVA16 <- 39

## Complete data set
data_both <- Est_UR_2sero_linear(data=data_raw, which_proportion=which_proportion, min_year=min_year, start_year_infer=start_year_infer, max_year=max_year, k_EVA71=k_EVA71, k_CVA16=k_CVA16)

## Training set
data_both_train <- Est_UR_2sero_linear(data=data_raw, which_proportion=which_proportion, min_year=min_year, start_year_infer=start_year_infer, max_year=last_year_included, k_EVA71=k_EVA71, k_CVA16=k_CVA16)

## Testing set
data_both_test <- Est_UR_2sero_linear(data=data_raw, which_proportion=which_proportion, min_year=min_year, start_year_infer=start_year_infer, max_year=max_year, k_EVA71=k_EVA71, k_CVA16=k_CVA16)

save(data_both, data_both_train, data_both_test, min_year, start_year_infer, last_year_included, max_year, k_EVA71, k_CVA16, file="data_2sero_1997.RData")
