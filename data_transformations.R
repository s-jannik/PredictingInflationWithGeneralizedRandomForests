## Author: Jannik Steenbergen

## libraries
library(readxl)
library(openxlsx)
library(urca)

## Clear workspace
cat("/014")        # clear command window
rm(list = ls())    # clear workspace
#dev.off()          # Clear plots

## Set working directory
setwd("C:/Users/janni/OneDrive/Documents/VU master/Thesis/code_and_datafiles")

source("data_transformation_functions.R")


## Load data
data <- read.csv("Quarterly_Data_EA_1990Q1_2022Q4.csv")


############################# Data preprocessing ###############################

data <- preprocess_data(data)
ts_data <- transform_data(data=data, take_log=FALSE, n_diff=0, year_on_year=FALSE, as_ts_object=TRUE)


########################### ADF and KPSS test results ##########################

## ADF test results
adf_results <- list("adf_level" = unitroot_test_func(data, take_log=FALSE, n_diff=0, test_type="adf"),
                    "adf_diff" = unitroot_test_func(data, take_log=FALSE, n_diff=1, test_type="adf"),
                    "adf_2diff" = unitroot_test_func(data, take_log=FALSE, n_diff=2, test_type="adf"),
                    "adf_log_level" = unitroot_test_func(data, take_log=TRUE, n_diff=0, test_type="adf"),
                    "adf_log_diff" = unitroot_test_func(data, take_log=TRUE, n_diff=1, test_type="adf"),
                    "adf_log_2diff" = unitroot_test_func(data, take_log=TRUE, n_diff=2, test_type="adf"))
# openxlsx::write.xlsx(adf_results, file="adf_tests_export.xlsx", rowNames=TRUE, showNA=TRUE)


## kpss (lshort=TRUE) test results
kpss_short_results <- list("kpss_level" = unitroot_test_func(data, take_log=FALSE, n_diff=0, test_type="kpss", lshort=TRUE), 
                           "kpss_diff" = unitroot_test_func(data, take_log=FALSE, n_diff=1, test_type="kpss", lshort=TRUE),
                           "kpss_2diff" = unitroot_test_func(data, take_log=FALSE, n_diff=2, test_type="kpss", lshort=TRUE),
                           "kpss_log_level" = unitroot_test_func(data, take_log=TRUE, n_diff=0, test_type="kpss", lshort=TRUE),
                           "kpss_log_diff" = unitroot_test_func(data, take_log=TRUE, n_diff=1, test_type="kpss", lshort=TRUE),
                           "kpss_log_2diff" = unitroot_test_func(data, take_log=TRUE, n_diff=2, test_type="kpss", lshort=TRUE))
# openxlsx::write.xlsx(kpss_short_results, file="kpss_short_tests_export.xlsx", rowNames=TRUE, showNA=TRUE)

## kpss (lshort=FALSE) test results
kpss_long_results <- list("kpss_level" = unitroot_test_func(data, take_log=FALSE, n_diff=0, test_type="kpss", lshort=FALSE), 
                          "kpss_diff" = unitroot_test_func(data, take_log=FALSE, n_diff=1, test_type="kpss", lshort=FALSE),
                          "kpss_2diff" = unitroot_test_func(data, take_log=FALSE, n_diff=2, test_type="kpss", lshort=FALSE),
                          "kpss_log_level" = unitroot_test_func(data, take_log=TRUE, n_diff=0, test_type="kpss", lshort=FALSE),
                          "kpss_log_diff" = unitroot_test_func(data, take_log=TRUE, n_diff=1, test_type="kpss", lshort=FALSE),
                          "kpss_log_2diff" = unitroot_test_func(data, take_log=TRUE, n_diff=2, test_type="kpss", lshort=FALSE))
# openxlsx::write.xlsx(kpss_long_results, file="kpss_long_tests_export.xlsx", rowNames=TRUE, showNA=TRUE)



################################### tcodes #####################################
tcode_results_kpss_short <- adf_kpss_tcode_func(adf_results=adf_results, kpss_results=kpss_short_results, sign_level=0.05, kpss_test_type="Level")
tcode_results_kpss_long <- adf_kpss_tcode_func(adf_results=adf_results, kpss_results=kpss_long_results, sign_level=0.05, kpss_test_type="Level")


################################### Plots ######################################
#plot_data(data=data, take_log=FALSE, n_diff=0, year_on_year=FALSE, purpose="analysis")
#plot_data(data=data, take_log=FALSE, n_diff=0, year_on_year=FALSE, purpose="presentation")


################################### tests for HICP #############################
HICP_data <- transform_data(data=data, take_log=TRUE, n_diff=1, year_on_year=TRUE, as_ts_object=FALSE)
HICP_adf_results <- unitroot_test_func(data=HICP_data, take_log=FALSE, n_diff=0, test_type="adf")
HICP_kpss_short_results <- unitroot_test_func(data=HICP_data, take_log=FALSE, n_diff=0, test_type="kpss", lshort=TRUE)
HICP_kpss_long_results <- unitroot_test_func(data=HICP_data, take_log=FALSE, n_diff=0, test_type="kpss", lshort=FALSE)


############################# pilot preparation ################################

# l12 tcodes
hybrid_tcodes <- tcode_results_kpss_long[,"tcode"]
names(hybrid_tcodes) <- colnames(data)
hybrid_tcodes["HICP"] = 9
hybrid_tcodes["RGDP"] = 8
hybrid_tcodes["LRHHSPR"] = 2

economic_tcodes <- hybrid_tcodes
economic_tcodes[c("M1", "M3")] <- 6
economic_tcodes[c("LRNFCSPR", "LRHHSPR", "LRNFPSSPR")] <- 1
economic_tcodes[c("EATEDSPR")] <- 1
economic_tcodes[c("USDEUROXRATE")] <- 5
economic_tcodes[c("UR")] <- 2

pilot_tcodes <- as.data.frame(cbind(hybrid_tcodes, economic_tcodes))

#openxlsx::write.xlsx(pilot_tcodes, file="pilot_tcodes.xlsx", rowNames=TRUE, showNA=TRUE)



############### Plots of all variables: level and transformed ##################

plot_data(data=data, take_log=FALSE, n_diff=0, year_on_year=FALSE, purpose="presentation")

plot_data_tcode(data=data, tcode=economic_tcodes)
