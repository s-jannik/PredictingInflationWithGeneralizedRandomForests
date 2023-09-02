## Clear workspace
cat("/014")        # clear command window
rm(list = ls())    # clear workspace
# dev.off()          # Clear plots


# libraries
library(openxlsx)

# set working directory
setwd("C:/Users/janni/OneDrive/Documents/VU master/Thesis/code_and_datafiles")

# import functions 
source("main_script_functions.R")
source("data_transformation_functions.R")

source("get_rf_gls.R")
source("parallel_get_grf.R")
source("parallel_get_rf_gls.R")



############################# 0 - Import data ##################################

# (original) data
o_data <- read.csv("Quarterly_Data_EA_1990Q1_2022Q4.csv")
o_data <- preprocess_data(data=o_data)

# tcodes (transformation codes)
pilot_tcodes <- openxlsx::read.xlsx("pilot_tcodes.xlsx", colNames=TRUE, rowNames=TRUE)
tcodes <- pilot_tcodes$economic_tcodes
names(tcodes) <- rownames(pilot_tcodes)


################################################################################

# Settings
horizons = c(1,2,4,8,12)   
max_lag = 4     # max number of lags to include in data

# Some data re-definitions
# apply variable transformation and obtain dataset including all lags
o_data <- transform_data_tcode(data=o_data, tcode=tcodes, use_which_data="all_available", as_ts_object=FALSE)
o_data <- get_prediction_data(data=o_data, max_lag=max_lag, h=0)    
target <- "HICP"

for (horizon in horizons){
  
  ############################ 1 - Define data #################################
  
  # save full sample target data
  df_y_full_sample <- na.omit(o_data[,target, drop=FALSE])
  
  # omit all quarters with NA values
  o_data <- na.omit(o_data)
  
  # group predictors: dependent variable lags vs. other predictors
  ylags_list <- c()
  for (i in 1:max_lag){
    ylags_list[i] <- paste("HICP_l", i, sep="")
  }
  other_predictors_list <- colnames(o_data)[!(names(o_data) %in% c(target, ylags_list))]
  
  # divide other predictors into groups 
  lvl_var_names <- (gsub('.{3}$', '', other_predictors_list))
  
  groups <- list()
  groups[["output"]] <- other_predictors_list[lvl_var_names %in% c("RGDP")]
  groups[["money_vol"]] <- other_predictors_list[lvl_var_names %in% c("M1", "M3")]
  groups[["credit_vol"]] <- other_predictors_list[lvl_var_names %in% c("LONFC", "LOHH", "LONFPS", "CRNFPS")]
  groups[["money_vol_ratio"]] <- other_predictors_list[lvl_var_names %in% c("M12GDP", "M32GDP")]
  groups[["credit_vol_ratio"]] <- other_predictors_list[lvl_var_names %in% c("LONFC2GDP", "LOHH2GDP", "LONFPS2GDP", "CRNFPS2GDP")]
  groups[["credit_spread"]] <- other_predictors_list[lvl_var_names %in% c("LRNFCSPR", "LRHHSPR", "LRNFPSSPR", "CRSPR")]
  groups[["yield_curve"]] <- other_predictors_list[lvl_var_names %in% c("YC3M", "YC1Y", "YC2Y", "YC3Y", "YC2YW", "YCEONIA")]
  groups[["interest_rate_spread"]] <- other_predictors_list[lvl_var_names %in% c("EATEDSPR")]
  groups[["asset_prices"]] <- other_predictors_list[lvl_var_names %in% c("STP", "HP")]
  groups[["exchange_rates"]] <- other_predictors_list[lvl_var_names %in% c("USDEUROXRATE", "NEER")]
  groups[["volatility_indicator"]] <- other_predictors_list[lvl_var_names %in% c("VSTOXX")]
  groups[["FCI_indicator"]] <- other_predictors_list[lvl_var_names %in% c("FCIBME", "FCIVAR", "FCIIPA", "FCIBBG", "FCIGS")]
  groups[["comp_stress_indicator"]] <- other_predictors_list[lvl_var_names %in% c("CISS", "CISSnew", "SRI")]
  groups[["commodity_prices"]] <- other_predictors_list[lvl_var_names %in% c("OILPUSD", "GASPEURUSD", "RAWMATPUSD")]
  groups[["industrial_production"]] <- other_predictors_list[lvl_var_names %in% c("IP", "CAPUTILMANUF")]
  groups[["employment"]] <- other_predictors_list[lvl_var_names %in% c("EMPL", "UR")]
  groups[["trade"]] <- other_predictors_list[lvl_var_names %in% c("RETTRADE")]
  groups[["business_cycle"]] <- other_predictors_list[lvl_var_names %in% c("EUROCOIN")]
  groups[["industrial_confidence"]] <- other_predictors_list[lvl_var_names %in% c("INDCONF")]
  groups[["consumer_confidence"]] <- other_predictors_list[lvl_var_names %in% c("CONCONF")]
  groups[["purch_man_indices"]] <- other_predictors_list[lvl_var_names %in% c("PMIMAN", "PMICONF")]
  
  
  ############################## 2 - POOS setup ##################################
  
  ######################## 2.1 User-specified settings ###########################
  # list of forests to compare
  list_forests <- c("regression_forest", "local_linear_forest", "rf_gls")
  list_models <- c("rw_benchmark", "ar1_benchmark", "ar4_benchmark", list_forests)
  
  # Parameters and hyperparameters
  to_be_tuned=NULL       # = NULL if no parameters should be tuned
  B = 2000
  subsample_frac = 0.75
  min_leaf_size = 1
  min_node_size = 5
  m_try = 7
  lambda_split = 0.1
  lambda_est = 0.5 
  tn = 100
  tc = (min_node_size - 1)
  get_ci = TRUE
  l = 2
  
  
  # Full sample start and end date
  data_start <- as.numeric(strsplit(gsub("\\D", " ", rownames(o_data)[1]), "  ")[[1]])
  data_end <- as.numeric(strsplit(gsub("\\D", " ", rownames(o_data)[nrow(o_data)]), "  ")[[1]])
  
  
  # Initialize time windows
  train_set_start <- data_start   # initialize
  target_start <- train_set_start   # date corresponding to first lead of target
  for (h_idx in 1:(horizon)) target_start <- update_date(target_start)
  
  #today <- c(2015, 4)
  today <- c(2013, 4)
  
  oos_test_set_start <- today     # initialization
  for (h_idx in 1:(horizon)) oos_test_set_start <- update_date(oos_test_set_start)
  target_end <- oos_test_set_start      # for clarity 
  
  termination_day <- data_end    # termination_day: last day at which a prediction is made
  for (h_idx in 1:(horizon)) termination_day <- downdate_date(termination_day)
  
  prediction_start_day <- oos_test_set_start        # for the purpose of plotting results
  
  # print initial dates 
  cat("train_set_start:", train_set_start,
      "\ntarget_start:", target_start,
      "\ntoday:", today,
      "\noos_test_set_start:", oos_test_set_start, 
      "\ntarget_end:", target_end, 
      "\ntermination_day:", termination_day,
      "\n")
  

  
  ############################ 2.2 POOS procedure ################################
  
  # Obtain all factor-based datasets
  factor_information <- get_all_factor_datasets(o_data=o_data, target=target, ylags_list=ylags_list, 
                                                other_predictors_list=other_predictors_list, 
                                                groups=groups, df_y_full_sample=df_y_full_sample, 
                                                train_set_start=train_set_start, today=today,
                                                target_start=target_start, target_end=target_end, 
                                                termination_day=termination_day)
  all_datasets <- factor_information[[1]]
  
  # Define full out-of-sample test set 
  oos_test_set <- factor_information[[2]]
  
  # Prepare variable importance
  permuted_oos_test_set <- apply(oos_test_set, MARGIN=2, FUN=sample, replace=FALSE)[,!(names(oos_test_set) %in% c(target))]
  rownames(permuted_oos_test_set) <- rownames(oos_test_set)
  
  
  # Initialize empty matrices for prediction results
  predictions <- matrix(NA, ncol=length(list_models)+1, nrow=nrow(oos_test_set))
  colnames(predictions) <- c(list_models, "true_val")
  rownames(predictions) <- rownames(oos_test_set)
  
  # Initialize empty matrix for variance estimates
  variances <- matrix(NA, ncol=length(list_models)+1, nrow=nrow(oos_test_set))
  colnames(variances) <- c(list_models, "true_val")
  rownames(variances) <- rownames(oos_test_set)
  
  # Initialize empty matrices for variable importance results
  vi_predictions <- list()  # variable importance predictions 
  vi_mat <- matrix(NA, nrow=nrow(permuted_oos_test_set), ncol=(ncol(permuted_oos_test_set)))
  colnames(vi_mat) <- colnames(permuted_oos_test_set)
  rownames(vi_mat) <- rownames(permuted_oos_test_set)
  for (forest in list_forests) vi_predictions[[forest]] <- vi_mat
  

  ############################ Start of POOS procedure ###########################
  
  # initialize counter
  counter <- 0 
  
  while (today[1] != update_date(termination_day)[1] | today[2] != update_date(termination_day)[2]){
  
    counter <- counter + 1
  
    print(paste("Horizon:", horizon, "/", "Date:", today[1], today[2], "/", "oos test point nr.", counter, "of", nrow(oos_test_set)), sep=" ")
    
    # define (the updated) full in-sample dataset
    is_data <- all_datasets[[get_dates(from=today, to=today)]]
    
    
    ############################### OOS prediction ###############################
    
    # define oos test point
    x <- window_func(oos_test_set[, !(names(oos_test_set) %in% c(target))], 
                     from=oos_test_set_start, to=oos_test_set_start)    
  
    # obtain and store predictions 
    if ("regression_forest" %in% list_models){
      
      rf <- get_grf(target = target, 
                    data = is_data, 
                    B = B,
                    subsample_size = round(subsample_frac*nrow(is_data)),
                    min_leaf_size = min_leaf_size, 
                    min_node_size = min_node_size, 
                    m_try = m_try,
                    forest_type = "regression_forest",
                    get_ci = get_ci,
                    l=l)
      
      pred_results <- predict_grf(grf_output = rf, x = x)
      predictions[counter, "regression_forest"] <- as.numeric(pred_results[1])
      variances[counter, "regression_forest"] <- as.numeric(pred_results[2])
    }
    
    if ("local_linear_forest" %in% list_models){
      
      llf <- get_grf(target = target, 
                    data = is_data, 
                    B = B,
                    subsample_size = round(subsample_frac*nrow(is_data)),
                    min_leaf_size = min_leaf_size, 
                    min_node_size = min_node_size, 
                    m_try = m_try,
                    forest_type = "local_linear_forest",
                    lambda_split = lambda_split,
                    get_ci=get_ci,
                    l=l)
      
      pred_results <- predict_grf(grf_output = llf, x = x, lambda_est = lambda_est)
      predictions[counter, "local_linear_forest"] <- as.numeric(pred_results[1])
      variances[counter, "local_linear_forest"] <- as.numeric(pred_results[2])
      
    }
    if ("rf_gls" %in% list_models){
      
      rf_gls <- get_rf_gls(target = target, 
                           data = is_data, 
                           B = B, 
                           subsample_size = round(subsample_frac*nrow(is_data)), 
                           tn = tn, 
                           tc = tc, 
                           m_try = m_try,
                           get_ci = get_ci,
                           l=l)
      
      pred_results <- predict_rf_gls(rf_gls_output = rf_gls, x = x)
      predictions[counter, "rf_gls"] <- as.numeric(pred_results[1])
      variances[counter, "rf_gls"] <- as.numeric(pred_results[2])
      
    }
    
    if ("rw_benchmark" %in% list_models){
      
      predictions[counter, "rw_benchmark"] <- is_data[nrow(is_data), target]
      
    }
    
    if ("ar1_benchmark" %in% list_models){
      
      ar_coefs <- lm(is_data[, target] ~ is_data[, ylags_list[1]])$coefficients
      intercept <- as.numeric(ar_coefs[1])
      phi_1 <- as.numeric(ar_coefs[2])
      predictions[counter, "ar1_benchmark"] <- intercept + phi_1*is_data[nrow(is_data),target]
      
    }
    
    if ("ar4_benchmark" %in% list_models){
      
      ar_coefs <- lm(is_data[, target] ~ is_data[, ylags_list[1]] 
                                              + is_data[, ylags_list[2]] 
                                              + is_data[, ylags_list[3]]
                                              + is_data[, ylags_list[4]])$coefficients
      intercept <- as.numeric(ar_coefs[1])
      phi_1 <- as.numeric(ar_coefs[2])
      phi_2 <- as.numeric(ar_coefs[3])
      phi_3 <- as.numeric(ar_coefs[4])
      phi_4 <- as.numeric(ar_coefs[5])
      predictions[counter, "ar4_benchmark"] <- intercept + (phi_1*is_data[nrow(is_data),target] 
                                                            + phi_2*is_data[(nrow(is_data)-1),target]
                                                            + phi_3*is_data[(nrow(is_data)-2),target]
                                                            + phi_4*is_data[(nrow(is_data)-3),target])
      
    }

    # true value 
    predictions[counter, "true_val"] <- as.numeric(window_func(oos_test_set[, target, drop=FALSE], from=oos_test_set_start, to=oos_test_set_start))  
    
    ######################## Variable Importance predictions #####################
    
    for (feature in 1:ncol(permuted_oos_test_set)){
      
      # adjust x such that 1 feature is permuted
      x_p <- x
      x_p[1, feature] <- permuted_oos_test_set[counter, feature]
      
      if ("regression_forest" %in% list_forests) vi_predictions[["regression_forest"]][counter, feature] <- as.numeric(predict_grf(grf_output = rf, x = x_p)[1])
      if ("local_linear_forest" %in% list_forests) vi_predictions[["local_linear_forest"]][counter, feature] <- as.numeric(predict_grf(grf_output = llf, x = x_p, lambda_est = lambda_est)[1])
      if ("rf_gls" %in% list_forests) vi_predictions[["rf_gls"]][counter, feature] <- as.numeric(predict_rf_gls(rf_gls_output = rf_gls, x = x_p)[1])
    }
  
    ################################ Update dates ################################ 
    
    today <- update_date(today)
    oos_test_set_start <- update_date(oos_test_set_start)
    target_end <- update_date(target_end)
  }
  
  ############################ End of POOS procedure #############################
  
  ## Export results 
  all_results <- list()
  all_results[["predictions"]] <- predictions
  all_results[["variances"]] <- variances
  all_results[["o_data"]] <- o_data
  all_results <- append(all_results, vi_predictions)
  
  if (horizon==1) openxlsx::write.xlsx(all_results, file="main_results/horizon=1_rf_ci_economictcodes.xlsx", rowNames=TRUE, showNA=TRUE)
  if (horizon==2) openxlsx::write.xlsx(all_results, file="main_results/horizon=2_rf_ci_economictcodes.xlsx", rowNames=TRUE, showNA=TRUE)
  if (horizon==4) openxlsx::write.xlsx(all_results, file="main_results/horizon=4_rf_ci_economictcodes.xlsx", rowNames=TRUE, showNA=TRUE)
  if (horizon==8) openxlsx::write.xlsx(all_results, file="main_results/horizon=8_rf_ci_economictcodes.xlsx", rowNames=TRUE, showNA=TRUE)
  if (horizon==12) openxlsx::write.xlsx(all_results, file="main_results/horizon=12_rf_ci_economictcodes.xlsx", rowNames=TRUE, showNA=TRUE)
  
}

