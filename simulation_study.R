rm(list = ls())

# set working directory
setwd("C:/Users/janni/OneDrive/Documents/VU master/Thesis/code_and_datafiles")

source("data_transformation_functions.R")
source("get_rf_gls.R")
source("parallel_get_grf.R")
source("parallel_get_rf_gls.R")

# load forest package libraries
library(randomForest)
library(grf)

############################## SIMULATED DATASET  ##############################

# load dataset
dataset <- get_sheet_list("dataset_friedmans_d10_s5_scalingfalse.xlsx")
target <- "target"

# initialize empty vectors for prediction results
all_prediction_results_rf <- c()
all_prediction_results_llf <- c()
all_prediction_results_rf_gls <- c()

# Parameters:
B=2000
subsample_frac =0.75 
min_leaf_size=1
min_node_size=5
tn=100
tc=min_node_size-1
m_try=3


# define forests to be tested
forests <- c("regression_forest", "local_linear_forest", "rf_gls")


for (i in 1:100){
  
  # define data
  data <- dataset[[i]]
  true_val <- data[nrow(data),target]
  x0 <- data[nrow(data),!(names(data) %in% c(target))]
  data <- data[-nrow(data),]
  
  print(paste("Simulation i=", i, sep=""))
  
  if ("regression_forest" %in% forests){
    # Obtain results
    grf_output <- get_grf(data=data, 
                          target=target, 
                          B=B, 
                          subsample_size=round(subsample_frac*nrow(data)), 
                          min_leaf_size=min_leaf_size, 
                          min_node_size=min_node_size,
                          forest_type="regression_forest",
                          m_try=m_try,
                          get_ci=TRUE,
                          l=2)
    
    
    results_jannik <- predict_grf(grf_output=grf_output, x=x0)
    pred_jannik <- as.numeric(results_jannik[[1]])
    var_jannik <- as.numeric(results_jannik[[2]])
    
    
    
    grf_athey <- regression_forest( X = data[,!(names(data) %in% c(target))],
                                       Y = data[, target],
                                       num.trees=B,
                                       sample.fraction=(subsample_frac/2),
                                       mtry=m_try,
                                       min.node.size=min_node_size,
                                       honesty=FALSE,
                                       alpha=0,                         # setting to 0 translates into making no imbalance restriction/penalty
                                       ci.group.size = 2,               
                                       tune.parameters = "none")
    results_grf_athey <- predict(grf_athey, newdata=x0, estimate.variance=TRUE)
    pred_grf_athey <- results_grf_athey[[1]]
    var_grf_athey <- results_grf_athey[[2]]
    
    prediction_results <- cbind(pred_jannik, pred_grf_athey, true_val, var_jannik, var_grf_athey)
    
    all_prediction_results_rf <- rbind(all_prediction_results_rf, prediction_results)
    colnames(all_prediction_results_rf) <- c("Jannik_rf", "Athey_rf", "true_value", "Var_jannik", "Var_athey")
  }
  
  if ("rf_rF" %in% forests){
    # Obtain results
    rf_rF <- randomForest(x = data[,!(names(data) %in% c(target))], 
                          y = data[, target], 
                          ntree = B, 
                          mtry = m_try, 
                          replace = FALSE, 
                          sampsize = round(subsample_frac*nrow(data)), 
                          nodesize = min_node_size)
    pred_rF <- as.numeric(predict(object = rf_rF, newdata = x0)) 
    
    prediction_results <- cbind(pred_rF, true_val)
    
    all_prediction_results_rf <- rbind(all_prediction_results_rf, prediction_results)
    colnames(all_prediction_results_rf) <- c("rf_rF", "true_value")
  }
  
  if ("local_linear_forest" %in% forests){
    # Obtain results
    grf_output <- get_grf(data=data, 
                          target=target, 
                          B=B, 
                          subsample_size=round(subsample_frac*nrow(data)), 
                          min_leaf_size=min_leaf_size, 
                          min_node_size=min_node_size,
                          forest_type="local_linear_forest",
                          m_try=m_try,
                          lambda_split=0.1,                          
                          get_ci=TRUE,
                          l=2)
    
    
    results_jannik <- predict_grf(grf_output=grf_output, x=x0, lambda_est=0.5)
    pred_jannik <- as.numeric(results_jannik[[1]])
    var_jannik <- as.numeric(results_jannik[[2]])
    
    
    llf_athey <- ll_regression_forest( X = data[,!(names(data) %in% c(target))],
                                       Y = data[, target],
                                       enable.ll.split = TRUE,
                                       ll.split.lambda = 0.1,           # equivalent to lambda_split
                                       ll.split.cutoff = 0,             # whether to use "local" betas even for small nodes. =0 means yes. 
                                       num.trees=B,
                                       sample.fraction=subsample_frac/2,
                                       mtry=m_try,
                                       min.node.size=min_node_size,
                                       honesty=FALSE,
                                       alpha=0,                         # setting to 0 translates into making no imbalance restriction/penalty
                                       ci.group.size = 2,               # ensures that confidence intervals are not computed, and allowing sample.fraction to be >0.5
                                       tune.parameters = "none")
    results_llf_athey <- predict(llf_athey, newdata=x0, ll.lambda=0.5, estimate.variance=TRUE)
    pred_llf_athey <- results_llf_athey[[1]]
    var_llf_athey <- results_llf_athey[[2]]
    
    prediction_results <- cbind(pred_jannik, pred_llf_athey, true_val, var_jannik, var_llf_athey)
    
    all_prediction_results_llf <- rbind(all_prediction_results_llf, prediction_results)
    colnames(all_prediction_results_llf) <- c("Jannik_LLF", "Athey_LLF", "true_value", "Var_jannik", "Var_athey")
  }
  if ("rf_gls" %in% forests){
    # Obtain results
    rf_gls_output <- get_rf_gls(data=data, 
                                target=target, 
                                B=B, 
                                subsample_size=round(subsample_frac*nrow(data)), 
                                tn=tn, 
                                tc=tc,
                                m_try=m_try,
                                get_ci=TRUE,
                                l=2)
    
    
    results_rf_gls <- predict_rf_gls(rf_gls_output=rf_gls_output, x=x0)
    pred_rf_gls <- as.numeric(results_rf_gls[[1]])
    var_rf_gls <- as.numeric(results_rf_gls[[2]])
    
    
    
    grf_athey <- regression_forest( X = data[,!(names(data) %in% c(target))],
                                    Y = data[, target],
                                    num.trees=B,
                                    sample.fraction=(subsample_frac/2),
                                    mtry=m_try,
                                    min.node.size=min_node_size,
                                    honesty=FALSE,
                                    alpha=0,                         
                                    ci.group.size = 2,               
                                    tune.parameters = "none")
    results_grf_athey <- predict(grf_athey, newdata=x0, estimate.variance=TRUE)
    pred_grf_athey <- results_grf_athey[[1]]
    var_grf_athey <- results_grf_athey[[2]]
    
    prediction_results <- cbind(pred_rf_gls, pred_grf_athey, true_val, var_rf_gls, var_grf_athey)
    
    all_prediction_results_rf_gls <- rbind(all_prediction_results_rf_gls, prediction_results)
    colnames(all_prediction_results_rf_gls) <- c("pred_rf_gls", "Athey_rf", "true_value", "var_rf_gls", "Var_athey")
  }

  if ("regression_forest" %in% forests){
    if (i == 5) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_05.csv", row.names=FALSE)
    if (i == 10) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_1.csv", row.names=FALSE)
    if (i == 20) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_2.csv", row.names=FALSE)
    if (i == 30) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_3.csv", row.names=FALSE)
    if (i == 40) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_4.csv", row.names=FALSE)
    if (i == 50) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_5.csv", row.names=FALSE)
    if (i == 60) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_6.csv", row.names=FALSE)
    if (i == 70) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_7.csv", row.names=FALSE)
    if (i == 80) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_8.csv", row.names=FALSE)
    if (i == 90) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_9.csv", row.names=FALSE)
    if (i == 100) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_10.csv", row.names=FALSE)
  }
  if ("rf_rF" %in% forests){
    if (i == 5) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_05.csv", row.names=FALSE)
    if (i == 10) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_1.csv", row.names=FALSE)
    if (i == 20) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_2.csv", row.names=FALSE)
    if (i == 30) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_3.csv", row.names=FALSE)
    if (i == 40) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_4.csv", row.names=FALSE)
    if (i == 50) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_5.csv", row.names=FALSE)
    if (i == 60) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_6.csv", row.names=FALSE)
    if (i == 70) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_7.csv", row.names=FALSE)
    if (i == 80) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_8.csv", row.names=FALSE)
    if (i == 90) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_9.csv", row.names=FALSE)
    if (i == 100) write.csv(all_prediction_results_rf, "simulation_results_ci_regression_forest/simulation_results_ci_rf_10.csv", row.names=FALSE)
  }
  if ("local_linear_forest" %in% forests){
    if (i == 5) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_05.csv", row.names=FALSE)
    if (i == 10) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_1.csv", row.names=FALSE)
    if (i == 20) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_2.csv", row.names=FALSE)
    if (i == 30) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_3.csv", row.names=FALSE)
    if (i == 40) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_4.csv", row.names=FALSE)
    if (i == 50) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_5.csv", row.names=FALSE)
    if (i == 60) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_6.csv", row.names=FALSE)
    if (i == 70) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_7.csv", row.names=FALSE)
    if (i == 80) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_8.csv", row.names=FALSE)
    if (i == 90) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_9.csv", row.names=FALSE)
    if (i == 100) write.csv(all_prediction_results_llf, "simulation_results_ci_llf/simulation_results_ci_lff_10.csv", row.names=FALSE)
  }
  if ("rf_gls" %in% forests){
    if (i == 5) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_05.csv", row.names=FALSE)
    if (i == 10) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_1.csv", row.names=FALSE)
    if (i == 20) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_2.csv", row.names=FALSE)
    if (i == 30) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_3.csv", row.names=FALSE)
    if (i == 40) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_4.csv", row.names=FALSE)
    if (i == 50) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_5.csv", row.names=FALSE)
    if (i == 60) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_6.csv", row.names=FALSE)
    if (i == 70) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_7.csv", row.names=FALSE)
    if (i == 80) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_8.csv", row.names=FALSE)
    if (i == 90) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_9.csv", row.names=FALSE)
    if (i == 100) write.csv(all_prediction_results_rf_gls, "simulation_results_ci_rf_gls/simulation_results_ci_rf_gls_10.csv", row.names=FALSE)
  }

}

