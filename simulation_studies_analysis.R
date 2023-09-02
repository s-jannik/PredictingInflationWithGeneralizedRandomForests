## Clear workspace
cat("/014")        # clear command window
rm(list = ls())    # clear workspace
dev.off()          # Clear plots

setwd("C:/Users/janni/OneDrive/Documents/VU master/Thesis/code_and_datafiles")
source("data_transformation_functions.R")

dataset <- get_sheet_list("dataset_friedmans_d10_s5_scalingfalse.xlsx")

plot(dataset[[1]][,1], type="l", ylab="", xlab="", ylim=c(-10,40))
for (i in 2:100){
  lines(dataset[[i]][,1])
}


################## Load results from simulations studies #######################

grf_results <- read.csv("simulation_results_ci_far/simulation_results_ci_rf_10_far.csv")
llf_results <- read.csv("simulation_results_ci_far/simulation_results_ci_llf_10_far.csv")
rf_gls_results <- read.csv("simulation_results_ci_far/simulation_results_ci_rf_gls_10_far.csv")
rf_gls_results <- rf_gls_results[, !names(rf_gls_results) %in% c("Athey_rf")]

prediction_results <- cbind(grf_results[1:3], rf_gls_results[,"pred_rf_gls"], llf_results[1:2], grf_results[,"true_value"])
colnames(prediction_results) <- c("rf_j", "rf_a", "rf_rF", "rf_gls", "llf_j", "llf_a", "true_value")

variance_results <- cbind(grf_results[c("Var_jannik", "Var_athey")], rf_gls_results[,"var_rf_gls"], llf_results[c("Var_jannik", "Var_athey")])
colnames(variance_results) <- c("rf_j", "rf_a", "rf_gls", "llf_j", "llf_a")



################################### MSE ########################################

MSE <- matrix(NA, nrow=1, ncol=(ncol(prediction_results)-1))
colnames(MSE) <- colnames(prediction_results)[-ncol(prediction_results)]

for(forest in colnames(MSE)) MSE[1,forest] <- mean((prediction_results[,forest]-prediction_results[,"true_value"])^2)



############################ predictions analysis ###############################

# mean absolute difference between standard forest models 
standard_forests <- c("rf_j", "rf_a", "rf_rF", "rf_gls")
mean_abs_diff <- matrix(NA, nrow=length(standard_forests), ncol=length(standard_forests))
colnames(mean_abs_diff) <- standard_forests
rownames(mean_abs_diff) <- standard_forests
mean_diff <- mean_abs_diff

for(i in standard_forests){
  for(j in standard_forests){
    
    mean_abs_diff[i,j] <- mean(abs(prediction_results[,i]-prediction_results[,j]))
    mean_diff[i,j] <- mean(prediction_results[,i]-prediction_results[,j])
    
  }
}

mean_abs_diff
mean_diff

# mean absolute difference between llf_j and llf_a
mean_abs_diff_llf <- mean(abs(prediction_results[,"llf_a"]-prediction_results[,"llf_j"]))
mean_diff_llf <- mean(prediction_results[,"llf_a"]-prediction_results[,"llf_j"])

mean_abs_diff_llf
mean_diff_llf

############################## variance analysis ###############################

# mean absolute difference between standard forest models 
standard_forests <- c("rf_j", "rf_a", "rf_gls")
var_mean_abs_diff <- matrix(NA, nrow=length(standard_forests), ncol=length(standard_forests))
colnames(var_mean_abs_diff) <- standard_forests
rownames(var_mean_abs_diff) <- standard_forests
var_mean_diff <- var_mean_abs_diff

for(i in standard_forests){
  for(j in standard_forests){
    
    var_mean_abs_diff[i,j] <- mean(abs(variance_results[,i]-variance_results[,j]))
    var_mean_diff[i,j] <- mean(variance_results[,i]-variance_results[,j])
    
  }
}

var_mean_abs_diff
var_mean_diff

# mean absolute difference between llf_j and llf_a
var_mean_abs_diff_llf <- mean(abs(variance_results[,"llf_a"]-variance_results[,"llf_j"]))
var_mean_diff_llf <- mean(variance_results[,"llf_a"]-variance_results[,"llf_j"])

var_mean_abs_diff_llf
var_mean_diff_llf



############################# confidence intervals #############################

is_inside <-  matrix(NA, nrow=nrow(variance_results), ncol=ncol(variance_results))  # true value is inside confidence interval
colnames(is_inside) <- colnames(variance_results)

width

for (forest in colnames(variance_results)){
  
  for (i in 1:nrow(is_inside)){
    
    ci_L <- prediction_results[i, forest] - sqrt(variance_results[i, forest])*1.96
    ci_U <- prediction_results[i, forest] + sqrt(variance_results[i, forest])*1.96
    
    over_lb <- prediction_results[i,"true_value"] > ci_L
    under_ub <- prediction_results[i,"true_value"] < ci_U
    
    is_inside[i, forest] <- over_lb*under_ub
    
  }
}

is_inside

# replace NAs by 0
is_inside[is.na(is_inside)] <- 0

coverage_rates <- apply(is_inside, MARGIN=2, FUN=mean)
coverage_rates



