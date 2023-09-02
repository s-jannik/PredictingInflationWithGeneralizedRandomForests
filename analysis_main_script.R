## Clear workspace
cat("/014")        # clear command window
rm(list = ls())    # clear workspace
dev.off()          # Clear plots

library(openxlsx)

# set working directory
setwd("C:/Users/janni/OneDrive/Documents/VU master/Thesis/code_and_datafiles")

# source files
source("data_transformation_functions.R")
source("main_script_functions.R")

library(forecast) # for diebold-mariano test

################################################################################

# import results 
all_results <- get_sheet_list(file="færdige_resultater/main_run_ci/horizon1_ci_samlet.xlsx")

# define results in environment
predictions <- all_results[["predictions"]]
variances <- all_results[["variances"]]
o_data <- all_results[["o_data"]]
vi_predictions <- all_results
vi_predictions[[c("predictions")]] <- NULL       # removes predictions from vi_predictions
vi_predictions[[c("variances")]] <- NULL       # removes predictions from vi_predictions
vi_predictions[[c("o_data")]] <- NULL            # removes o_data from vi_predictions


################################################################################
## re-define some objects
list_models <- colnames(predictions)[!colnames(predictions)=="true_val"]
list_forests <- names(vi_predictions)
target="HICP"



############################# 1 - Predictions ##################################

# MSE 
MSE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(MSE) <- list_models

for (forest in list_models){
  MSE[1, forest] <- mean((predictions[,"true_val"] - predictions[, forest])^2)
}

# RMSE 
RMSE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(RMSE) <- list_models

for (forest in list_models){
  RMSE[1, forest] <- sqrt(mean((predictions[,"true_val"] - predictions[, forest])^2))
}

# Mean error
ME <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(ME) <- list_models

for (forest in list_models){
  ME[1, forest] <- mean(predictions[,"true_val"] - predictions[, forest])
}

# Mean absolute error
MAE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(MAE) <- list_models

for (forest in list_models){
  MAE[1, forest] <- mean(abs(predictions[,"true_val"] - predictions[, forest]))
}


################################# RMSE ratios ##################################

#s <- c("regression_forest", "local_linear_forest", "rf_gls", "ar1_benchmark", "ar4_benchmark")
#RMSE[,c(s, "rw_benchmark")]

RMSE[, list_models]/RMSE[,"rw_benchmark"]
MAE[, list_models]/MAE[,"rw_benchmark"]


################################# fractions ####################################

errors <- predictions[,1:6] - predictions[,"true_val"]
sq_errors <- errors^2
sq_errors <- sq_errors[,list_models]


# fraction of time model is worst performing model

fractions <- matrix(NA, ncol=ncol(sq_errors), nrow=1)
colnames(fractions) <- colnames(sq_errors)
for (i in colnames(sq_errors)){
  
  fractions[1, i] <- sum(colnames(sq_errors)[apply(sq_errors, MARGIN=1, FUN=which.max)] == i)/nrow(sq_errors)

}


# fraction of time model is best performing model

fractions_2 <- matrix(NA, ncol=ncol(sq_errors), nrow=1)
colnames(fractions_2) <- colnames(sq_errors)
for (i in colnames(sq_errors)){
  
  fractions_2[1, i] <- sum(colnames(sq_errors)[apply(sq_errors, MARGIN=1, FUN=which.min)] == i)/nrow(sq_errors)
  
}

# fraction of time forest is best performing forest
forests <- c("regression_forest", "local_linear_forest", "rf_gls")
fractions_4 <- matrix(NA, ncol=length(forests), nrow=1)
colnames(fractions_4) <- forests
sq_errors_forests <- sq_errors[,forests]
for (i in forests){
  
  fractions_4[1, i] <- sum(forests[apply(sq_errors_forests, MARGIN=1, FUN=which.min)] == i)/nrow(sq_errors_forests)
  
}

# fraction of times forest outperforms standard forest

fractions_5 <- matrix(NA, ncol=(ncol(sq_errors_forests)-1), nrow=1)
colnames(fractions_5) <- colnames(sq_errors_forests)[!colnames(sq_errors_forests) %in% "regression_forest"]
for (i in colnames(fractions_5)){
  
  sq_errors_short <- sq_errors_forests[,c(i, "regression_forest")]
  #print(sq_errors)
  
  fractions_5[1, i] <- sum(colnames(sq_errors_short)[apply(sq_errors_short, MARGIN=1, FUN=which.min)] == i)/nrow(sq_errors_short)
  
}

# fraction of times model outperforms RW

fractions_3 <- matrix(NA, ncol=(ncol(sq_errors)-1), nrow=1)
colnames(fractions_3) <- colnames(sq_errors)[!colnames(sq_errors) %in% "rw_benchmark"]
for (i in colnames(fractions_3)){
  
  sq_errors_short <- sq_errors[,c(i, "rw_benchmark")]
  #print(sq_errors)
  
  fractions_3[1, i] <- sum(colnames(sq_errors_short)[apply(sq_errors_short, MARGIN=1, FUN=which.min)] == i)/nrow(sq_errors_short)
  
}

fractions
fractions_2
fractions_3
fractions_4
fractions_5

################################## DM tests ####################################


DM_stats <- c()

for (i in c(1,2,3,4,6)){
  DM_stats[i] <- dm.test(e1=errors[,i], e2=errors[,5], alternative = "greater", h=1, power=1, varestimator="acf")$p.value
}
names(DM_stats) <- colnames(errors)

DM_stats


############################# Predictions plot #################################

# prediction start and end date
pred_start <- as.numeric(strsplit(gsub("\\D", " ", rownames(predictions)[1]), "  ")[[1]])
pred_end <- as.numeric(strsplit(gsub("\\D", " ", rownames(predictions)[nrow(predictions)]), "  ")[[1]])

ts_predictions <- ts(predictions, frequency=4, start=pred_start, end=pred_end)

# convert original data to ts object for plotting
ts_o_data <- transform_data(data=o_data, take_log=FALSE, n_diff=FALSE, year_on_year=FALSE, as_ts_object=TRUE)
plot(ts_o_data[,target], ylab="", xlab="")

par(mfrow=c(2,1))

colours <- c("regression_forest"="red", "local_linear_forest"="blue", "rf_gls"="green", 
             "rw_benchmark"="cyan", "ar1_benchmark"="brown", "ar4_benchmark"="orange")

par(mfrow=c(1,1))
# plot from specific date
plot(window(ts_o_data[,target], start=c(2012, 4)), ylab="", xlab="")
for(model in (list_models[!list_models %in% c("rw_benchmark", "ar4_benchmark")])){
  lines(ts_predictions[,model], col=colours[model])
}



############################ Cumulative sum plot ###############################

par(mfrow=c(1,1))

benchmark <- "rw_benchmark"

cumsum_mat <- matrix(NA, nrow=nrow(predictions), ncol=(ncol(predictions)-2))
colnames(cumsum_mat) <- list_models[!(list_models %in% "rw_benchmark")]

for(model in list_models[!(list_models %in% benchmark)]){
  
  # compute cumulative sum of squared errors
  se <- (predictions[,"true_val"] - predictions[,model])^2
  cumsum <- cumsum(se)
  
  # compute for benchmark
  se_bm <- (predictions[,"true_val"] - predictions[,benchmark])^2
  cumsum_bm <- cumsum(se_bm)
  
  cumsum_mat[,model] <- cumsum-cumsum_bm
}
  
ts_cumsum <- ts(cumsum_mat, start=pred_start, frequency=4)
  
ylim <- max(abs(max(ts_cumsum)), abs(min(ts_cumsum)))
ylim <- c((1-ylim*1.1), (1+ylim*1.1))

plot(ts_cumsum[,"regression_forest"], ylim=ylim, main="", xlab="", ylab="", col=colours["regression_forest"])
abline(h=c(1), col=c("black"), lty=c(1), lwd=c(1))
for(model in list_models[!(list_models %in% c("rw_benchmark", "regression_forest", "ar4_benchmark"))]){
 lines(ts_cumsum[,model], col=colours[model]) 
}

  



################################################################################
######################### 2 - Confidence intervals #############################
################################################################################

# Negative variances
negative_var_indices <- is.na(sqrt(variances[,list_forests]))
negative_var_indicator <- apply(negative_var_indices, MARGIN=1, FUN=sum) != 0

####### Coverage rates (recording negative variance as 0)
coverage_rates <- matrix(NA, nrow=1, ncol=length(list_forests))
colnames(coverage_rates)<-list_forests

for (forest in list_forests){
  
  variance <- variances[,forest]
  lb <- predictions[,forest] - sqrt(variance)*1.95
  ub <- predictions[,forest] + sqrt(variance)*1.95
  
  over_lb <- predictions[,"true_val"] > lb
  under_ub <- predictions[,"true_val"] < ub
  
  in_interval <- over_lb*under_ub
  
  in_interval[is.na(in_interval)] <- 0
  
  coverage_rates[,forest] <- sum(in_interval)/length(in_interval)
  
}


####### Coverage rates (only recording for all obs with no negative variances)
coverage_rates <- matrix(NA, nrow=1, ncol=length(list_forests))
colnames(coverage_rates)<-list_forests

for (forest in list_forests){
  
  variance <- variances[,forest]
  lb <- predictions[,forest] - sqrt(variance)*1.95
  ub <- predictions[,forest] + sqrt(variance)*1.95
  
  over_lb <- predictions[,"true_val"] > lb
  under_ub <- predictions[,"true_val"] < ub
  
  in_interval <- over_lb*under_ub
  
  in_interval <- in_interval[!negative_var_indicator]
  
  coverage_rates[,forest] <- sum(in_interval)/length(in_interval)
  
}


####### interval width
interval_width <- matrix(NA, nrow=1, ncol=length(list_forests))
colnames(interval_width)<-list_forests

negative_var_indices <- is.na(sqrt(variances[,list_forests]))
negative_var_indicator <- apply(negative_var_indices, MARGIN=1, FUN=sum) != 0


for (forest in list_forests){
  
  variance <- variances[,forest]
  lb <- predictions[,forest] - sqrt(variance)*1.95
  ub <- predictions[,forest] + sqrt(variance)*1.95
  
  width <- ub-lb
  width <- width[!negative_var_indicator]
  
  
  interval_width[,forest] <- mean(width)
  
}

coverage_rates
interval_width


####### Confidence interval plots
plot(variances[,"regression_forest"], type="l")

par(mfrow=c(1,1))
# plot from specific date
plot(window(ts_o_data[,target], start=c(2012, 4)), ylab="", xlab="")
for(model in (list_models[!list_models %in% c("rw_benchmark", "ar4_benchmark")])){
  
  # add the predictions
  lines(ts_predictions[,model], col=colours[model])
  
  if (model %in% c("regression_forest", "rf_gls")){
    
    # add the confidence intervals 
    variance <- variances[,model]
    lb <- ts_predictions[,model] - sqrt(variance)*1.95
    ub <- ts_predictions[,model] + sqrt(variance)*1.95
    
    lines(lb, col=colours[model], lty=3)
    lines(ub, col=colours[model], lty=3)
    
  }
  
}


################################################################################
########################### 3 - Variable importance ############################
################################################################################

################### Plot version 1 - individual factor importance ##############

# RMSE vi_predictions 
vi_RMSE <- matrix(NA, nrow=length(list_forests), ncol=ncol(vi_predictions[[1]]))
rownames(vi_RMSE) <- list_forests
colnames(vi_RMSE) <- colnames(vi_predictions[[1]])

for (feature in 1:ncol(vi_RMSE)){
  for (forest in list_forests){
    vi_RMSE[forest, feature] <- sqrt(mean((predictions[,"true_val"] - vi_predictions[[forest]][, feature])^2))
  }
}

# Final variable importance measure 
vi_percentage_gains <- matrix(NA, nrow=length(list_forests), ncol=ncol(vi_predictions[[1]]))
rownames(vi_percentage_gains) <- list_forests
colnames(vi_percentage_gains) <- colnames(vi_predictions[[1]])

for (feature in 1:ncol(vi_RMSE)){
  
  # interpreted as "percentage decrease in RMSE from including the predictor as opposed to not to include it"
  for (forest in list_forests){
    vi_percentage_gains[forest, feature] <- -1 * (RMSE[1, forest] - vi_RMSE[forest, feature]) / vi_RMSE[forest, feature]*100
  }
}

# change variable names
colnames(vi_percentage_gains)[colnames(vi_percentage_gains) == "HICP_l1"] <- "Y lag 1"
colnames(vi_percentage_gains)[colnames(vi_percentage_gains) == "HICP_l2"] <- "Y lag 2"
colnames(vi_percentage_gains)[colnames(vi_percentage_gains) == "HICP_l3"] <- "Y lag 3"
colnames(vi_percentage_gains)[colnames(vi_percentage_gains) == "HICP_l4"] <- "Y lag 4"

# Variable importance plots 
par(mfrow=c(1,3))
titles <- c("RF", "LLF", "RF-GLS")
for (forest_nr in 1:length(list_forests)){
  
  sorted_vi_gains <- sort(vi_percentage_gains[forest_nr,], decreasing=TRUE)[1:10]
  sorted_vi_gains <- sort(sorted_vi_gains)
  
  barplot(sorted_vi_gains, horiz=TRUE, col="red", cex.names=1.1, las=1, main=titles[forest_nr], xlim=range(pretty(c(0, sorted_vi_gains))))
}



###################### Variable importance plots version 2 #####################
vi_groups <- list()

vi_groups[["AR"]] <- c("Y lag 1", "Y lag 2", "Y lag 3", "Y lag 4")
vi_groups[["output"]] <- c("f1")
vi_groups[["money - credit"]] <- c("f2", "f3", "f4", "f5")
vi_groups[["interest rates"]] <- c("f6", "f7", "f8")
vi_groups[["exchange rates"]] <- c("f10")
vi_groups[["stock market"]] <- c("f9", "f11")
vi_groups[["financial conditions"]] <- c("f12")
vi_groups[["stress"]] <- c("f13")
vi_groups[["prices"]] <- c("f14")
vi_groups[["industrial production"]] <- c("f15")
vi_groups[["employment"]] <- c("f16")
vi_groups[["confidence"]] <- c("f19", "f20")
vi_groups[["PMI"]] <- c("f21")
vi_groups[["other"]] <- c("f17", "f18")

bar_colors <- c("AR"="red", "output"="blue", "money - credit"="steelblue1", 
                "interest rates"="green", "exchange rates"="darkgreen", "stock market"="yellow", 
                "financial conditions"="orange", "stress"="cyan", "prices"="pink", "industrial production"="darkorchid", 
                "employment"="brown", "confidence"="gray87", "PMI"="gray45", "other"="black")



vi_group_percentage_gains <- c()

for (i in 1:length(vi_groups)){
  
  factors <- vi_groups[[i]]
  
  if (length(factors) == 1){
    
    new_vi <- vi_percentage_gains[, factors]
    
  }
  else {
    new_vi <- apply(vi_percentage_gains[, factors], 1, sum)
  }

  vi_group_percentage_gains <- cbind(vi_group_percentage_gains, new_vi)
  
}
colnames(vi_group_percentage_gains) <- names(vi_groups)


# set negative gains to 0
vi_group_percentage_gains[sign(vi_group_percentage_gains)==-1] <- 0



titles <- c("RF", "LLF", "RF-GLS")
for (forest_nr in 1:length(list_forests)){
  
  sorted_vi_gains <- sort(vi_group_percentage_gains[forest_nr,], decreasing=TRUE)[1:10]
  sorted_vi_gains <- sort(sorted_vi_gains)
  
  cols <- bar_colors[names(sorted_vi_gains)]
  
  barplot(sorted_vi_gains, names.arg="", horiz=TRUE, col=cols, cex.names=1.1, las=1, main=titles[forest_nr], xlim=range(pretty(c(0, sorted_vi_gains))))
}


###################### Variable importance plots version 3 #####################


# standardize such that sum is 1
for (forest in rownames(vi_group_percentage_gains)){
  
  vi_group_percentage_gains[forest,] <- vi_group_percentage_gains[forest,]/sum(vi_group_percentage_gains[forest,])
  
}

# sort 
sorted_vi_group_percentage_gains <- vi_group_percentage_gains

par(mfrow=c(1,1))
barplot(t(sorted_vi_group_percentage_gains), border=bar_colors, col=bar_colors, names.arg=c("RF", "LLF", "RF-GLS"), args.legend = list(bty = "n", x = "bottom", ncol = 3))

# Obtain legend separately 
plot(0,0, ylim=c(2,100))
legend(x = -1, y=50, legend = colnames(sorted_vi_group_percentage_gains), fill = bar_colors, ncol = 4, inset = 1, cex=0.89)






################################################################################
######################### 4- Treating inflation is I(1) ########################
################################################################################

predictions_new <- matrix(NA, nrow=nrow(predictions), ncol=ncol(predictions))
colnames(predictions_new) <- colnames(predictions)
rownames(predictions_new) <- rownames(predictions)

inflation_raw <- read.csv("inflation_data.csv")
inflation <- as.data.frame(inflation_raw[,"HICP"])
rownames(inflation) <- inflation_raw[,1]
colnames(inflation) <- "HICP"

# Correct predictions 

today <- c(2013,4)    # initialize today

for (t in 1:nrow(predictions)){
  for (model in list_models){
  
    last_inflation_obs <- window_func(inflation, from=today, to=today)[1,1]
    predictions_new[t, model] <- last_inflation_obs + predictions[t, model]
    
  }
  
  oos_date <- as.numeric(strsplit(gsub("\\D", " ", rownames(predictions)[t]), "  ")[[1]])
  
  if (oos_date[1] != 2023) predictions_new[t,"true_val"] <- window_func(inflation, from=oos_date, to=oos_date)[1,1]
  
  today <- update_date(today)
}



## rename predictions_new predictions
predictions <- predictions_new
predictions <- as.data.frame(predictions)


################################# Predictions ##################################

# MSE 
MSE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(MSE) <- list_models

for (forest in list_models){
  MSE[1, forest] <- mean((predictions[,"true_val"] - predictions[, forest])^2)
}

# RMSE 
RMSE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(RMSE) <- list_models

for (forest in list_models){
  RMSE[1, forest] <- sqrt(mean((predictions[,"true_val"] - predictions[, forest])^2))
}

# Mean error
ME <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(ME) <- list_models

for (forest in list_models){
  ME[1, forest] <- mean(predictions[,"true_val"] - predictions[, forest])
}

# Mean absolute error
MAE <- matrix(NA, ncol=length(list_models), nrow=1)
colnames(MAE) <- list_models

for (forest in list_models){
  MAE[1, forest] <- mean(abs(predictions[,"true_val"] - predictions[, forest]))
}


s <- c("regression_forest", "local_linear_forest", "rf_gls", "ar1_benchmark", "ar4_benchmark")

# I(1) results
RMSE[, list_models]/RMSE[,"rw_benchmark"]
RMSE[,c(list_models, "rw_benchmark")]
MAE[, list_models]/MAE[,"rw_benchmark"]
