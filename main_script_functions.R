source("pca.R")


update_date <- function(old_date){
  
  old_year <- old_date[1]
  old_quarter <- old_date[2]
  
  if (old_quarter == 4) return(c(old_year+1, 1))
  
  else if (old_quarter < 4) return(c(old_year, old_quarter+1))
  
  else print("Warning. Date not updated correctly.")
  
}

downdate_date <- function(old_date){
  
  old_year <- old_date[1]
  old_quarter <- old_date[2]
  
  if (old_quarter == 1) return(c(old_year-1, 4))
  
  else if (old_quarter > 1) return(c(old_year, old_quarter-1))
  
  else print("Warning. Date not updated correctly.")
  
}


get_dates <- function(from, to){
  # Description
    # simple function returning a vector of all dates from X to X at quarterly frequency
  # Arguments
    # from: start date in form c(1990, 1)
    # to: end date in form c(1990, 1)
  
  if (from[1] > to[1]){
    print("ERROR: from is before to.")
    return(NA)
  }
  if (from[1] == to[1]){
    if (from[2] > to[2]){
      print("ERROR: from is before to.")
      return(NA)
    }
  }
  
  year_quarter <- matrix(from, ncol=2, nrow=1)
  
  # initialize running date at from
  running_date <- from 
  
  while(running_date[1] != to[1] | running_date[2] != to[2]){
    
    running_date <- update_date(running_date)
    
    year_quarter <- rbind(year_quarter, running_date)
    
    
  }
  
  dates <- c()
  for (i in 1:nrow(year_quarter)){
    dates[i] <- paste(year_quarter[i,1], " Q", year_quarter[i,2] , sep="")
  }
  return(dates)
}

window_func <- function(data, from, to){
  # Description
    # function extracting a window of data from a dataframe with rownames of form "1990 Q1"
  # Arguments
    # data: dataframe with rownames of form "1990 Q1"
    # from: starting quarter of form c(1990, 1)
    # to: ending quarter of form c(1990, 1)
  
  year_quarter <- matrix(NA, ncol=2, nrow=nrow(data))
  
  for (i in 1:nrow(data)){
    date <- as.numeric(strsplit(gsub("\\D", " ", rownames(data)[i]), "  ")[[1]])
    year_quarter[i,1] <- date[1]
    year_quarter[i,2] <- date[2]
  }
  
  # if desired window is not in the dataset, print error message
  if (from[1] < year_quarter[1,1]){
    print("ERROR in window_func: Desired window not available in data.")
  }
  if (from[1] == year_quarter[1,1]){
    if (from[2] < year_quarter[1,2]){
      print("ERROR in window_func: Desired window not available in data.")
    }
  }
  if (to[1] > year_quarter[nrow(year_quarter),1]){
    print("ERROR in window_func: Desired window not available in data.")
  }
  if (to[1] == year_quarter[nrow(year_quarter),1]){
    if (to[2] > year_quarter[nrow(year_quarter),2]){
      print("ERROR in window_func: Desired window not available in data.")
    }
  }
  
  from_idx <- which((from[1]==year_quarter[,1])*(from[2]==year_quarter[,2]) == 1)
  to_idx <- which((to[1]==year_quarter[,1])*(to[2]==year_quarter[,2]) == 1)
  
  data_new <- data[from_idx:to_idx, , drop=FALSE]
  
  return(data_new)
}


dimensionality_reduction <- function(o_data, target, ylags_list, other_predictors_list, groups, df_y_full_sample, today, target_start, target_end){
  
  # Arguments
    # o_data: all data available at a given date
    # other_predictors_list: vector of the names of other predictors in o_data
    # groups: list of grouped variable names
  
  ######### Center and standardize data
  
  # Centering and standardizing predictors excl. ylags
  X <- as.matrix(o_data[,c(other_predictors_list)])
  X <- scale(X, center=TRUE, scale=TRUE)
  
  
  ######### Group predictors
  factors <- c() 
  loadings <- list()
  factor_names <- c()
  for (group_nr in 1:length(groups)){
    
    pca_results <- dfm_pca(X=X[,groups[[group_nr]]], r=1)
    group_factor <- pca_results$F
    
    group_loadings <- pca_results$L
    names(group_loadings) <- groups[[group_nr]]
    
    factors <- cbind(factors, group_factor)
    loadings[[group_nr]] <- group_loadings
    
    # add factor name
    factor_names <- c(factor_names, paste("f", group_nr, sep=""))
    
  }
  colnames(factors) <- factor_names
  names(loadings) <- factor_names
  
  ## construct full dataset with lagged regressors up until and including first oos test point
  
  # obtain target leads
  target_data_leads <- window_func(df_y_full_sample, from=target_start, to = target_end)
  
  # data including an out-of-sample part
  data_incl_oos <- cbind(target_data_leads, o_data[, ylags_list], factors)
  
  # extract test point 
  oos_point <- data_incl_oos[nrow(data_incl_oos), , drop=FALSE]   # note: all features in oos_point are in-sample
  
  # obtain final dataset to be fed to forest 
  data <- window_func(data_incl_oos, from=target_start, to=today)  
    # note: features include information starting from train_set_start
    # note: the dates in "data" refer to the data at which HICP is observed 
  
  return(list("data"=data, "oos_point"=oos_point, "loadings"=loadings))
}


get_all_factor_datasets <- function(o_data, target, ylags_list, other_predictors_list, 
                                    groups, df_y_full_sample, train_set_start, today,
                                    target_start, target_end, termination_day){
  
  all_datasets <- list()
  oos_points <- c()
  all_loadings <- list()
  preselected_var_list <- list()
  
  while (today[1] != update_date(termination_day)[1] | today[2] != update_date(termination_day)[2]){    # note: updating termination_day ensures also termination_day data is included
    
    o_data_today <- window_func(data=o_data, from=train_set_start, to=today)  # original data available tomorrow
    
    results <- dimensionality_reduction(o_data=o_data_today, target=target, ylags_list=ylags_list, 
                                        other_predictors_list=other_predictors_list, groups=groups, 
                                        df_y_full_sample=df_y_full_sample, today=today,
                                        target_start=target_start, target_end=target_end)
    # store results
    
    data <- results[[1]]
    oos_point <- results[[2]] 
    loadings <- results[[3]]
    
    all_datasets[[paste(today[1], " Q", today[2], sep="")]] <- data
    oos_points <- rbind(oos_points, oos_point)
    all_loadings[[paste(today[1], " Q", today[2], sep="")]] <- loadings
    
    # Update dates
    today <- update_date(today)
    target_end <- update_date(target_end)

  }
  
  return(list(all_datasets, oos_points, all_loadings))
}