## Author: Jannik Steenbergen

## libraries
library(readxl)
library(openxlsx)
library(urca)

preprocess_data <- function(data){
  # Remove a set of variables from the dataset 
  data <- data[,!colnames(data) %in% c("BONFC", "ESI", "CSTRCONF", "STOASSESS", "BCI", "RETCONF", "SERVCONF")]
  
  # create quarterly sequence
  q_seq <- seq(as.Date("1990-01-01"), by="quarter", length.out = 132)
  
  # convert quarterly sequence to preferred formatting
  q_seq <- paste(as.POSIXlt(q_seq)$year+1900, quarters(q_seq))
  
  # Remove column with dates and add rownames
  data <- data[,-1]
  row.names(data) <- q_seq
  
  # Put inflation (target) as first variable
  data <- cbind("HICP"=data[,"HICP"], data[,!colnames(data) %in% c("HICP")])
}


lag_func <- function(y, lag_order){
  
  y <- c(rep(NA, lag_order), y[1:(length(y)-lag_order)])
  
  return(y)
}


get_prediction_data <- function(data, max_lag, h){
  # description:
    # returns prediction data
  # inputs:
    # data: data as dataframe
    # max_lag: maximum lag order
    # h: forecasting horizon  
  
  pred_data <- data.frame(row.names=rownames(data))
  
  # HICP
  pred_data[,"HICP"] <-data[,"HICP"]
  
  for (v in colnames(data)){
    
    for (lag in 1:max_lag){
      
      v_name <- paste(v, "_l", lag, sep="")
      
      pred_data[,v_name] <- lag_func(y=data[,v], lag_order=(lag+(h-1)))
      
    }
  }
  
  return(pred_data)
}

transform_data <- function(data, take_log, n_diff, year_on_year, as_ts_object){
  # Arguments:
    # data: dataframe 
    # take_log: TRUE or FALSE
    # n_diff: number of times do do differencing 
    # year_on_year: TRUE or FALSE. Whether to do year_on_year differencing on the first round of differencing
    # as_ts_object: if TRUE, returns data as ts object
  # Output: data transformed in accordance with inputs
  
  
  if (take_log == TRUE){
    
    # indices of variables with negative or zero values
    indices_all <- 1:ncol(data)
    indices_neg <- indices_all[apply(apply(data, MARGIN=2, sign)==-1, MARGIN=2, any, na.rm=TRUE)]
    indices_zero <- indices_all[apply(apply(data, MARGIN=2, sign)==0, MARGIN=2, any, na.rm=TRUE)]
    indices_neg_or_zero <- sort(unique(c(indices_neg, indices_zero)))
    
    # remove variables with negative values
    data[,indices_neg_or_zero] <- NA
    
    # take log 
    data <- log(data)
    
  }
  
  if (n_diff>0){
    if (year_on_year == TRUE){
      data <- diff(as.matrix(data), lag=4, differences=1)
      if (n_diff == 2){
        data <- diff(data, lag=1, differences=1)
      }
    } else{
      data <- diff(as.matrix(data), differences=n_diff)
    }
  }
  
  # return as ts object or dataframe
  if (as_ts_object == TRUE){
    
    # obtain start year and start quarter
    start <- as.numeric(strsplit(gsub("\\D", " ", rownames(data)[1]), "  ")[[1]])
    start_year <- start[1]
    start_quarter <- start[2]
    
    # define as ts
    ts_data <- ts(data, frequency=4, start=c(start_year,start_quarter))
    
    return(ts_data)
  }
  
  else{
    data <- as.data.frame(data)
    return(data)
  }
}

transform_data_tcode <- function(data, tcode, use_which_data="all_available", as_ts_object){
  # Arguments:
    # data: dataframe with qseq as rownames 
    # tcode: a named vector of tcodes for each variable
    # use_which_data: either "all_available" or "largest_balanced"
    # "all_available" uses the largest dataset available for each variable. The size differs across variables
    # "largest_balanced" uses the length of the shortest TRANSFORMED variable for all TRANSFORMED variables.
      # this means that inflation uses more observations than other series since it is year-on-year difference, 
      # but transformed inflation is just as long as the shortest transformed variable. 
  # Output
    # transformed data as ts object 
  
  
  # transform data to ts object
  start <- as.numeric(strsplit(gsub("\\D", " ", rownames(data)[1]), "  ")[[1]])
  start_year <- start[1]
  start_quarter <- start[2]
  
  ts_data <- ts(data, frequency=4, start=c(start_year,start_quarter))
  
  # initialize empty ts object for transformed data
  ts_t_data <- ts(frequency=4, start=start(ts_data))
  
  for (var in colnames(ts_data)){
    
    if (tcode[var] == 1){
      
      t_var <- ts_data[,var]                # transformed variable
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 2){
      
      t_var <- diff(ts_data[,var])
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 3){
      
      t_var <- diff(ts_data[,var], differences=2)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 4){
      
      t_var <- log(ts_data[,var])
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 5){
      
      t_var <- diff(log(ts_data[,var]))
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 6){
      
      t_var <- diff(log(ts_data[,var]), differences=2)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 7){
      
      t_var <- diff(ts_data[,var]/lag(ts_data[,var]) - 1)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 8){
      
      t_var <- diff(log(ts_data[,var]))*400
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 9){
      
      t_var <- diff(log(ts_data[,var]), lag=4)*100
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    if (tcode[var] == 10){
      
      t_var <- diff(diff(log(ts_data[,var]), lag=4)*100)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
  }
  
  # remove support column
  ts_t_data <- ts_t_data[,-1]
  
  colnames(ts_t_data) <- colnames(ts_data)
  
  # return
  if(as_ts_object == TRUE) {
    return(ts_t_data)
  }
  
  if(as_ts_object == FALSE){
    
    t_data <- as.data.frame(ts_t_data)
    rownames(t_data) <- rownames(data)
    return(t_data)
  }
  
  
}

plot_data <- function(data, take_log, n_diff, year_on_year, purpose){
  # Arguments: 
  # see function transform_data
  
  # obtain transformed data 
  ts_data <- transform_data(data=data, take_log=take_log, n_diff=n_diff, year_on_year=year_on_year, as_ts_object=TRUE)
  
  ## plot
  
  # title labels
  titles1 <- c("", "log ")
  titles2 <- c("level", "1diff", "2diff")
  titles3 <- c("", "yoy ")
  
  #par(mfrow=c(3,2))
  par(mfrow=c(3,3))
  
  # only plot the log-transformed variables
  indices <- 1:ncol(ts_data)
  indices <- indices[apply(!is.na(ts_data), MARGIN=2, FUN=any)]
  
  for (i in indices){
    
    if (purpose == "analysis"){
      ts.plot(ts_data[,i], gpars=list(ylab="", xlab=""), 
              main=paste(colnames(ts_data)[i], " (", titles1[take_log+1], titles3[year_on_year+1], titles2[n_diff+1], ")", sep=""))
    }
    else if (purpose == "presentation"){
      ts.plot(ts_data[,i], gpars=list(ylab="", xlab=""), 
              main=colnames(ts_data)[i])
    }
  }
}

plot_data_tcode <- function(data, tcode){
  # Arguments: 
  # see function transform_data_tcode
  
  # obtain transformed data 
  ts_data <- transform_data_tcode(data=data, tcode=tcode, as_ts_object=TRUE)
  
  ## plot
  
  # title labels
  titles1 <- c("", "log ")
  titles2 <- c("level", "1diff", "2diff")
  titles3 <- c("", "yoy ")
  
  #par(mfrow=c(3,2))
  par(mfrow=c(3,3))
  
  # only plot the log-transformed variables
  indices <- 1:ncol(ts_data)
  indices <- indices[apply(!is.na(ts_data), MARGIN=2, FUN=any)]
  
  for (i in indices){
    
    ts.plot(ts_data[,i], gpars=list(ylab="", xlab=""), 
            main=colnames(ts_data)[i])
  }
  
}

unitroot_test_func <- function(data, take_log, n_diff, test_type, lshort=TRUE){
  # Arguments:
    # data: data as dataframe
    # take_log: takes TRUE or FALSE
    # n_diff: number of times to take differnces 
    # test_type: either "adf" or "kpss"
    # lshort: if test_type=="kpss", lshort determines lag truncation parameter
  # Output: 
    # ADF test statistics
  
  # define deterministic components / null hypothesis
  
  if (test_type == "adf") det_comps <- c("none", "drift", "trend")
  if (test_type == "kpss") det_comps <- c("Level", "Trend")
  
  # initialize dataframe to store results 
  if (test_type == "adf") res_colnames <- c("ADF(n)", "ADF(c)", "ADF(ct)", "n_lagged_diff")
  if (test_type == "kpss") res_colnames <- c("KPSS(Level)", "KPSS(Trend)", "pvalue(Level)", "pvalue(Trend)", "l_trunc")
  
  res <- matrix(NA, nrow=ncol(data), ncol=length(res_colnames))
  colnames(res) <- res_colnames
  rownames(res) <- colnames(data)
  
  # obtain transformed data
  data <- transform_data(data=data, take_log=take_log, n_diff=n_diff, year_on_year=FALSE, as_ts_object=FALSE)
  
  # in case of log-transformation, test only the variables for which the log is defined for all values
  indices <- 1:ncol(data)
  indices <- indices[apply(!is.na(data), MARGIN=2, FUN=any)]
  
  # for each variable
  for (i in indices){
    
    # define data
    y <- na.omit(data[,i])
    
    if (!any(is.nan(y))){
      
      # for each det. comp. combination / null hypothesis
      for (j in 1:length(det_comps)){
        
        # determine n_lagged_diff
        if (test_type == "adf") res[i, "n_lagged_diff"] <- ceiling(length(y)^(1/3))
        
        # perform unit root test
        if (test_type == "adf"){
          ur.df_object <- ur.df(y=y, type=det_comps[j], lags=res[i, "n_lagged_diff"], selectlags="AIC")
          res[i, j]  <- attr(ur.df_object, "teststat")[1]
        }
        if (test_type == "kpss"){
          
          suppressWarnings(kpss_results <- tseries::kpss.test(x=y, null=det_comps[j], lshort=lshort)) # supresses warning regarding p-values. Well-known problem that the package has trouble reporting them. Doesn't affect test statistics.
          
          res[i, j] <- as.numeric(kpss_results[[1]])
          res[i, (j+2)] <- as.numeric(kpss_results[[3]])
          res[i, "l_trunc"] <- as.numeric(kpss_results[[2]])
          
        }
      }
    }
  }
  return(res)
}

get_sheet_list <- function(file){
  # Description:
    # function importing xlsx file with multiple sheets
  
  sheet_list <- list()
  sheet_names <- openxlsx::getSheetNames(file)
  for (sn in sheet_names){
    sheet_list[[sn]] <- openxlsx::read.xlsx(file, sheet=sn, colNames=TRUE, rowNames=TRUE)
  }
  return(sheet_list)
}

adf_kpss_tcode_func <- function(adf_results, kpss_results, sign_level, kpss_test_type){
  # Arguments:
    # adf_results: list of outputs from unitroot_test_func for each possible transformation
    # kpss_results: list of outputs from unitroot_test_func for each possible transformation
    # sign_level: significance level as float 
  # kpss_test_type: "Level" or "Trend"
  
  # Import chosen deterministic components
  det_comps <- get_sheet_list(file = "unitroot_tests_det_comps.xlsx")
  
  # Import log transformation decision
  log_decision <- get_sheet_list(file = "unitroot_tests_log_decision.xlsx")[[1]]
  
  # Define list of variables
  variables_list <- rownames(log_decision)
  
  # Critical values
  cv_adf <- cbind(c(-2.60, -2.24, -1.95), c(-3.51, -3.17, -2.89), c(-4.04, -3.73, -3.45))
  colnames(cv_adf) <- c("n", "c", "ct")
  rownames(cv_adf) <- c("0.01", "0.05", "0.1")
  
  cv_kpss <- cbind(c(0.739, 0.463, 0.347), c(0.216, 0.146, 0.119))
  colnames(cv_kpss) <- c("Level", "Trend")
  rownames(cv_kpss) <- c("0.01", "0.05", "0.1")
  
  # initialize matrix to store tcodes
  tcode_res <- matrix(NA, nrow=length(variables_list), ncol=4)
  rownames(tcode_res) <- variables_list
  colnames(tcode_res) <- c("tcode", "adf_rejected", "kpss_rejected", "H_alternative")  # H_alternative is alternative hyp. (from adf) which we reject the null in favor of
  tcode_res <- as.data.frame(tcode_res)
  
  # run adf-kpss-search algorithm for each variable
  for (var in variables_list){
    
    # obtain log decision
    log_dec_var <- log_decision[var, ]
    
    # extract ADF and KPSS test results for log decision
    if (log_dec_var == 0){
      
      adf_results_var <- adf_results[1:3]
      kpss_results_var <- kpss_results [1:3]
    }
    if (log_dec_var == 1){
      
      adf_results_var <- adf_results[4:6]
      kpss_results_var <- kpss_results [4:6]
      
    }
    
    # initialize test order
    test_order <- 0         # (initialize at test for I(0))
    
    # initialize rejection decisions
    reject_adf = "initialization"
    reject_kpss = "initialization"
    
    while (reject_adf != TRUE & reject_kpss != FALSE){
      
      ## ADF test for I(0)
      
      # det. comps
      if (log_dec_var == 0) det_comps_var <- det_comps[[test_order+1]][var, "no_log"]
      if (log_dec_var == 1) det_comps_var <- det_comps[[test_order+1]][var, "log"]
      
      # test stat
      if (det_comps_var == "n") adf_stat <- adf_results_var[[test_order+1]][var, "ADF(n)"]
      if (det_comps_var == "c") adf_stat <- adf_results_var[[test_order+1]][var, "ADF(c)"]
      if (det_comps_var == "ct") adf_stat <- adf_results_var[[test_order+1]][var, "ADF(ct)"]
      
      # rejection dec.
      reject_adf <- adf_stat < cv_adf[as.character(sign_level), det_comps_var]
      
      ## KPSS test for I(0)
      
      # test stat
      if (kpss_test_type == "Level") kpss_stat <- kpss_results_var[[test_order+1]][var, "KPSS(Level)"]
      if (kpss_test_type == "Trend") kpss_stat <- kpss_results_var[[test_order+1]][var, "KPSS(Trend)"]
      
      # rejection dec.
      reject_kpss <- kpss_stat > cv_kpss[as.character(sign_level), kpss_test_type]
      
      
      ## if procedure has reached termination
      if (reject_adf == TRUE | reject_kpss == FALSE){
        
        # Rejection decision
        tcode_res[var, "adf_rejected"] <- reject_adf
        tcode_res[var, "kpss_rejected"] <- reject_kpss
        
        # tcode
        if (log_dec_var == 0){
          
          # tcode
          if (test_order == 0) tcode_res[var, "tcode"] <- 1
          if (test_order == 1) tcode_res[var, "tcode"] <- 2
          if (test_order == 2) tcode_res[var, "tcode"] <- 3
          
        }
        
        if (log_dec_var == 1){
          
          if (test_order == 0) tcode_res[var, "tcode"] <- 4
          if (test_order == 1) tcode_res[var, "tcode"] <- 5
          if (test_order == 2) tcode_res[var, "tcode"] <- 6
          
        }
        
        # H_alternative
        if (reject_adf == TRUE){
          
          if (det_comps_var == "n") tcode_res[var, "H_alternative"] <- "stationary with zero mean"
          if (det_comps_var == "c") tcode_res[var, "H_alternative"] <- "stationary with nonzero mean"
          if (det_comps_var == "ct") tcode_res[var, "H_alternative"] <- "trend-stationary"
          
        }
        if (reject_adf == FALSE){
          tcode_res[var, "H_alternative"] <- "ADF did not reject. Conclusion comes from KPSS."
        }
      }
      else{
        
        # update test order
        test_order <- test_order + 1
        
      }
    }
  }
  return(tcode_res)
}

