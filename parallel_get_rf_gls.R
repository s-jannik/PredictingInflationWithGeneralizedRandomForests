library(foreach)
library(doParallel)
source("parallel_get_grf.R")

get_rf_gls <- function(target, data, B, subsample_size, tn, tc, m_try, get_ci=FALSE, l=2){
  
  ## Obtain Sigma2
  grf_output <- get_grf(target=target, 
                        data=data, 
                        B=B, 
                        subsample_size=subsample_size, 
                        min_leaf_size=1, 
                        min_node_size=(tc+1),
                        m_try=m_try, 
                        forest_type="regression_forest")
  
  # get in-sample predictions
  is_pred <- c()
  for (obs in 1:nrow(data)){
    
    # test point
    x0 <- data[obs, !(names(data) %in% c(target))]
    
    # prediction
    is_pred[obs] <- predict_grf(grf_output=grf_output, x=x0)
  }
  
  is_errors <- data[,target] - is_pred
  is_errors <- is_errors
  
  # Get Sigma
  ar_results <- ar(is_errors, aic=FALSE, order.max=1)
  phi <- ar_results$ar
  
  print(paste("phi", phi, sep=" "))
  
  phi_powers <- c()
  for (idx in 1:nrow(data)){
    
    phi_powers <- c(phi_powers, phi^(idx-1))
    
  }
  
  Sigma<-c()
  for (idx in 1:length(phi_powers)){
    
    lower_part <- phi_powers[1:(length(phi_powers)-idx+1)]
    
    Sigma <- cbind(Sigma, c(rep(0, (idx-1)), lower_part))
    
  }
  Sigma <- Sigma+t(Sigma)-diag(nrow(Sigma))
  
  
  # Calculate Sigma^-(1/2)
  Sigma_inv2 <- Sigma
  Sigma_inv2 <- Matpow(Sigma_inv2, numer=(-0.5))

  
  ########################## start of parallel computing #######################
  
  # Setup of parallel backend
  n.cores <- parallel::detectCores() - 1
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  ## Part 1 - Regular forest predictions 
  
  forest <- foreach(b = 1:B) %dopar% {          # forest will be list of B trees
    
    source("get_rf_gls.R")

    # print(paste("RF-GLS tree number", b, sep=" "))
    
    # Randomly draw subsample
    subsample_indicator <- c(rep(0, nrow(data)))
    to_be_included <- sample.int(n=nrow(data), size=subsample_size, replace=FALSE)
    subsample_indicator[to_be_included] <- 1
    
    # Define subsample selection matrix and Q
    P <- diag(subsample_indicator)
    Q <- t(Sigma_inv2) %*% t(P) %*% P %*% Sigma_inv2
    
    # Obtain tree
    get_DART_tree(target=target, data=data, subsample_indicator=subsample_indicator, tn=tn, tc=tc, m_try=m_try, Q=Q)
    
  }
  
  ## Part 2 - Confidence intervals 
  
  if (get_ci == TRUE){
    
    # initialize list to store all little bags of trees
    forest_blb <- list()              # blb = bootstrap of little bags
    
    forest_blb <- foreach (g = 1:(B/l)) %dopar% {
      
      source("get_rf_gls.R")
      
      #print(paste("bag nr. ", g, " of ", (B/l), sep=""))
      
      # draw random half-sample H_g
      halfsample_idx <- sort(sample.int(n=nrow(data), size=floor(nrow(data)*l^(-1)), replace=FALSE))
      
      # initialize list to store all trees in the g-th bag
      bag <- list()
      
      for (l_nr in 1:l){
        
        # draw subsample from half-sample H_g
        subsample_indicator <- c(rep(0, nrow(data)))
        to_be_included <- sample(halfsample_idx, size=floor(subsample_size*l^(-1)), replace=FALSE)
        subsample_indicator[to_be_included] <- 1
        
        # Define subsample selection matrix and Q
        P <- diag(subsample_indicator)
        Q <- t(Sigma_inv2) %*% t(P) %*% P %*% Sigma_inv2
        
        # Obtain tree
        tree <- get_DART_tree(target=target, data=data, subsample_indicator=subsample_indicator, tn=tn, tc=tc, m_try=m_try, Q=Q)
        bag[[l_nr]] <- tree
        
      }
      bag
    }
  }
  if (get_ci == FALSE){
    
    forest_blb <- NULL
    
  }
  
  # stop cluster again 
  parallel::stopCluster(cl = my.cluster)
  
  ########################## end of parallel computing #########################
  
  # return forest and further needed information
  rf_gls_output <- list(B, forest, get_ci, forest_blb, l)
  return(rf_gls_output)
}