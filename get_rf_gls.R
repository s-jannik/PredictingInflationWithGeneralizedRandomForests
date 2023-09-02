source("get_grf.R")
library(powerplus)
library(Matrix)


# igain function
igain_GLS <- function(y, Z, Z0, pre_split_SSE, Q){
  # Description:
    # This function performs GLS estimation and evaluates eq (5) of the RF-GLS paper. 
  # Output: 
    # The information gain obtained from the given split. More specifically, the reduction in total MSE from performing the given split. 
  
  term <- t(Z) %*% Q  # evaluated first for higher speed
  beta_GLS_Z <- Matrix::chol2inv(chol(term %*% Z)) %*% (term %*% y)
  
  # post-split total sum of squared errors
  errors <- y - Z %*% beta_GLS_Z
  post_split_SSE <- t(errors) %*% Q %*% errors
  
  # obtain information gain 
  igain <- 1/length(y) * (pre_split_SSE - post_split_SSE)
  
  return(as.numeric(igain))
}

max_igain_split_GLS <- function(y, data, C_K, split_var_name, Q){
  # Description:
    # This function obtains the optimal split point for 1 variable. 
  # Details on arguments: 
    # C_K: node to be split 
    # Q: covariance-related matrix
    # split_var_name: splitting variable name as character
  
  ## obtain Z0 and pre-split total SSE (sum of squared errors): 
  Z0 <- C_K
  
  # pre-split beta
  term0 <- t(Z0) %*% Q
  beta_GLS_Z0 <- Matrix::chol2inv(chol(term0 %*% Z0)) %*% (term0 %*% y)

  # pre-split total SSE
  errors <- y - Z0 %*% beta_GLS_Z0
  pre_split_SSE <- t(errors) %*% Q %*% errors
  
  
  ## Prepare elements necessary to obtain Z: 
  
  # obtain instances in node K
  K <- ncol(Z0)
  x <- data[as.logical(Z0[,K]), split_var_name]       # splitting variable data
  
  # Define split options (all possible splits)
  x_sorted = sort(unique(x))
  options = apply(cbind(x_sorted[-length(x_sorted)], x_sorted[-1]), 1, mean) 
  
  # empty vector to store results
  igain <- c()
  
  
  ## Compute information gain for every possible split
  for (option in options){
    
    # Define left and right child node 
    C_K_L <- Z0[,K] * as.numeric(data[, split_var_name] < option)
    C_K_R <- Z0[,K] * as.numeric(!(data[, split_var_name] < option))
                # ensures that only the instances that are actually in node K are split into left and right child
    
    # Define Z
    Z <- cbind(Z0[,-K], C_K_L, C_K_R) 
    
    # Compute information gain
    igain <- cbind(igain, igain_GLS(y=y, Z=Z, Z0=Z0, pre_split_SSE=pre_split_SSE, Q=Q))
  }
  
  ## Obtain max igain and optimal split point 
  max_igain <- max(igain)
  max_igain_idx <- which(igain==max_igain)
  if (length(max_igain_idx) > 1) max_igain_idx = sample(max_igain_idx, 1)   # if multiple splitpoints provide same igain, randomly select 1
  
  max_igain_split = options[max_igain_idx]
  
  return(c("max_igain"=max_igain, "max_igain_split"=max_igain_split))
}


best_split_DART <- function(y, data, target, C_K, Q, m_try){
  # Description:
    # This function obtains the optimal split variable and split point combination.
  
  feature_names <- colnames(data[,!(names(data) %in% c(target))])
  
  # Random split variables selection
  var_try <- sample(feature_names, size=m_try, replace=FALSE)
  
  # empty vector for data
  best_split_per_var <- c()
  
  for(split_var_name in var_try){
    
    best_split_per_var <- rbind(best_split_per_var, max_igain_split_GLS(y=y, data=data, C_K=C_K, split_var_name=split_var_name, Q=Q))
    
  }
  

  best_split_var_idx <- which.max(best_split_per_var[,"max_igain"])
  split_variable <- var_try[best_split_var_idx]
  split_value <- best_split_per_var[best_split_var_idx,"max_igain_split"]

  return(list("split_variable" = split_variable, "split_value"=split_value))
}


get_DART_tree <- function(target, data, subsample_indicator, tn, tc, m_try, Q){
  # Description:
    # This function builds a regression tree using DART principles. 
        # The elements of a tree frame are described in line 132 f.
  # Details on arguments:
    # All arguments follow exactly the definition in the RF-GLS paper (Saha et. al., 2021). 
  # Output:
    # list containing
      # frame: the tree in form of a dataframe
      # C_frak: matrix indicating the complete set of observations in each leaf
  
  
  ################### Pre-definitions and initializations ######################
  
  # Define target variable data as y
  y <- data[, target]
  
  # initialize dataframe that will store tree
  columns <- c("left_child", "right_child", "split_variable", "split_point", 
               "is_leaf" , "n", "prediction", "in_queue", "k", "col_next_C_frak")
  frame <- as.data.frame(matrix(NA, nrow=1, ncol=length(columns)))
  colnames(frame) <- columns
    # Descriptions:
      # left_child & right_child: row number of the node's left at right child. If the node is a leaf, both are 0.
      # split_variable: name of the split variable 
      # split_point: the value at which the split is performed
      # is_leaf: 1 if yes, 0 if no. 
      # n: number of elements in the node
      # prediction: prediction value, only nonzero if node is a leaf
      # k: level at which the node occurs in the tree
      # col_next_C_frak: the column of the node in C_frak at level k+1
  
  
  # initialize root
  frame[1, "n"] <- nrow(data)
  frame[1, "prediction"] <- mean(y)
  frame[1, "in_queue"] <- 1
  frame[1, "k"] <- 1
  frame[1, "col_next_C_frak"] <- 1
  
  # initializations concerning tree algorithm
  k <- 1              # level
  g <- 1              # vector of total number of nodes up until each level. Element 1 in g is level 1, element 2 is level 2, and so on
  
  C_frak <- list(matrix(subsample_indicator, nrow=nrow(data), ncol=1))          # initialize root such that it contains only the subsample
  
  
  
  ################## Recursive binary partitioning procedure ###################
  
  # While max_num_nodes and min_node_size conditions are fulfilled
  while ( (g[k] < tn) & (any(apply(C_frak[[k]], MARGIN = 2, FUN=sum) > tc) )){
    
    # update
    k <- k + 1
    
    # initialize C_frak[[k]]
    C_frak <- append(C_frak, list(c()))
    g[k] <- 0
    
    for (l1 in 1:g[k-1]){
      
      C_l1_ksub1 <- C_frak[[k-1]][,l1]
      
      # if stopping conditions violated
      if (sum(C_l1_ksub1) <= tc | g[k] >= tn) {
        
        # don't split, just add to the new complete set of nodes 
        C_frak[[k]] <- cbind(C_frak[[k]], C_l1_ksub1)
        g[k] <- g[k] + 1
      
        
        ## Save new information in frame
        # identify leaf among all nodes up until level (k-1)
        leaf_row <- which(frame[frame$k != k, ]$col_next_C_frak == l1)
        
        # remove from queue and store column number for next C_frak
        frame[leaf_row, "in_queue"] = 0 
        frame[leaf_row, "col_next_C_frak"] = ncol(C_frak[[k]])
      
  
        # Re-evaluate stopping condition (max_num_nodes). Break if violated.
        if (g[k] >= tn) {
          
          print("MANUAL BREAK")
          
          break  
        }
        
      }
      else {
        
        # obtain membership matrix with C_l1_ksub1 (node to be split) pushed to the last column
        C_K <- cbind(C_frak[[k-1]][,-l1], C_frak[[k-1]][,l1])
        
        # Perform DART split
        best_split_result <- best_split_DART(y=y, data=data, target=target, C_K=C_K, Q=Q, m_try=m_try)
        split_variable <- best_split_result[[1]]
        split_value <- best_split_result[[2]]
        
        # Define data in left and right child node
        left_child_mask <- as.numeric(data[, split_variable] < split_value)
        right_child_mask <- as.numeric(!(data[, split_variable] < split_value))
        C_l1_1_k <- C_l1_ksub1*left_child_mask
        C_l1_2_k <- C_l1_ksub1*right_child_mask
        
        # Append children to complete set of nodes in level k
        C_frak[[k]] <- cbind(C_frak[[k]], C_l1_1_k, C_l1_2_k)
        
        # update g[k]
        g[k] <- g[k] + 2
        
        
        ## Save new information in frame 
        
        # Identify node-to-be-split (node_tbs) (as the oldest node in the queue)
        node_tbs_indicator <- which(frame[frame$k==(k-1),"in_queue"]==1)[1]
        node_tbs <- as.numeric(rownames(frame[frame$k==(k-1),])[node_tbs_indicator])
        
        ## Fill out information in node_tbs
        frame[node_tbs, "left_child"] <- nrow(frame)+1
        frame[node_tbs, "right_child"] <- nrow(frame)+2
        frame[node_tbs, "split_variable"] <- split_variable
        frame[node_tbs, "split_point"] <- split_value
        frame[node_tbs, "is_leaf"] <- 0
        frame[node_tbs, "prediction"] <- 0
        frame[node_tbs, "in_queue"] <- 0                                          # define as evaluated
        
        # Store column number of node_tbs in the next C_frak
        frame[node_tbs, "col_next_C_frak"] <- 0
        
        
        
        ##### Add child nodes and add to queue
        
        # Append 2 new rows to frame: left child and right child
        frame[(nrow(frame) + 1),] = c(rep(NA, ncol(frame)))
        frame[(nrow(frame) + 1),] = c(rep(NA, ncol(frame)))
        
        # Add left and right child row to queue
        frame[(nrow(frame) -1), "in_queue"] = 1
        frame[nrow(frame), "in_queue"] = 1
        
        
        
        # Fill out remaining information (under the assumption that the two children are leaves)
        
        # LEFT CHILD
        frame[(nrow(frame) -1), "left_child"] <- 0                                # has no children (yet)
        frame[(nrow(frame) -1), "right_child"] <- 0
        frame[(nrow(frame) -1), "split_variable"] <- 0                            # has not been split (yet)
        frame[(nrow(frame) -1), "split_point"] <- 0                               # has not been split (yet)
        frame[(nrow(frame) -1), "is_leaf"] <- 1                                   # is assumed leaf (until proven splittable)
        frame[(nrow(frame) -1), "n"] <- sum(C_l1_1_k)
        frame[(nrow(frame) -1), "prediction"] <- mean(y[as.logical(C_l1_1_k)])
        frame[(nrow(frame) -1), "k"] <- k
        frame[(nrow(frame) -1), "col_next_C_frak"] <- ncol(C_frak[[k]])-1
        
        # RIGHT CHILD
        frame[nrow(frame), "left_child"] <- 0                                     # has no children (yet)
        frame[nrow(frame), "right_child"] <- 0
        frame[nrow(frame), "split_variable"] <- 0                                 # has not been split (yet)
        frame[nrow(frame), "split_point"] <- 0                                    # has not been split (yet)
        frame[nrow(frame), "is_leaf"] <- 1                                        # is assumed leaf (until proven splittable)
        frame[nrow(frame), "n"] <- sum(C_l1_2_k)
        frame[nrow(frame), "prediction"] <- mean(y[as.logical(C_l1_2_k)])
        frame[nrow(frame), "k"] <- k
        frame[nrow(frame), "col_next_C_frak"] <- ncol(C_frak[[k]])
      }
    }
  }
  
  # Remove any remaining leafs from queue
  frame[frame$in_queue==1, "in_queue"] <- 0
  
  # Return only final C_frak
  C_frak <- C_frak[[k]]
  
  return(list("frame"=frame, "C_frak"=C_frak))
}



get_rf_gls <- function(target, data, B, subsample_size, tn, tc, m_try, get_ci=FALSE, l=2){
  # Description:
    # This function runs the full RF-GLS algorithm 
  # Output: 
    # Returns a list of containing the forest, B and information needed to produce confidence intervals. The forest is a list of tree frames.
  
  
  ## Obtain Sigma2 by running a standard RF
  grf_output <- get_grf(target=target, 
                        data=data, 
                        B=B, 
                        subsample_size=subsample_size, 
                        min_leaf_size=1, 
                        min_node_size=(tc+1),
                        m_try=m_try, 
                        forest_type="regression_forest")
  
  # get in-sample (is) predictions
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
  
  Sigma_inv2 <- Sigma
  Sigma_inv2 <- Matpow(Sigma_inv2, numer=(-0.5))
  
  
  ## Part 1 - Regular forest predictions
  
  # initialize list to store all trees
  forest <- list()
  
  for (b in 1:B){
    
    print(paste("RF-GLS tree number", b, sep=" "))
    
    # Randomly draw subsample
    subsample_indicator <- c(rep(0, nrow(data)))
    to_be_included <- sample.int(n=nrow(data), size=subsample_size, replace=FALSE)
    subsample_indicator[to_be_included] <- 1
    
    # Define subsample selection matrix and Q
    P <- diag(subsample_indicator)
    Q <- t(Sigma_inv2) %*% t(P) %*% P %*% Sigma_inv2

    
    # Obtain tree
    tree <- get_DART_tree(target=target, data=data, subsample_indicator=subsample_indicator, tn=tn, tc=tc, m_try=m_try, Q=Q)
    forest[[b]] <- tree

  }
  
  ## Part 2 - Confidence intervals 
  
  if (get_ci == TRUE){
    
    # initialize list to store all little bags of trees
    forest_blb <- list()              # blb = bootstrap of little bags
    
    for (g in 1:(B/l)) {
      
      print(paste("bag nr. ", g, " of ", (B/l), sep=""))
      
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
      forest_blb[[g]] <- bag
    }
  }
  if(get_ci == FALSE){
    
    forest_blb <- NULL
    
  }
  
  
  # return forest and further needed information
  rf_gls_output <- list(B, forest, get_ci, forest_blb, l)
  return(rf_gls_output)
}

predict_rf_gls <- function(rf_gls_output, x){
  # Description:
    # This function takes the output from the get_rf_gls function and returns a prediction (potentially also a variance estimate).
  # Arguments
    # rf_gls_output: the output obtained from get_rf_gls function
    # x: test point, a (1xd)-dimensional dataframe with d as the number of features
  
  # extract information from  rf_gls_output
  B <- rf_gls_output[[1]]
  forest <- rf_gls_output[[2]]
  get_ci <- rf_gls_output[[3]]
  forest_blb <- rf_gls_output[[4]]
  l <- rf_gls_output[[5]]
  
  ### Part 1 - Regular forest predictions
  
  # initialize vector to hold tree predictions 
  tree_prediction <- c()
  
  # obtain prediction from each tree 
  for (b in 1:B){
    
    # obtain leaf into which x falls 
    x_leaf <- get_neighborhood(tree_frame = forest[[b]]$frame, instance = x)
    
    # obtain prediction of tree b
    tree_prediction[b] <- forest[[b]]$frame[x_leaf, "prediction"]
  }
  RF_prediction <- mean(tree_prediction)
  
  
  ### Part 2 - Confidence intervals 
  if (get_ci == TRUE){
   
    ## obtain H_blb (bootstrap of little bags estimator)
    
    RF_bag_predictions <- matrix(NA, nrow=(B/l), ncol=l)
    colnames <- c()
    for (l_nr in 1:l){
      colnames <- c(colnames, paste("pred_bag_", l_nr))
    }
    colnames(RF_bag_predictions) <- colnames
    
    # obtain prediction for each bag g
    for (g in 1:(B/l)){
      
      # get tree prediction from each bag 
      for (l_nr in 1:l){
        
        # obtain leaf into which x falls 
        x_leaf <- get_neighborhood(tree_frame = forest_blb[[g]][[l_nr]]$frame, instance = x)
 
        # obtain prediction
        RF_bag_predictions[g, l_nr] = forest_blb[[g]][[l_nr]]$frame[x_leaf, "prediction"]
        
      }
    }
    
    between_group_var <- mean((apply(RF_bag_predictions, MARGIN=1, FUN=mean) - RF_prediction)^2)
    within_group_var_terms <-  c()
    for (g in 1:(B/l)){
      
      support_sum <- c()
      for (l_nr in 1:l){
        support_sum <- (RF_bag_predictions[g, l_nr] - mean(RF_bag_predictions[g,]))^2
      }
      within_group_var_terms[g] <- mean(support_sum)
    }
    within_group_var <- mean(within_group_var_terms)
    
    H_blb <- between_group_var -1/(l-1)*within_group_var
    
    sigma2_hat <- H_blb
  }
  
  if (get_ci == TRUE) return(c("prediction"=RF_prediction, "variance"=sigma2_hat))
  if (get_ci == FALSE) return(RF_prediction)
}





