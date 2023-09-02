############################## Metric functions ################################

# For Breimann RF
variance <- function(y, X_f, beta_ridge){
  # Description:
    # This function computes the MSE given a target variable. 
  # Args:
    # y: target variable (data)
    # X_f: unused argument
    # beta_ridge: unused argument

  if (length(y) <=1) return(0)
  else{
    
    var = sum((y-mean(y))^2)*(length(y))^(-1)
    
    return(var)              
  }
}

# For LLF
MSE_LLF <- function(y, X_f, beta_ridge){
  # Description:
  # This function returns the MSE of a ridge split.
  # Args:
    # y: target data as vector
    # X_f: data for all features (as matrix)
    # beta_ridge: ridge estimate obtained in parent node
  # Output:
    # MSE: Mean squared error obtained from ridge regression
  
  # Add back intercept
  X_intercept <- cbind(1, X_f)

  # obtain MSE
  MSE <- mean((y - X_intercept%*%beta_ridge)[,1]^2)
  
  return(MSE)
}


igain <- function (y, X_f, mask, func, beta_ridge){
  # Description:
    # This function computes the information gain given a mask. 
  # Args: 
    # y: target variable (data)
    # mask: boolean vector specifiying instances falling below/above the split value
    # func: metric function
    # beta_ridge: ridge estimate obtained in parent node
  # Output: 
    # igain: information gain (AKA split criterion value)
  
  
  R1_size = sum(mask)                 # Region (child) 1 size
  R2_size = sum(!mask)                # Region (child) 2 size
  
  if (R1_size==0 | R2_size==0){       
    igain = 0
  } else {
    igain = func(y=y, X_f=X_f, beta_ridge=beta_ridge) - R1_size/(R1_size+R2_size)*func(y=y[mask], X_f=X_f[mask,, drop=FALSE], beta_ridge=beta_ridge) - R2_size/(R1_size+R2_size)*func(y=y[!mask], X_f=X_f[!mask, ,drop=FALSE], beta_ridge=beta_ridge)
  }

  return(igain)
}

max_igain_split <- function(x, y, X_f, min_leaf_size, min_node_size, func, beta_ridge){
  # Description:
    # This function obtains the igain-maximizing split conditional on the splitting conditions being fulfilled.
  # Args: 
    # x: splitting variable data (a single feature)
    # y: target variable data (as vector)
    # X_f: data for all features (as dataframe)
    # func: metric function
    # beta_ridge: ridge estimate obtained in parent node
  # Output: 
    # max_igain: max information gain value
    # max_igain_split: information gain-maximizing splitpoint
  
  
  splitpoints = c()         # splitpoints
  igain = c()               # information gain of splitpoints
  
  # Define split options (all possible splits)
  x_sorted = sort(unique(x))
  options = apply(cbind(x_sorted[-length(x_sorted)], x_sorted[-1]), 1, mean) 

  # Compute information gain for all split options
  for (option in options){
    
    mask = x < option
    
    # compute igain and store results if splitting conditions fulfilled (min_leaf_size and min_node_size)
    if (sum(mask)>=min_leaf_size & sum(!mask)>=min_leaf_size & length(mask)>=min_node_size){
      
      option_igain = igain(y=y, X_f=X_f, mask=mask, func=func, beta_ridge=beta_ridge)     # information gain of given option
      splitpoints = c(splitpoints, option)
      igain = c(igain, option_igain)
      
    }
  }
  
  if (length(igain)==0) return(c(NA, NA))         # if no split fulfilling conditions was found 
  
  else {
    max_igain = max(igain)
    max_igain_idx = which(igain==max_igain)
    if (length(max_igain_idx) > 1) max_igain_idx = sample(max_igain_idx, 1)   # if multiple splitpoints provide same igain, randomly select 1
    
    max_igain_split = splitpoints[max_igain_idx]
    return(c(max_igain, max_igain_split))
  }
}

best_split <- function(target, data, min_leaf_size, min_node_size, m_try, func, lambda_split){
  # Description:
    # This function obtains the best split variable, the best split point of that variable and the maximized information gain of that variable. 
  #Args: 
    # target: target variable name as character object
    # data: dataframe including target variable and all features
    # m_try: number of variables to be randomly selected and tried for each split
    # func: metric function
    # lambda_split: the ridge penalty parameter incorporated in the splitting rules
  # Output:
    # split_variable: split variable name (as character object)
    # split_value: split value (as float)
    # split_igain: split information gain (as float)
  
  # Obtain dataframe containing only features and vector containing only target
  X_f <- data[,!(names(data) %in% c(target))]
  y <- data[,target]
  
  # if min_node_size violated, or min_leaf_size cannot be satisfied terminate
  if (length(y) < min_node_size | (length(y)/2) < min_leaf_size){
    return("min_leaf_size and/or min_node_size violated for all splits for all variables")
  }
  else {
    
    # Random split variables selection
    var_try <- sort(sample(ncol(X_f), size=m_try, replace=FALSE))
    X_f <- X_f[,var_try]
    
    # if we are dealing with local_linear_forest
    if(!is.null(lambda_split)){
    
      ## Run ridge regression in parent
      # Define X_f as X
      X <- as.matrix(X_f)
      
      # Ridge regression w. unpenalized intercept 
      X_intercept <- cbind(1, X)
      diagonal <- diag(ncol(X_intercept))
      diagonal[1,1] <- 0
      beta_ridge <- chol2inv(chol(t(X_intercept)%*%X_intercept + lambda_split*diagonal)) %*%t(X_intercept)%*%y

      
      # Define centered and standardized X_f (if standardization is applied)
      X_f_std <- X
      
    }
    else{ #(i.e. if regression_forest is used)

      # define as null
      X_f_std <- NULL
      beta_ridge <- NULL
      
    }
    
    ## Obtain best split for each feature
    best_split_per_var = apply(X=X_f, MARGIN=2, FUN=max_igain_split, y=y, X_f=X_f_std, min_leaf_size=min_leaf_size, min_node_size=min_node_size, func=func, beta_ridge=beta_ridge)
  
    # if any variable can be split while fulfilling splitting criteria
    if (any(!is.na(best_split_per_var))){
      
      # Obtain max_igain_split results from variable with highest information gain
      split_variable_idx = which.max(best_split_per_var[1,])   
      split_variable = names(X_f)[split_variable_idx]
      split_value = best_split_per_var[[2, split_variable]]    
      split_igain = best_split_per_var[[1, split_variable]]
  
      # Return results
      return(list(split_variable, split_value, split_igain))
      
    }
    else {
      return("min_leaf_size and/or min_node_size violated for all splits for all variables")
    }
  }
}

get_prediction = function(y){
  # Description:
    # This function obtains a prediction as the mean of the input data.
  # Args: 
    # y: target variable data for a given region
  # Output:
    # prediction: prediction in the form of the mean
  
  prediction = mean(y)
  
  return(prediction)
}

get_tree <- function(target, data, min_leaf_size, min_node_size, m_try, func, lambda_split){
  # Description:
    # This function builds a regression tree using CART principles. 
  # Arguments:
    # target: target variable (as character string)
    # data: target variable and features (as dataframe)
    # min_leaf_size: minimum leaf size (as integer)
    # min_node_size: minimum size of a node for a split to be allowed
    # func: metric function / local MSE function
    # m_try: number of variables to be randomly selected and tried for each split
    # lambda_split: parameters involved in local estimation function. i.e. lambda for Local Linear Forest
  # Output:
    #  list containing
        # frame: the tree in form of a dataframe
        # nodes_data: a list of dataframes containing the data for each node
  
  # initialize dataframe that will store results
  columns <- c("left_child", "right_child", "split_variable", "split_point", "is_leaf" , "n", "prediction", "in_queue")
  frame <- as.data.frame(matrix(NA, nrow=1, ncol=length(columns)))
  colnames(frame) <- columns
        # note: "in_queue" is 0 if node has been evaluated (split or left as leaf) and 1 if node is unevaluated (in queue)
  
  # initialize list that will store results
  nodes_data <- list()
  
  # initialize root: add root to queue
  frame$in_queue[1]=1
  nodes_data[[1]] <- data
  
  # While queue is not empty (/while there are unevaluated nodes left in the tree)
  while(sum(frame$in_queue)!=0){
    
    # obtain row indicices for nodes in queue
    indices <- which(frame$in_queue!=0)
    
    # evaluate each node in queue
    for (idx in indices){
    
      # obtain best split based on node data
      node_data <- nodes_data[[idx]]
      best_split_result = best_split(target=target, data=node_data, min_leaf_size=min_leaf_size, min_node_size=min_node_size, m_try=m_try, func=func, lambda_split=lambda_split)
      
      # if stopping conditions are not violated: split
      if (best_split_result[1] != "min_leaf_size and/or min_node_size violated for all splits for all variables"){
        
        split_variable = best_split_result[[1]]
        split_value = best_split_result[[2]]
        split_igain = best_split_result[[3]]
        
        # obtain split data
        mask <- node_data[,split_variable] < split_value
        left_child_data <- node_data[mask,]
        right_child_data <- node_data[!mask,]
        
        # fill out information in row idx
        frame[idx, "left_child"] = nrow(frame)+1
        frame[idx, "right_child"] = nrow(frame)+2
        frame[idx, "split_variable"] = split_variable
        frame[idx, "split_point"] = split_value
        frame[idx, "is_leaf"] = 0
        frame[idx, "n"] = nrow(node_data)
        frame[idx, "in_queue"] = 0
        frame[idx, "prediction"] = 0
        
        # append 2 new dataframes to nodes_data: left child and right child
        nodes_data[[(nrow(frame)+1)]] <- left_child_data
        nodes_data[[(nrow(frame)+2)]] <- right_child_data
        
        # append 2 new rows to frame: left child and right child
        frame[nrow(frame) + 1,] = c(rep(NA, ncol(frame)))
        frame[nrow(frame) + 1,] = c(rep(NA, ncol(frame)))
        
        # Add left and right child row to queue
        frame[nrow(frame) -1, "in_queue"] = 1
        frame[nrow(frame), "in_queue"] = 1
        
      }
      
      # if stopping conditions are violated: terminate
      else {
        
        # obtain prediction
        prediction = get_prediction(node_data[,target])
        
        # fill out information in row idx
        frame[idx, "left_child"] = 0
        frame[idx, "right_child"] = 0
        frame[idx, "split_variable"] = 0
        frame[idx, "split_point"] = 0
        frame[idx, "is_leaf"] = 1
        frame[idx, "n"] = nrow(node_data)
        frame[idx, "in_queue"] = 0
        frame[idx, "prediction"] = prediction
  
      }
    }
  }
  
  output <- list("frame" = frame, "nodes_data" = nodes_data)
  return(output)
  
}

############################### GRF function ###################################

# Neighborhood function needed for get_grf function
get_neighborhood <- function(tree_frame, instance){
  # Description:
    # A function obtaining the leaf into which a single instance falls. 
  # Args:
    # tree_frame: tree frame in the form of a dataframe
    # instance: a (1xd)-dimensional dataframe with d as the number of features
  # Output:
    # index of the leaf into which the supplied instance falls
  
  # initialize
  node_idx = 1
  at_leaf = 0
  
  # while not at a leaf
  while (at_leaf != 1){
    
    # if instance has arrived at a leaf
    if (tree_frame[node_idx, "is_leaf"] == 1){
      
      # terminate while-loop
      at_leaf = 1
      
      return(c("leaf_instance" = node_idx))
    }
    
    # if instance falls into left child
    else if (instance[1, tree_frame[node_idx, "split_variable"]] < tree_frame[node_idx, "split_point"]){
      
      # set node_idx to left child
      node_idx = tree_frame[node_idx, "left_child"]
    }
    # if instance falls into right child child
    else if (instance[1, tree_frame[node_idx, "split_variable"]] >= tree_frame[node_idx, "split_point"]){
      
      # set node_idx to right child
      node_idx = tree_frame[node_idx, "right_child"]
    }
  }
}

# Neighborhood function needed for get_grf function
get_neighbors <- function(tree_frame, x, data){
  # Description:
    # A function obtaining the instances in data that fall into the same leaf as x.
  # Args: 
    # tree_frame: tree frame in the form of a dataframe
    # x: test point, a (1xd)-dimensional dataframe with d as the number of features
    # data: target and feature data in form of dataframe
  # Output: 
  
  # Obtain leaf into which x falls
  leaf_x = get_neighborhood(tree_frame=tree_frame, instance=x)
  
  # Remove target data from data
  data_no_target = data[,!(names(data) %in% c(target))]
  
  # Initialize empty vector for results
  in_leaf_x <- c()
  
  # Obtain instances falling into leaf_x
  for (row_instance in 1:nrow(data_no_target)){
    
    leaf_instance=get_neighborhood(tree_frame=tree_frame, instance=data_no_target[row_instance,])
    
    # If instances falls into leaf_x, store 
    if (leaf_instance == leaf_x) in_leaf_x <- c(in_leaf_x, row_instance)
  }
  
  # Obtain neighbors
  neighbors = data[in_leaf_x, ]
  
  return (neighbors)
}


get_grf <- function(target, data, B, subsample_size, min_leaf_size, min_node_size, m_try, forest_type="regression_forest", lambda_split=NULL, get_ci=FALSE, l=2){
  # Args:
    # target: target variable as character
    # data: target and feature data in form of dataframe
    # B: number of trees
    # subsample_size: size of subsample to grow trees on
    # min_leaf_size: minimum leaf size
    # min_node_size: minimum size of a node for a split to be allowed
    # m_try: number of features to randomly select as split candidates
    # forest_type: the type of forest to grow. Options are "regression_forest" and "local_linear_forest"
    # lambda_split: the ridge penalty parameter incorporated in the splitting rules
    # get_ci: whether confidence intervals should be obtained
    # l: little bag size. Only becomes effective if get_ci=TRUE
  # Output: 
    # Returns a list of containing the target, data, forest_type, B and the forest. The forest is a list of frames and nodes_data lists, thus
    # a list of get_tree outputs. 
  
  # Select appropriate metric function for given forest type
  if (forest_type == "regression_forest") func = variance
  else if (forest_type == "local_linear_forest") func = MSE_LLF
  else return("ERROR: no valid forest_type argument provided")
  
  # Ensure that lambda_split is set to NULL for regression forest
  if(forest_type == "regression_forest") lambda_split = NULL
  
  # if data has rownames, remove
  rownames(data) <- c()
  
  ## Part 1 - Regular forest predictions
  
  # initialize list to store all trees
  forest <- list()
  
  for (b in 1:B){
    
    print(paste(forest_type, "tree number", b, sep=" "))
    
    # draw subsample without replacement
    subsample_idx <- sort(sample.int(n=nrow(data), size=subsample_size, replace=FALSE))
    I <- data[subsample_idx,]
    
    # obtain tree using I
    tree <- get_tree(target=target, data=I, min_leaf_size=min_leaf_size, min_node_size=min_node_size, m_try=m_try, func=func, lambda_split=lambda_split)
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
        subsample_idx <- sort(sample(halfsample_idx, size=floor(subsample_size*l^(-1)), replace=FALSE))
        I <- data[subsample_idx, ]
        
        # obtain tree using I
        tree <- get_tree(target=target, data=I, min_leaf_size=min_leaf_size, min_node_size=min_node_size, m_try=m_try, func=func, lambda_split=lambda_split)
        bag[[l_nr]] <- tree
      }
      forest_blb[[g]] <- bag
    }
  }
  
  if (get_ci == FALSE) forest_blb = NULL
  
  # return forest and further needed information
  grf_output <- list(target, data, forest_type, B, forest, get_ci, forest_blb, l)
  
  return(grf_output)
}

predict_grf <- function(grf_output, x, lambda_est=NULL){
# Arguments
  # grf_output: the output obtained from get_grf function
  # x: test point, a (1xd)-dimensional dataframe with d as the number of features
  # lambda_est: the ridge penalty parameter incorporated in the final estimation function for local_linear_forest
# Output:
  # Forest prediction and potentially variance estimate
  
  # extract information from grf_output
  target <- grf_output[[1]]
  data <- grf_output[[2]]
  forest_type <- grf_output[[3]]
  B <- grf_output[[4]]
  forest <- grf_output[[5]]
  get_ci <- grf_output[[6]]
  forest_blb <- grf_output[[7]]
  l <- grf_output[[8]]
  
  ### Part 1 - Regular forest predictions
  
  # initialize weight vector with zeros
  alpha <- c(rep(0, nrow(data)))
  
  # obtain weights for each observation
  for (b in 1:B){
  
    # obtain leaf into which x falls 
    x_leaf <- get_neighborhood(tree_frame = forest[[b]]$frame, instance = x)
    
    # obtain its neighbors in x_leaf
    neighbors <- forest[[b]]$nodes_data[[x_leaf]]
    
    # assign/add weights to neighbors 
    alpha[as.numeric(rownames(neighbors))] <- alpha[as.numeric(rownames(neighbors))] + 1/nrow(neighbors)
  }
  
  # obtain "final/real" weights
  alpha = alpha/B
  
  # obtain prediction
  if (forest_type == "regression_forest"){
    
    # obtain theta(x) as solution to eq (2) from GRF paper
    RF_prediction = sum(alpha * data[,target]) # (Breiman RF)
  }
  
  if (forest_type == "local_linear_forest"){
    
    # Define matrix A
    A <- diag(alpha)
    
    # Define data including only features
    data_f <- data[, !(names(data) %in% c(target))]
    
    # Define target variable data
    y <- data[, target]
    
    # define (d+1)x(d+1) diagonal matrix J
    J <- diag(1, (ncol(data_f)+1))
    J[1,1] <- 0           # to not penalize intercept
    
    # Construct the centered (at test point x) matrix with intercept
    x_matrix <- matrix(as.numeric(x), ncol=ncol(x), nrow=nrow(data_f), byrow=TRUE)
    Delta <- cbind(1, (data_f-x_matrix))
    #Delta <- cbind(1, (data_f - do.call("rbind", replicate(nrow(data_f), x, simplify = FALSE))))
    
    # Convert to matrix objects
    Delta <- as.matrix(Delta)
    
    
    # Obtain local linear forest estimate
    RF_prediction <- (chol2inv(chol(t(Delta) %*% A %*% Delta + lambda_est*J)) %*% t(Delta) %*% A %*% y)[1]
  }
  
  
  ### Part 2 - Confidence intervals
  
  if (get_ci == TRUE){
    
    ## obtain xi (LLF paper notation)
    
    if (forest_type == "local_linear_forest"){
      
      # estimate M (LLF paper eq 22)
      M <- 0
      
      for (i in 1:length(alpha)){
        
        M <- M + alpha[i]*as.numeric(crossprod(Delta[i,], Delta[i,])) + lambda_est*J 
        
      }
      
      e1 <- c(1, rep(0, (nrow(M)-1)))
      xi <- t(chol2inv(chol(M))) %*% e1
    }
    
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
        
        # initialize weight vector with zeros
        alpha <- c(rep(0, nrow(data)))
        
        # obtain leaf into which x falls 
        x_leaf <- get_neighborhood(tree_frame = forest_blb[[g]][[l_nr]]$frame, instance = x)
        
        # obtain its neighbors in x_leaf
        neighbors <- forest_blb[[g]][[l_nr]]$nodes_data[[x_leaf]]
        
        # assign/add weights to neighbors 
        alpha[as.numeric(rownames(neighbors))] <- alpha[as.numeric(rownames(neighbors))] + 1/nrow(neighbors)

        # obtain prediction
        if (forest_type == "regression_forest"){
          
          RF_bag_predictions[g, l_nr] = sum(alpha * data[,target])
        }
        
        if (forest_type == "local_linear_forest"){
          
          # Define new matrix A
          A <- diag(alpha)
          
          
          # Obtain local linear forest estimate
          RF_bag_predictions[g, l_nr] <- (chol2inv(chol(t(Delta) %*% A %*% Delta + lambda_est*J)) %*% t(Delta) %*% A %*% y)[1]
        }
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
    
    if (forest_type == "regression_forest") sigma2_hat <- H_blb
    if (forest_type == "local_linear_forest") sigma2_hat <- H_blb * t(xi) %*% xi
    
    
  }
  
  if (get_ci == TRUE) return(c("prediction"=RF_prediction, "variance"=sigma2_hat))
  if (get_ci == FALSE) return(RF_prediction)
}

