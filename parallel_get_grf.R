library(foreach)
library(doParallel)

get_grf <- function(data, target, B, subsample_size, min_leaf_size, min_node_size, m_try, forest_type="regression_forest", lambda_split=NULL, get_ci=FALSE, l=2){
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
  # Output: 
    # Returns a list of containing the target, data, forest_type, B and the forest. The forest is a list of list of frames and nodes_data lists, thus
    # a list of get_tree outputs. 
  
  # Select appropriate metric function for given forest type
  if (forest_type == "regression_forest") func = variance
  else if (forest_type == "local_linear_forest") func = MSE_LLF
  else return("ERROR: no valid forest_type argument provided")
  
  # Ensure that lambda_split is set to NULL for regression forest
  if(forest_type == "regression_forest") lambda_split = NULL
  
  # if data has rownames, remove
  rownames(data) <- c()
  
  
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
  
  
  forest <- foreach(b = 1:B) %dopar% {          # forest will be list of B trees
    
    source("get_grf.R")
    
    #print(paste(forest_type, "tree number", b, sep=" "))
    
    # draw subsample without replacement
    subsample_idx <- sort(sample.int(n=nrow(data), size=subsample_size, replace=FALSE))
    I <- data[subsample_idx,]
    
    # obtain tree using I
    get_tree(target=target, data=I, min_leaf_size=min_leaf_size, min_node_size=min_node_size, m_try=m_try, func=func, lambda_split=lambda_split)
    
  }
  
  if (get_ci == TRUE){
    
    
    forest_blb <- foreach(g = 1:(B/l)) %dopar% {
      
      source("get_grf.R")
      
      #print(paste("bag nr. ", g, " of ", (B/l), sep=""))
      
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
      bag
    }
  }
  
  if (get_ci == FALSE) forest_blb = NULL
  
  # stop cluster again 
  parallel::stopCluster(cl = my.cluster)
  
  ########################## end of parallel computing #########################
  
  # return forest and further needed information
  grf_output <- list(target, data, forest_type, B, forest, get_ci, forest_blb, l)
  
  return(grf_output)
}