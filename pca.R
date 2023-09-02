#### LARS algorithm

lars_func <- function(y, X, k_max, center_y, center_and_std_X){
  
  # Centering and standardizing
  if (center_y == TRUE) y <- y - mean(y)
  if (center_and_std_X == TRUE){
    
    X_demeaned <- (X - apply(X, MARGIN=2, FUN=mean))    # demeaning
    sd_X <- apply((X_demeaned^2), MARGIN=2, FUN=sqrt)   # standard deviation
    X <- X_demeaned/sd_X
    
  }
  
  # predefinitions
  m <- ncol(X)
  n <- nrow(X)
  k <- 0                  # step counter
  
  # initializations
  A <- c()                # active set 
  A_c <- c(1:ncol(X))     # inactive set 
  
  mu <- 0                 # initial prediction
  X_k <- c()              # initial active set regressor matrix 
  
  while (k < k_max & k < min(m, (n-1))){
    
    # update k
    k <- k + 1
    
    # obtain current residuals 
    r_current <- y - mu
    
    
    ## identify next active variable 
    
    c <- c()                              # current correlations (for variables in A_c)
    
    for (j in A_c){
      
      c <- c(c, X[,j] %*% r_current)
      
    }

    j_next <- A_c[which.max(abs(c))]      # index of next active x_j
    C <- max(abs(c))                      # greatest absolute current correlation at step k 
    
    # store sign 
    s <- sign(max(c))
    
    # Update sets
    A <- append(A, j_next)
    A_c <- A_c[A_c != j_next]
    
    ## Obtain components (2.4)-(2.6) in Efron et. al. (2004)
    
    X_k <- cbind(X_k, (X[, j_next]*s))
    
    G_k <- t(X_k) %*% X_k
    G_k_inv <- solve(G_k)         # alternatively, chol2inv(G_k)
    
    ones_k <- rep(1, length(A))
    A_k <- as.numeric((t(ones_k) %*% G_k_inv %*% ones_k)^(-1/2))
    
    u_k = A_k * X_k %*% G_k_inv %*% ones_k              # equiangular vector
    a <- t(X) %*% u_k                                   # inner product vector
    
    ## update mu
    
    # obtain gamma
    gamma_candidates <- c()
    for (j in A_c){
      
      c_j <- X[,j] %*% r_current
      a_j <- a[j]
      
      gamma_candidates <- c(gamma_candidates, (C - c_j)/(A_k - a_j), (C + c_j)/(A_k + a_j))
      
    }
    gamma <- min(gamma_candidates[gamma_candidates > 0])
    
    # update estimate mu
    mu <- mu + gamma*u_k
  
  }
  return(A)
}


#### Factor extraction - Stock and Watson (2002)

dfm_pca <- function(X, r){
  # X: (TxN) matrix of variables
  # r: number of factors to extract
  
  T <- nrow(X)  
  N <- ncol(X)  # dimension of multiple time series
  
  # Obtain Lambda_hat
  Sigma_X <- (1/T) * t(X) %*% X
  
  eigenvectors <- eigen(Sigma_X)$vectors[,1:r]

  Lambda_hat <- eigenvectors  * sqrt(N)            # * sqrt(N) is following dfmpc package approach

  # Obtain F_hat 
  F_hat <- X %*% Lambda_hat * (1/N)
  
  colnames <- c()
  for (f in 1:r) colnames <- c(colnames, paste("f", f, sep=""))
  colnames(F_hat) <- colnames
  
  return(list("F"=F_hat, "L"=Lambda_hat))
}
