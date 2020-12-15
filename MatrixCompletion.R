# Sparse Bayesian Methods of Low- Rank Matrix Estimation

#######################
####### Function of Matrix Completion
#######################

# Y is the incomplete matrix which contains NA
# a and b are hyperparameters. they should be small.
# tau: criteria for stop iteration

VBMatrixCompletion <- function(Y, a=10^-5, b=10^-5, tau = 10^-5, tol = 10^2){
  # Compute filter matrix and corrupted rate
  m <- nrow(Y)
  n <- ncol(Y)
  P <- matrix(1, nrow = m, ncol = n)
  P[is.na(Y)] <- 0
  p <- 1 - sum(P)/(m*n)
  Y_0 <- Y
  Y_0[is.na(Y_0)] <- 0
  svd_value <- svd(Y_0)
  index <- which(svd_value$d>sum(svd_value$d)*10^-3)
  # initialize A and B
  A <- svd_value$u[,index] %*% diag(sqrt(svd_value$d[index]))
  B <- svd_value$v[,index] %*% diag(sqrt(svd_value$d[index]))
  
  # initialize parameters
  k <- ncol(A)
  Gamma_ab <- matrix(0, nrow = k, ncol = k)
  diag(Gamma_ab) <- 1
  beta <- 1
  Sigma_a <- array(0, dim = c(k, k, m))
  Sigma_b <- array(0, dim = c(k, k, n))
  for (i in 1:m) {
    diag(Sigma_a[,,i])<-1
  }
  for (i in 1:n) {
    diag(Sigma_b[,,i])<-1
  }
  
  X_old <- A %*% t(B)
  times <- 1
  repeat {
    times <- times +1
    temp1 <- matrix(0, nrow = 1, ncol = k)
    temp2 <- matrix(0, nrow = 1, ncol = k)
    for (i in 1:m) {
      index1 <- which(Y_0[i,]!=0)
      Bi <- B[index1,]
      temp <- matrix(0, nrow = k, ncol = k)
      for (l in 1:length(index1)) {
        temp <- temp + Sigma_b[,,index1[l]]
      }
      BiTBi <- t(Bi) %*% Bi + temp
      Sigma_a[,,i] <- solve(beta * BiTBi + Gamma_ab)
      A[i,] <- beta * Y[i,index1] %*% Bi %*% Sigma_a[,,i]
      temp1 <- temp1 + diag(Sigma_a[,,i])
    }
    
    for (j in 1:n) {
      index2 <- which(Y_0[,j]!=0)
      Aj <- A[index2,]
      temp <- matrix(0, nrow = k, ncol = k)
      for (l in 1:length(index2)) {
        temp <- temp + Sigma_a[,,index2[l]]
      }
      AjTAj <- t(Aj) %*% Aj + temp
      Sigma_b[,,j] <- solve(beta * AjTAj + Gamma_ab)
      B[j,] <- beta * t(Y[index2,j]) %*% Aj %*% Sigma_b[,,j]
      temp2 <- temp2 + diag(Sigma_b[,,j])
    }
    
    diag(Gamma_ab) <- (2 *a + m + n)/(2 *b + temp1 + temp2)
    X_new <- A %*% t(B)
    err <- 0
    for (s in 1:m) {
      observed <- which(P[s,]!=0)
      err <- err + sum(diag(t(B[observed,]) %*% B[observed,] %*% Sigma_a[,,s]))
      for (s1 in 1:length(observed)) {
        err <- err + sum(diag(matrix(A[s,],ncol=1) %*% matrix(A[s,],nrow = 1) %*% Sigma_b[,,observed[s1]]))
        err <- err + sum(diag(Sigma_a[,,s] %*% Sigma_b[,,observed[s1]]))
      }
    }
    beta <- (p * m * n)/(sum((Y_0-X_new * P)^2) + err)
    if(sqrt(sum((X_new-X_old)^2)/sum((X_old)^2))<tau) break
    X_old <- X_new
  }
  estimate_index <- which(diag(Gamma_ab)<tol)
  estimate_rank <- length(estimate_index)
  A <- A[,estimate_index]
  B <- B[,estimate_index]
  X_estimate <- A %*% t(B)
  result <- list(X_estimate = X_estimate, A = A, B = B, estimate_rank = estimate_rank)
  return(result)
}

####################
### Functions of generating simulations
####################

GenerateX <- function(m, n, k, standard_error = 1){
  A <- matrix(rnorm(m*k, mean = 0, sd = standard_error), nrow = m, ncol = k)
  B <- matrix(rnorm(n*k, mean = 0, sd = standard_error), nrow = n, ncol = k)
  X <- A %*% t(B)
  result <- list(X = X, A = A, B= B)
  return(result)
}

FilterMatrix <- function(m, n, p){
  P <- matrix(rbinom(m*n, 1, 1-p), nrow = m, ncol = n)
  P_NA <- P
  P_NA[P_NA==0] <- NA
  result <- list(P = P, P_NA = P_NA)
  return(result)
}

GenerateNoise <- function(m, n, standard_error){
  noise <- matrix(rnorm(m*n, mean = 0, sd = standard_error), nrow = m, ncol = n)
  return(noise)
}


###########################
####### Simulation1: matrix completion without noise
###########################

# Compare with OptSpace
library(ROptSpace)
library(filling)
m <- 100
n <- 100
r <- c(3,5,7,9,10,11,12,13,15)
p <- 0.2



set.seed(10)

l <- length(r)
recovery_error <- matrix(0, nrow = l, ncol = 3)
colnames(recovery_error) <- c("VSBL","OutSpace","SVT")
estimate_rank <- matrix(0, nrow = l, ncol = 2)
colnames(estimate_rank) <- c("VSBL","OutSpace")

iter <- 20
for (i in 1:l) {
  for (j in 1:iter) {
    GenerateX_value <- GenerateX(m, n, r[i], standard_error = 1)
    Filter_value <- FilterMatrix(m, n, p)
    Y <- GenerateX_value$X * Filter_value$P_NA
    Estimate_value <- VBMatrixCompletion(Y)
    recovery_error[i,1] <- recovery_error[i,1] + sqrt(sum((Estimate_value$X_estimate-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))                    
    estimate_rank[i,1] <- estimate_rank[i,1] + Estimate_value$estimate_rank
    
    opt_result <- OptSpace(Y, showprogress = FALSE)
    estimate_rank[i, 2] <- estimate_rank[i,1] + nrow(opt_result$S)
    opt_matrix <- opt_result$X %*% opt_result$S %*% t(opt_result$Y)
    recovery_error[i, 2] <- recovery_error[i, 2] + sqrt(sum((GenerateX_value$X-opt_matrix)^2))/sqrt(sum((GenerateX_value$X)^2))  
    
    svt_result <- fill.SVT(Y)$X
    recovery_error[i, 3] <- recovery_error[i, 3] + sqrt(sum((svt_result-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))
  }
  recovery_error[i,] <- recovery_error[i,]/iter
  estimate_rank[i,] <- estimate_rank[i,]/iter
}
library(ggplot2)

recovery_error_data <- data.frame(method = rep(c("VSBL","OptSpace","SVT"), each = l), rank = rep(r, 3), recovery_error=as.vector(recovery_error))

ggplot(recovery_error_data, aes(x= rank, y=recovery_error, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

estimate_rank_data <- data.frame(method = rep(c("VSBL","OptSpace"), each = l), rank = rep(r, 2), estimate_rank=as.vector(estimate_rank))

ggplot(estimate_rank_data, aes(x= rank, y=estimate_rank, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

#############################
########## Simulation2: matrix completion
#############################

m <- 100
n <- 100
r <- c(3,5,7,8,10,12,14)
p <- 0.2

l <- length(r)
recovery_error <- matrix(0, nrow = l, ncol = 3)
colnames(recovery_error) <- c("VSBL","OutSpace","SVT")
estimate_rank <- matrix(0, nrow = l, ncol = 2)
colnames(estimate_rank) <- c("VSBL","OutSpace")

iter <- 10
for (i in 1:l) {
  for (j in 1:iter) {
    GenerateX_value <- GenerateX(m, n, r[i], standard_error = 1)
    Filter_value <- FilterMatrix(m, n, p)
    Y <- (GenerateX_value$X + GenerateNoise(m,n, standard_error = 0.05))* Filter_value$P_NA
    Estimate_value <- VBMatrixCompletion(Y)
    recovery_error[i,1] <- recovery_error[i,1] + sqrt(sum((Estimate_value$X_estimate-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))                    
    estimate_rank[i,1] <- estimate_rank[i,1] + Estimate_value$estimate_rank
    
    opt_result <- OptSpace(Y, showprogress = FALSE)
    estimate_rank[i, 2] <- estimate_rank[i,1] + nrow(opt_result$S)
    opt_matrix <- opt_result$X %*% opt_result$S %*% t(opt_result$Y)
    recovery_error[i, 2] <- recovery_error[i, 2] + sqrt(sum((GenerateX_value$X-opt_matrix)^2))/sqrt(sum((GenerateX_value$X)^2))  
    
    svt_result <- fill.SVT(Y)$X
    recovery_error[i, 3] <- recovery_error[i, 3] + sqrt(sum((svt_result-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))
  }
  recovery_error[i,] <- recovery_error[i,]/iter
  estimate_rank[i,] <- estimate_rank[i,]/iter
}
library(ggplot2)

recovery_error_data <- data.frame(method = rep(c("VSBL","OptSpace","SVT"), each = l), rank = rep(r, 3), recovery_error=as.vector(recovery_error))

ggplot(recovery_error_data, aes(x= rank, y=recovery_error, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

estimate_rank_data <- data.frame(method = rep(c("VSBL","OptSpace"), each = l), rank = rep(r, 2), estimate_rank=as.vector(estimate_rank))

ggplot(estimate_rank_data, aes(x= rank, y=estimate_rank, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))


#######################
###### Simulation
#######################

dimension <- seq(from = 10, to = 30, by = 5)
r <- 2
p <- 0.2

l <- length(dimension)
recovery_error <- matrix(0, nrow = l, ncol = 3)
colnames(recovery_error) <- c("VSBL","OutSpace","SVT")
estimate_rank <- matrix(0, nrow = l, ncol = 2)
colnames(estimate_rank) <- c("VSBL","OutSpace")

iter <- 1
for (i in 1:l) {
  for (j in 1:iter) {
    GenerateX_value <- GenerateX(dimension[i], dimension[i], r, standard_error = 1)
    Filter_value <- FilterMatrix(dimension[i], dimension[i], p)
    Y <- (GenerateX_value$X + GenerateNoise(dimension[i],dimension[i], standard_error = 0.05))* Filter_value$P_NA
    Estimate_value <- VBMatrixCompletion(Y)
    recovery_error[i,1] <- recovery_error[i,1] + sqrt(sum((Estimate_value$X_estimate-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))                    
    estimate_rank[i,1] <- estimate_rank[i,1] + Estimate_value$estimate_rank
    
    opt_result <- OptSpace(Y, showprogress = FALSE)
    estimate_rank[i, 2] <- estimate_rank[i,1] + nrow(opt_result$S)
    opt_matrix <- opt_result$X %*% opt_result$S %*% t(opt_result$Y)
    recovery_error[i, 2] <- recovery_error[i, 2] + sqrt(sum((GenerateX_value$X-opt_matrix)^2))/sqrt(sum((GenerateX_value$X)^2))  
    
    svt_result <- fill.SVT(Y)$X
    recovery_error[i, 3] <- recovery_error[i, 3] + sqrt(sum((svt_result-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))
  }
  recovery_error[i,] <- recovery_error[i,]/iter
  estimate_rank[i,] <- estimate_rank[i,]/iter
}

recovery_error_data <- data.frame(method = rep(c("VSBL","OptSpace","SVT"), each = l), rank = rep(r, 3), recovery_error=as.vector(recovery_error))

ggplot(recovery_error_data, aes(x= rank, y=recovery_error, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

estimate_rank_data <- data.frame(method = rep(c("VSBL","OptSpace"), each = l), rank = rep(r, 2), estimate_rank=as.vector(estimate_rank))

ggplot(estimate_rank_data, aes(x= rank, y=estimate_rank, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))




