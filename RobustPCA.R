# Sparse Bayesian Methods of Low- Rank Matrix Estimation


##########################
####### Robust PCA
##########################

VBRobustPCA <- function(Y, a=10^-5, b=10^-5, tau = 10^-5, tol = 100){
  # Compute filter matrix and corrupted rate
  m <- nrow(Y)
  n <- ncol(Y)
  svd_value <- svd(Y)
  index <- which(svd_value$d>sum(svd_value$d)*10^-3)
  # initialize A and B
  A <- svd_value$u[,index] %*% diag(sqrt(svd_value$d[index]))
  B <- svd_value$v[,index] %*% diag(sqrt(svd_value$d[index]))
  
  # initialize parameters
  k <- ncol(A)
  Gamma_ab <- matrix(0, nrow = k, ncol = k)
  diag(Gamma_ab) <- 1
  beta <- 1
  Sigma_a <- matrix(0, nrow = k, ncol = k)
  Sigma_b <- matrix(0, nrow = k, ncol = k)
  diag(Sigma_a) <- 1
  diag(Sigma_b) <- 1
  E <- Y - A %*% t(B)
  alpha <- matrix(runif(m*n), nrow = m, ncol = n)
  
  X_old <- A %*% t(B)
  repeat {
    BTB <- t(B) %*% B + n * Sigma_b
    Sigma_a <- solve(beta * BTB + Gamma_ab)
    A <- beta * (Y - E) %*% B %*% Sigma_a
    
    ATA <- t(A) %*% A + m * Sigma_a
    Sigma_b <- solve(beta * ATA + Gamma_ab)
    B <- beta * t(Y- E) %*% A %*% Sigma_b
    Sigma_e <- 1/(beta + alpha)
    E <- beta * Sigma_e * (Y - A %*% t(B))
    
    alpha <- 1/(E^2 + Sigma_e)
    
    temp1 <- m * diag(Sigma_a) + diag(t(A) %*% A)
    temp2 <- n * diag(Sigma_b) + diag(t(B) %*% B)
    diag(Gamma_ab) <- (2 *a + m + n)/(2 *b + temp1 + temp2)
    
    X_new <- A %*% t(B)
    beta <- (m *n)/(sum((Y-X_new-E)^2) + n * sum(diag(t(A)%*%A%*%Sigma_b)) + m * sum(diag(t(B)%*%B%*%Sigma_a)) + m*n*sum(diag(Sigma_a%*%Sigma_b)) + sum(Sigma_e))
    
    if(sqrt(sum((X_new-X_old)^2)/sum((X_old)^2))<tau) break
    X_old <- X_new
  }
  estimate_index <- which(diag(Gamma_ab)<tol)
  estimate_rank <- length(estimate_index)
  A <- A[,estimate_index]
  B <- B[,estimate_index]
  X_estimate <- A %*% t(B)
  E <- E * (alpha<1000)
  result <- list(X_estimate = X_estimate, E=E,A = A, B = B, estimate_rank = estimate_rank)
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

SparseNoise <- function(m, n, p=0.05){
  x <- floor(m*n*p)
  P <- matrix(0, nrow = m, ncol = n)
  P[sample(m*n,x)]<- runif(x, min = -10, max = 10)
  return(P)
}

DenseNoise <- function(m, n, standard_error){
  noise <- matrix(rnorm(m*n, mean = 0, sd = standard_error), nrow = m, ncol = n)
  return(noise)
}


###########################
####### Simulation1: matrix completion
###########################

# Compare with OptSpace

m <- 200
n <- 200
r <- seq(from = 5, to = 25, by = 2)
p <- 0.05

set.seed(10)

l <- length(r)
recovery_error_X <- matrix(0, nrow = l)
recovery_error_E <- matrix(0, nrow = l)
estimate_rank <- matrix(0, nrow = l)

iter <- 10
for (i in 1:l) {
  for (j in 1:iter) {
    GenerateX_value <- GenerateX(m, n, r[i], standard_error = 1)
    Sparse_value <- SparseNoise(m, n)
    # Add white noise or not
    #Y <- Sparse_value + GenerateX_value$X
    Y <- Sparse_value + GenerateX_value$X + GenerateNoise(m,n,0.01)
    Estimate_value <- VBRobustPCA(Y)
    recovery_error_X[i] <- recovery_error_X[i] + sqrt(sum((Estimate_value$X_estimate-GenerateX_value$X)^2))/sqrt(sum((GenerateX_value$X)^2))
    recovery_error_E[i] <- recovery_error_E[i] + sqrt(sum((Estimate_value$E-Sparse_value)^2))/sqrt(sum((GenerateX_value$X)^2))
    estimate_rank[i] <- estimate_rank[i] + Estimate_value$estimate_rank
  }
  recovery_error_X[i] <- recovery_error_X[i]/iter
  recovery_error_E[i] <- recovery_error_E[i]/iter
  estimate_rank[i] <- estimate_rank[i]/iter
}


library(ggplot2)

recovery_error_X <- data.frame(method = rep("VSBL", each = l), rank = r, recovery_error_X=as.vector(recovery_error_X))

ggplot(recovery_error_X, aes(x= rank, y=recovery_error_X)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

recovery_error_E <- data.frame(method = rep("VSBL", each = l), rank = r, recovery_error_E=as.vector(recovery_error_E))

ggplot(recovery_error_E, aes(x= rank, y=recovery_error_E)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

estimate_rank_data <- data.frame(method = rep("VSBL", each = l), rank = r, estimate_rank=as.vector(estimate_rank))

ggplot(estimate_rank_data, aes(x= rank, y=estimate_rank, group = method)) + geom_line(aes(linetype=method, color = method)) + geom_point(aes(color=method))

###########################
##########Generating in other distributions
###########################

GenerateX_unif <- function(m, n, k, width = 5){
  A <- matrix(runif(m*k, from = -width, to = width), nrow = m, ncol = k)
  B <- matrix(rnorm(n*k, mean = -width, to = width), nrow = n, ncol = k)
  X <- A %*% t(B)
  result <- list(X = X, A = A, B= B)
  return(result)
}









