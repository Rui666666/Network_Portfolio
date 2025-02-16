# ===================================================
# Numeric study for optimization (11) and (12)
# Table 1 in Appendix 1: 
# Mean and standard deviation of relative errors 
# ===================================================

rm(list = ls())
library(quadprog)
library(clusterGeneration)
library(igraph)
library(Matrix)


rep = 1000
Obj_1 = matrix(NA,rep,1)
Obj_2 = matrix(NA,rep,1)

for (i in c(1:rep)){
  # Generate a random positive definite matrix
  random_pd_matrix <- genPositiveDefMat(dim = 8)$Sigma
  # Convert the positive definite matrix to a correlation matrix
  random_pd_matrix = abs(random_pd_matrix)
  A <- cov2cor(random_pd_matrix)
  diag(A) = 0
  
  # Define diagonal matrix Gamma 
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  eigen_centrality(g)$vector
  Gamma <- diag(eigen_centrality(g)$vector)  
  
  ## Optimization problem (11) 
  Sigma_xi = random_pd_matrix
  Q = A %*% Gamma %*% Sigma_xi %*% t(Gamma) %*% t(A)+ A %*% Gamma%*% Sigma_xi+Sigma_xi%*%t(Gamma)%*% t(A)
  
  # Define matrices
  Dmat <- Sigma_xi + Q  # Dmat is the sum of Sigma_xi and Q
  dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
  Amat <- matrix(1, 8, 1)  # Constraint matrix for the sum of weights equal to 1
  bvec <- 1  # Right-hand side of the constraint (1 in this case)
  # Solve quadratic programming problem
  #sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  repeat{
    result = tryCatch({
      sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
      sol
      }, error = function(e){
        return(NULL)
        })
    break
  }

  # Extract solution vector
  w_1 <- sol$solution
  Obj_1[i,1]= t(w_1)%*%(Sigma_xi+Q)%*%w_1
  
  ## Optimization Problem (12)
  M =  t(w_1) %*%A%*% Gamma %*%matrix(1, 8, 1) 
  Dmat <- Sigma_xi  
  dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
  Amat <- cbind(matrix(1, 8, 1), -A %*% Gamma %*%matrix(1, 8, 1))  
  bvec <- c(1,-M)  # Right-hand side of the constraint (1 and M in this case), lets set M = 100 for example
  # Solve quadratic programming problem
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  w_2 <- sol$solution
  problem2_min_value = sol$value
  Obj_2[i,1] = t(w_2)%*%(Sigma_xi+Q)%*%w_2
}

diff = abs(Obj_1-Obj_2)
rl_err = abs(Obj_1-Obj_2)/Obj_1

# mean of relative errors
mean(rl_err)
# standard deviation of relative errors
sd(rl_err)




