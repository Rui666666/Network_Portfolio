# Portfolio without short selling

# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("/Users/ruiren/Documents/Github/Network_Portfolio/Review_Response/Reviewer1_comment4_Shortselling")

p=80
select = read.csv("stock_select_80.csv",header = TRUE)[[1]]
stock_select =as.numeric(select)

# Number of observations to generate
n <- 500
n_outsample <- 250

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# load data
prices<-read.csv("SP500 securities_up_20230306.csv")
ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))

#return
return<- Return.calculate(ZOO, method="log")
return<- return[-1, ]
returnstd<-xts(return)


data = returnstd[1:n, stock_select]
data_out = returnstd[(n+1):(n+n_outsample), stock_select]

# Estimate mean vector
estimated_mu <- colMeans(data)

# Estimate correlation matrix
estimated_A <- cor(data)

# Estimate covariance matrix
estimated_Sigma <- cov(data)

# Print results
cat("Estimated mean vector:\n", estimated_mu, "\n")
cat("Estimated correlation matrix:\n")
print(estimated_A)

# Estimate Eigenvector centrality
EC_DS =linfun3_1(estimated_A-diag(1,p,p)-diag(max(eigen(estimated_A)$value),p,p),
                 rep(0,p),
                 lambda=0.45, #Danzig-type hyperparameter
                 abs(eigen(estimated_A-diag(1,p,p)-diag(max(eigen(estimated_A)$value),p,p))$vector[1,1])
)
EC_DS=EC_DS/max(EC_DS)

###### minimum variance portfolio  #####
portf_minVar =globalMin.portfolio(estimated_mu,estimated_Sigma,FALSE)
w =portf_minVar$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_minVar<-rowSums(aus)+1
cumureturn_minVar<-cumprod(return_minVar)
w_minVar<-w

###### mean variance portfolio  ######
portf_meanVar =efficient.portfolio(estimated_mu,estimated_Sigma,0,FALSE)
w =portf_meanVar$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_meanVar<-rowSums(aus)+1
cumureturn_meanVar<-cumprod(return_meanVar)
w_meanVar<-w

###### equally weighted portfolio #####
w =matrix(1/p,1,p)
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_equal<-rowSums(aus)+1
cumureturn_equal<-cumprod(return_equal)
w_equal<-w

###### 2 constraints network portfolio #####
phi_star = quantile(EC_DS,0.35) 
net.gmin.port = network.efficient.portfolio(EC_DS, estimated_Sigma, phi_star,FALSE)
w =net.gmin.port$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_vary_with_phi<-rowSums(aus)+1
cumureturn_network_vary_with_phi<-cumprod(return_network_vary_with_phi)
w_network_vary_with_phi<-w

###### 3 constraints network portfolio #####
mu_star = mean(estimated_mu)

net.gmin.port = network.3constraint.portfolio(EC_DS, estimated_mu, estimated_Sigma, phi_star, mu_star, FALSE)
w =net.gmin.port$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_3constraint_plugin<-rowSums(aus)+1
cumureturn_network_3constraint_plugin<-cumprod(return_network_3constraint_plugin)
w_network_3constraint_plugin<-w

##### Compare results ######
tb = rbind(cbind(cumureturn_minVar[n_outsample]/n_outsample*252*100,
                 cumureturn_meanVar[n_outsample]/n_outsample*252*100,
                 cumureturn_equal[n_outsample]/n_outsample*252*100,
                 cumureturn_network_vary_with_phi[n_outsample]/n_outsample*252*100,
                 cumureturn_network_3constraint_plugin[n_outsample]/n_outsample*252*100
),
cbind(std(return_minVar)*sqrt(252)*100,
      std(return_meanVar)*sqrt(252)*100,
      std(return_equal)*sqrt(252)*100,
      std(return_network_vary_with_phi)*sqrt(252)*100,
      std(return_network_3constraint_plugin)*sqrt(252)*100
),
cbind(cumureturn_minVar[n_outsample]/n_outsample*252/(std(return_minVar)*sqrt(252))*100,
      cumureturn_meanVar[n_outsample]/n_outsample*252/(std(return_meanVar)*sqrt(252))*100,
      cumureturn_equal[n_outsample]/n_outsample*252/(std(return_equal)*sqrt(252))*100,
      cumureturn_network_vary_with_phi[n_outsample]/n_outsample*252/(std(return_network_vary_with_phi)*sqrt(252))*100,
      cumureturn_network_3constraint_plugin[n_outsample]/n_outsample*252/(std(return_network_3constraint_plugin)*sqrt(252))*100),
cbind(maxDrawdown(return_minVar-1),
      maxDrawdown(return_meanVar-1),
      maxDrawdown(return_equal-1),
      maxDrawdown(return_network_vary_with_phi-1),
      maxDrawdown(return_network_3constraint_plugin-1)
)
)
rownames(tb) = c()
xtable(tb,digits = 2)

