library(mvtnorm)

rain <- read.table("rainfall.dat", sep = "")

#--------------- 1a --------------# Normal model
  
# Setting initial values
mu0 = 1
sigma0 = 1
tao0 = sigma0
nDraws = 500
gibbs_data = matrix(0, nDraws, 2)
n = length(rain[,1])
v0 = 1

sigma = 1
for (i in 1:nDraws){
  
  # Update mu
  tao_n = 1/(n/sigma0+1/sigma)
  w = (n/sigma)/((n/sigma)+(1/tao0))
  mu_n = w*mean(rain[,1]) + (1-w)*mu0
  mu = rnorm(1, mu_n, tao_n)
  gibbs_data[i,1] = mu
  
  # Update sigma^2
  # Making draws from Inv-chi2
  v_n = v0 + n 
  draws = rchisq(n = 1, df=v_n)
  #Convert to Inv-chi2 dist.
  var = (v0*sigma0+sum((rain-mu)^2))/(n+v0)
  sigma=(n)*var/draws
  gibbs_data[i,2] = sigma
  
}

# Analyze these plots
#plot(gibbs_data[,1], gibbs_data[,2])
par(mfrow=c(1,2))
plot(gibbs_data[,1], xlab="Iteration", ylab="Mu sample")
plot(gibbs_data[,2], xlab="Iteration", ylab="Sigma sample")
last_gibbs = gibbs_data[500,]

#--------------- 1b --------------# Mixture normal model

# Done in separate file: NormalMixtureGibbs

#--------------- 1c --------------# Graphical comparison

# Done in separate file: NormalMixtureGibbs

#--------------- 2a --------------# Metropolis Random Walk for Poisson regression

ebay <- read.table("eBayNumberOfBidderData.dat", sep = "", header =TRUE)

# Getting Maximum likelihood values for the model
model = glm(ebay$nBids~., poisson, data=ebay[,-2])
opti_betas = model$coefficients

#--------------- 2b --------------# Bayesian analysis of the Poisson regression

nPara=dim(ebay[-1])[2]
X=as.matrix(ebay[-1])

LogPostLogistic <- function(beta,y,X,mu,sigma){
  nPara <- length(beta);
  datamatrix <- X%*%beta;
  lambda = exp(datamatrix)
  #Ln of likelihood
  logLik <- sum(y*datamatrix - lambda) #factorial not needed, since it's only used as a constant (Per has confirmed this :D)
  print(logLik)
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  #Ln of prior
  logPrior <- dmvnorm(x=beta, mean = rep(0, length(beta)), sigma=100* solve(t(X) %*% X), log=TRUE);
  return(logLik + logPrior)
}
#Starting values before optimization
startVal <- rep(0,nPara)
mu = rep(0,nPara)
sigma = 100* solve(t(X) %*% X)
#This function returns optimal parameters and corresponding value of the function from these opti parameters
#fnscale = -1 means we have max-problem (and not min)
#Hessian = TRUE meaning we return the second derivate matrix:)
OptimResults<-optim(startVal,LogPostLogistic,gr=NULL,as.vector(ebay$nBids),as.matrix(ebay[,-1]),mu,sigma,method=c("BFGS"), 
                    control=list(fnscale=-1),hessian=TRUE)

opti_beta = OptimResults$par
opti_val = OptimResults$value
opti_hes = -solve(OptimResults$hessian)
opti_beta #same values as earlier
opti_hes

#--------------- 2c --------------# Posterior draws using MCMC

hes = opti_hes

LogPostFunc = function(beta,y,X){
  nPara <- length(beta);
  datamatrix <- X %*% beta;
  lambda = exp(datamatrix)
  
  #Ln of likelihood
  logLik <- sum(y*datamatrix - lambda) #factorial not needed (PERRA)
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  #Ln of prior (OBS tanspose of beta here!)
  logPrior <- dmvnorm(x=t(beta), mean = rep(0, nPara), sigma=100* solve(t(X) %*% X), log=TRUE);
  return(logLik + logPrior)
  
}
sample_beta = t(rmvnorm(1, opti_beta, opti_hes))

AP = exp(LogPostFunc(sample_beta, ebay$nBids, X) - LogPostFunc(opti_beta, ebay$nBids, X) )

RMWSampler=function(func, init_beta, hessian, c, ...){
  
  # Initial values for beta (Sample proposal)
  sample_beta = t(rmvnorm(1, init_beta, hessian * c))
  
  #Get acceptance probability
  AP = exp(func(sample_beta, ...) - func(init_beta, ...))
  
  #Get decision value
  u = runif(1, 0, 1)
  if (u<AP){
    print("New sample!")
    return(sample_beta)
  } else {
    print("No new sample")
    return(init_beta)
  }
}

init_beta = rep(8,9)
beta_matrix=matrix(0,9,5000)
beta_matrix[,1]=init_beta
last_beta = init_beta
tuning=1.3

for (i in 2:5000){
  last_beta = RMWSampler(func = LogPostFunc, init_beta = last_beta, hessian = opti_hes, c = tuning, ebay$nBids, X)
  beta_matrix[,i]=last_beta
}

#Plotting all betas
par(mfrow=c(1,2))
plot(beta_matrix[1,], type='l', xlab='Iteration', ylab = 'Beta 1')
plot(beta_matrix[2,], type='l', xlab='Iteration', ylab = 'Beta 2')
plot(beta_matrix[3,], type='l', xlab='Iteration', ylab = 'Beta 3')
plot(beta_matrix[4,], type='l', xlab='Iteration', ylab = 'Beta 4')
plot(beta_matrix[5,], type='l', xlab='Iteration', ylab = 'Beta 5')
plot(beta_matrix[6,], type='l', xlab='Iteration', ylab = 'Beta 6')
plot(beta_matrix[7,], type='l', xlab='Iteration', ylab = 'Beta 7')
plot(beta_matrix[8,], type='l', xlab='Iteration', ylab = 'Beta 8')
plot(beta_matrix[9,], type='l', xlab='Iteration', ylab = 'Beta 9')

#Compare last betas 
#convergence is close to previous analysis
last_beta

#--------------- 2d --------------# Predictive distribution
#Removing burn-in values when calculating values for convergence by only taking the last 1000 betas
mini_x = c(1,1,1,1,0,0,0,1,0.5)
draws = beta_matrix[,4001:5000]

lambdas = exp(mini_x %*% draws) #exponential since it is our model
distribution = dpois(1, lambdas)
hist(lambdas, breaks = 40, xlab='Lambdas')

#dist = rpois(100, lambdas)
#hist(dist)

# Plot distribution
dist = rep(0,4)
for (i in 1:length(lambdas)){
  dist = dist + dpois(c(0,1,2,3), lambdas[i])
}
dist=dist/length(lambdas) #mean of probability 
# Think this is the best
plot(c(0,1,2,3), dist ) #plot output for probability

# Probability of no bidders
dist[1] #= 0.3579




