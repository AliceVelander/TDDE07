library(mvtnorm)

data <- read.csv("Temp_Link.txt", sep = "")
time= data[,1]
temp= data[,2]
n <- 366

##-----------------1------------------##
#1a
# Getting data from the file
# mu is the mu_0 for the beta values
# mu1 affects the height, mu2 shifts the temperature graph, mu3 affects the curvature
mu = c(0,100,-100)
# Larger omega creates small volatility
omega = 0.1
I =diag(3)
omega= I *omega
# nu sets the amount of influence of the prior
v= 6
sigma_2= 1

#100 sample från inv chi 2 - sätter in i ovan uttryck - multiplicerar över en normalfördelning.
#100 sigma 
nDraws=100

# Making draws from Inv-chi2 
draws= rchisq(n = nDraws, df=v) 
#Convert to Inv-chi2 dist. 
sample_var = sigma_2 
sigma_sq=(v)*sample_var/draws

# Making draw som normal dist. beta given sigma2 
#loop for each sigma_sq
#See if the temperatures looks reasonable 
#don't compare with data given 
dev.off()
for (sigma in sigma_sq) {
  beta = rmvnorm(n = 1, mean = mu, sigma*solve(omega))
  time_series = beta[,1]+beta[,2]*time+beta[,3]*time^2 
  plot(x=time, y=time_series,xlim = c(0,1),ylim= c(-40,60), xlab ="", ylab = "", type = 'l') 
  par(new = TRUE)
}


#1b
# Calculating "weights", mu_n and nu_n X = cbind(1, time, time^2)
Y = temp
beta_hat = solve(t(X)%*%X)%*%t(X)%*%Y
mu_n = solve(t(X) %*% X+omega) %*% (t(X) %*% X %*% beta_hat+omega%*%mu)
omega_n = t(X) %*% X+omega
v_n = v + length(time)
sigma_2_n = (v*sigma_2 + (t(Y) %*% Y + t(mu) %*% omega %*% mu - t(mu_n) %*% omega_n %*% mu_n))/v_n

# Making draws from Inv-chi2
draws_post= rchisq(n = nDraws, df=v_n)
#Convert to Inv-chi2 dist. 
sigma_sq_post=(v_n)*sigma_2_n/draws_post

#Calculating marginal distributions for beta and sigma2 and showing them in histograms
beta_post_1 = c()
beta_post_2 = c()
beta_post_3 = c()

for (sigma in sigma_sq_post) {
  beta_post = rmvnorm(n = 1, mean = mu_n, sigma*solve(omega_n)) 
  beta_post_1 = c(beta_post_1, beta_post[,1])
  beta_post_2 = c(beta_post_2, beta_post[,2])
  beta_post_3 = c(beta_post_3, beta_post[,3])
}

par(mfrow=c(2,2)) 
hist(beta_post_1) 
hist(beta_post_2) 
hist(beta_post_3) 
hist(sigma_sq_post)

# Scatter plot of the temperature 
time_series_median = c()
tails = matrix(0,2, length(time))

for (t in time) {
  time_series_post =c()
  for (i in 1:length(beta_post_1)) {
    series_day = beta_post_1[i] + beta_post_2[i]*t + beta_post_3[i]*t^2
    time_series_post = c(time_series_post, series_day) 
    }
  tails[,(match((t), time))] = quantile(time_series_post, c(0.025, 0.975))
  time_series_median = c(time_series_median, median(time_series_post))
}

plot(x= time, y = time_series_median, ylim = c(-15, 30), type = 'l')

par(new = TRUE)
plot(x = time, y = temp, xlab = "", ylim = c(-15, 30), ylab = "")
par(new = TRUE)
plot(x = time, y = tails[1,], "", type = "l", col = "blue") 
par(new = TRUE)
plot(x = time, y = tails[2,], "", type = "l", col = "red")

#1c
#Max value when derivative is zero. --> temp ' = 0 = -B1 /(2* B2)
x_hat = -1 * beta_post_2 / (2* beta_post_3) 
hist(x_hat, breaks=20)
x_hat

#1d
#Set ridge regression on prior, meaning mu0 = 0 and omega0 = large (lambda) in order to get more shrinkage and less variance.
# This is equivalent to a penalized likelihood and less overfit.

##-----------------2------------------##
women <- read.table("WomenWork.dat", header=TRUE) 
y = as.vector(women[,1])
X = as.matrix(women[,-1])
nPara=dim(X)[2]
LogPostLogistic <- function(beta,y,X,mu,sigma){ 
  nPara <- length(beta);
  datamatrix <- X%*%beta;
  
  #Ln of likelihood
  logLik <- sum( datamatrix*y -log(1 + exp(datamatrix)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  #Ln of prior
  logPrior <- dmvnorm(beta, matrix(0,nPara,1), sigma, log=TRUE); 
  return(logLik + logPrior)
}

#Starting values before optimization 
startVal <- rep(0,nPara)
mu = rep(0,nPara)
sigma = 10^2* diag(nPara)

#This function returns optimal parameters and corresponding value of the function from these opti parameters
#fnscale = -1 means we have max-problem (and not min)
#Hessian = TRUE meaning we return the second derivate matrix 
OptimResults<-optim(startVal,LogPostLogistic,gr=NULL,y,X,mu,sigma,met hod=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)
opti_beta = OptimResults$par
opti_val = OptimResults$value
opti_hes = -solve(OptimResults$hessian)

# Calculating credible inverval for NSmallChild
nDraws = 1000
NSmallChild_beta_draws = rmvnorm(n = nDraws, mean = opti_beta, sigma = opti_hes)[,7]
NSmallChilld_quantile = quantile(NSmallChild_beta_draws, c(0.025, 0.975)) 
glmModel <- glm(Work ~ 0 + ., data = women, family = binomial)

#2b
# Making the predictive function
predict_work = function(feature_values, Ndraws, beta, sigma){
  draws = rmvnorm(n= Ndraws, mean = beta, sigma = sigma)
  work = exp(feature_values %*% t(draws))/(1+exp(feature_values %*% t(draws)))
  return(work) 
  }
# Making predictive distribution for the sample women
women = c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
pred = predict_work(women, 1000, opti_beta, opti_hes)
hist(pred, xlab = "Probability to work", main = "Histogram of the probability distribution")

#2c
# Calculating distribution for posterior of binomial = 0
for (i in 1:10){
  pred = ifelse(predict_work(women, 1000, opti_beta, opti_hes)>0.5, 1, 0)
  binomial = binomial + pred 
  }
hist(binomial, xlim=c(0,10), main="Histogram for the number out of 10 women working")

