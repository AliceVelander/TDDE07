library(mvtnorm)

#--------------- 1a --------------# AR(1)-process

# Writing a function that simulates from the AR(1)-process
ar_process = function(mu, fi, sigma, T){
  x = rep(0,T)
  x[1] = mu
  for (i in 2:T){
    x[i] = mu + fi* (x[i-1]-mu) + rnorm(1,0, sqrt(sigma))
  }
  return(x)
}


data = ar_process(20, 0.1, 2, 200)
plot(data, main="Data with fi=0.1")
par(mfrow = c(1,2))
data = ar_process(20, 1, 2, 200)
plot(data, main="1: Data with fi=1")
data = ar_process(20, 1, 2, 200)
plot(data, main="2: Data with fi=1")

# Effect of fi: Small fi reduces the effect of the previous value. Meaning that the plots seems more random than with a
# higher fi. A large fi increases the effect of the previous value. The graph then looks like a coherent line

#--------------- 1b --------------# AR(1)-process

data_1 = ar_process(10, 0.3, 2, 200)
data_2 = ar_process(10, 0.95, 2, 200)

# Writing the Stan model
StanModel = '
data {
int<lower=0> T; // Number of observations
real y[T];
}
parameters {
real mu;
real<lower=0> sigma2;
real<lower=-1, upper=1> fi;
}
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2

y[1] ~ normal(mu,sqrt(sigma2));
for(i in 2:T)
y[i] ~ normal(mu + fi* (y[i-1] - mu), sqrt(sigma2));
}'

## MCMC for fi = 0.3
data = list(T=200, y=data_1)
niter = 2000
burnin = 1000
MCMC1 = stan(model_code=StanModel, data=data, warmup=burnin, iter=niter, chains=4)

# Print the fitted model
print(MCMC1,digits_summary=3)
# Extract posterior samples
postDraws <- extract(MCMC1)
# Do traceplots of the first chain
par(mfrow = c(1,3))
plot(postDraws$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot of mu")
plot(postDraws$fi[1:(niter-burnin)],type="l",ylab="fi",main="Traceplot of fi")
plot(postDraws$sigma2[1:(niter-burnin)],type="l",ylab="sigma",main="Traceplot of sigma")

# Do automatic traceplots of all chains
traceplot(MCMC1)

# 95% credible interval
sort_mu = sort(postDraws$mu[1:(niter-burnin)])
interval_mu1 = c(sort_mu[(length(sort_mu)*0.025)],sort_mu[(length(sort_mu)*0.975)])
mean_mu1 = mean(sort_mu)

sort_fi1 = sort(postDraws$fi[1:(niter-burnin)])
interval_fi1 = c(sort_fi1[(length(sort_fi1)*0.025)],sort_fi1[(length(sort_fi1)*0.975)])
mean_fi1= mean(sort_fi1)

sort_sigma = sort(postDraws$sigma2[1:(niter-burnin)])
interval_sigma1 = c(sort_sigma[(length(sort_sigma)*0.025)],sort_sigma[(length(sort_sigma)*0.975)])
mean_sigma1 = mean(sort_sigma)

  ## MCMC for fi = 0.95
data = list(T=200, y=data_2)
MCMC2 = stan(model_code=StanModel, data=data, warmup=burnin, iter=niter, chains=4)

# Print the fitted model
print(MCMC2,digits_summary=3)
# Extract posterior samples
postDraws <- extract(MCMC2)
# Do traceplots of the first chain
par(mfrow = c(1,3))
plot(postDraws$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot of mu")
plot(postDraws$fi[1:(niter-burnin)],type="l",ylab="fi",main="Traceplot of fi")
plot(postDraws$sigma2[1:(niter-burnin)],type="l",ylab="sigma",main="Traceplot of sigma")

# Do automatic traceplots of all chains
traceplot(MCMC2)

# 95% credible interval
sort_mu = sort(postDraws$mu[1:(niter-burnin)])
interval_mu2 = c(sort_mu[(length(sort_mu)*0.025)],sort_mu[(length(sort_mu)*0.975)])
mean_mu2 = mean(sort_mu)

sort_fi = sort(postDraws$fi[1:(niter-burnin)])
interval_fi2 = c(sort_fi[(length(sort_fi)*0.025)],sort_fi[(length(sort_fi)*0.975)])
mean_fi2= mean(sort_fi)

sort_sigma = sort(postDraws$sigma2[1:(niter-burnin)])
interval_sigma2 = c(sort_sigma[(length(sort_sigma)*0.025)],sort_sigma[(length(sort_sigma)*0.975)])
mean_sigma2 = mean(sort_sigma)

#--------------- 1c --------------# Campylobacter
camp <- read.csv("campy.dat", sep = "")

# Writing the Stan model
CampModel = '
data {
int<lower=0> T; // Number of observations
int y[T];
}
parameters {
real mu;
real<lower=0> sigma2;
real<lower=-1, upper=1> fi;
real norm;
vector[T] x;
}
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(140,2); // Scaled-inv-chi2 with nu 1, sigma 2

x[1] ~ normal(mu,sqrt(sigma2));
y[1] ~ poisson(exp(x[1]));

for(i in 2:T){
x[i] ~ normal(mu + fi* (y[i-1] - mu), sqrt(sigma2));
y[i] ~ poisson(exp(x[i]));
}

}'

data = list(T=140, y=camp[,1])
niter = 2000
burnin = 1000
MCMC_camp = stan(model_code=CampModel, data=data, warmup=burnin, iter=niter, chains=4)

summary = summary(MCMC_camp)
mean = summary$summary[5:144,1]
lower = summary$summary[5:144,4]
upper = summary$summary[5:144,8]

# Print the fitted model
print(MCMC_camp,digits_summary=3)
# Extract posterior samples
postDraws <- extract(MCMC_camp)
plot(postDraws$mu, type='l', main ='Traceplot for my', xlab='mu')
plot(postDraws$fi, type='l', main ='Traceplot for fi', xlab='fi')
plot(postDraws$sigma2, type='l', main ='Traceplot for sigma', xlab='sigma')

plot(camp[,1], type='l')
lines(x = seq(1,140), exp(mean), lwd=1, lty=2, col = "yellow")
lines(x = seq(1,140), exp(lower), lwd=1, lty=2, col = "red")
lines(x = seq(1,140), exp(upper), lwd=1, lty=2, col = "green")

#--------------- 1d --------------# New prior

# New Stan Model
CampModel_prior = '
data {
int<lower=0> T; // Number of observations
int y[T];
}
parameters {
real mu;
real<lower=0> sigma2;
real<lower=-1, upper=1> fi;
real norm;
vector[T] x;
}
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(140,0.1); // Adjust sigma for increments to be small

x[1] ~ normal(mu,sqrt(sigma2));
y[1] ~ poisson(exp(x[1]));

for(i in 2:T){
x[i] ~ normal(mu + fi* (y[i-1] - mu), sqrt(sigma2));
y[i] ~ poisson(exp(x[i]));
}

}'

data = list(T=140, y=camp[,1])
niter = 2000
burnin = 1000
MCMC_camp_prior = stan(model_code=CampModel_prior, data=data, warmup=burnin, iter=niter, chains=4)

summary_prior = summary(MCMC_camp_prior)
mean_prior = summary_prior$summary[5:144,1]
lower_prior = summary_prior$summary[5:144,4]
upper_prior = summary_prior$summary[5:144,8]

plot(camp[,1], type='l', main='Series plot with new prior', ylab='Value')
lines(x = seq(1,140), exp(mean_prior), lwd=1, lty=2, col = "yellow")
lines(x = seq(1,140), exp(lower_prior), lwd=1, lty=2, col = "red")
lines(x = seq(1,140), exp(upper_prior), lwd=1, lty=2, col = "green")


