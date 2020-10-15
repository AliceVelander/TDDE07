setwd("~/Studier+Studentliv/labb Bay")

##-----------------1------------------##
#1a
#Likelihood with bernoulli distribution
likelihood = function(theta, s, n){
  like=theta^s*(1-theta)^(n-s)
  return(like)
}

#prior is Beta(2,2) distributed
prior= function(theta, alpha, beta){
  prio = theta^(alpha-1)*(1-theta)^(beta-1)
  return(prio)
}

theta = seq(0,1,0.01)
posterior = prior(theta,2,2) * likelihood(theta,5,20) 
plot(x=theta, y=posterior, type="l")

alpha=2+5
beta=2+20
expected=alpha/(alpha+beta) #=0.2413793
variance=alpha*beta/((alpha+beta)^2*(alpha+beta+1)) #=0.006103845

vec_mean=c()
vec_sd=c()
for (n in 1:10000){
  func = rbeta(n,2+5,2+20,ncp=0)
  mean=mean(func)
  var = sd(func)^2
  vec_mean=c(vec_mean, mean)
  vec_sd=c(vec_sd, var)
}

plot(x=1:10000, y=vec_mean, type='b', xlab="Number of random draws", ylab="Standard mean")
plot(x=1:10000, y=vec_sd, type='b', xlab="Number of random draws", ylab="Varaince")
 

#1b
#calculate the probability that alpha > 0.3
#We have distribution for beta 
theta=rbeta(1000,2+5,2+20,ncp=0)
larger=theta>0.3
mean(larger) #=0.245

prob = 1 - pbeta(q=0.3,shape1=7,shape2=22)
prob 

#1c 
new_log = log(theta/(1-theta))
hist(new_log, xlab="log-odds values") 
plot(density(new_log), xlab="log-odds values") 

##-----------------2------------------##
#install.packages("invgamma")
library(invgamma)
#install.packages("LaplacesDemon")
library(LaplacesDemon)
#2a
#Known/Given data y
y=c(44,25,45,52,30,63,19,50,34,67)
n= length(y)
mean=3.7
y_mean =mean(y)
nDraws=10000

draws= rchisq(nDraws, df=n, ncp = 0)
#Convert to Inv-chi2 dist.
sample_var=sum((log(y)-mean)^2)/n
sigma_sq=(n)*sample_var/draws

#Plot histograms from draws
hist(sigma_sq, xlim=c(0,1), breaks=32, xlab="sigma squared")

#Compare with theoretical
interval = seq(min(sigma_sq), max(sigma_sq), 0.001)
invschi = dinvchisq(interval, n, sample_var)
hist(invschi, breaks=40, xlab="sigma squared")

#2b
#Gini coefficient 
normal=sqrt(sigma_sq)/sum(sqrt(sigma_sq)) 
gini=2*pnorm(normal/sqrt(2),0,1)-1 
hist(gini)
mean_gini=mean(gini)
plot(density(gini))

#Very equal income. mean = 5.641896e-05

#2c
sort_gini=sort(gini)
#plot(density(sort_gini[501:9500]))
lower=sort_gini[500] #=3.839527e-05
upper=sort_gini[9500] #=8.335225e-05

#View whole distribution with marks of lower and upper tail
hist(sort_gini, breaks=40, xlab="Gini values with marked lower and higher tail")
abline(v=lower, lwd=1, lty=2, col = "red")
abline(v=upper, lwd=1, lty=2, col = "red")

#View 90% of distribution (without tails)
hist(sort_gini[500:9500], breaks=30)

#Do a kernel density estimation to compute a 90% HDP for gini 
sort_density_gini=density(sort_gini)
norm_gini=cumsum(sort_density_gini$y)/sum(sort_density_gini$y)
low_index = min(which(norm_gini>0.05))
low_val=sort_density_gini$x[low_index]
large_index=min(which(norm_gini>0.95))
high_val=sort_density_gini$x[large_index] 

hist(gini, breaks=40, xlab="Gini values with marked HPD interval")
abline(v=low_val, lwd=1, lty=2, col = "blue")
abline(v=high_val, lwd=1, lty=2, col = "blue")
#abline(v=lower, lwd=1, lty=2, col = "blue")
#abline(v=upper, lwd=1, lty=2, col = "blue")

##-----------------3------------------##
#3a
data_radian = c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
mu=2.39

#distribution von mises
von_mises = function(kappa, y, mu){
  tall=exp(kappa*cos(y-mu))
  num=2*pi*besselI(x=kappa, nu=0)
  return(tall/num)
}
#Kappa deside how wuch wind outline 
#try out values that grasp all essentails in graph below :)
kappa_vals = seq(0,7, 0.01) 

#Alpha = 1 given for prior
prior = dexp(kappa_vals, 1)

likelihood = function(kappa) {
  return (prod(von_mises(kappa, data_radian, mu)))
}

#calculate and plot the posterior for all of the k-values (grid)
like_values = sapply(kappa_vals, likelihood)
posterior = like_values * prior
plot(x=kappa_vals, y= posterior, type="l", xlab="values for kappa", ylab= "posterior")

#3b
#Get the mode of the posterior
mode= (which.max(posterior))
kappa_vals[mode] 

