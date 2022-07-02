setwd("C:/Users/jihan/downloads")

library(rstan)
library(MASS)
library(base)
library(shinystan)
library(MCMCpack)
library(stats)
library(lavaan)
library(psych)

N = 200                    # Number of participants in a pilot study
C = 5                     # The number of categories for each item
M = 4                     # Number of subtests
P1 = P2 = P3 = P4 = 5     # Number of items in each subtest
P = P1 + P2 + P3 + P4     # Number of total items
K = 3                     # Number of experts
eta = sqrt(1/(5*K))                         # Variance of mu
T.rho_1 = T.rho_2 = T.rho_3 = T.rho_4 = c(0.5,0.5,0.5,0.3,0.7)  # True item-to-domain correlations
T.rho = c(T.rho_1,T.rho_2,T.rho_3,T.rho_4)
rho0_1 = rho0_2 = rho0_3 = rho0_4 = c(0.5,0.5,0.5,0.3,0.7)
rho0 = c(rho0_1,rho0_2,rho0_3,rho0_4)
item = paste0("Item",1:P)
`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}
equation1 = item[1]
equation2 = item[P1+1]
equation3 = item[P1+P2+1]
equation4 = item[P1+P2+P3+1]
for (Item in item[2:P1]) {equation1 = equation1 + '+' +Item}
for (Item in item[(P1+2):(P1+P2)]) {equation2 = equation2 + '+' +Item}
for (Item in item[(P1+P2+2):(P1+P2+P3)]) {equation3 = equation3 + '+' +Item}
for (Item in item[(P1+P2+P3+2):P]) {equation4 = equation4 + '+' +Item}


#number of simulations
S=20
    
# Priors for MCMC
Zero <- rep(0,M) # a vector of Zeros (fixed means for person parameter)
W <- matrix(c(1,0,0,0,
              0,1,0,0,
              0,0,1,0,
              0,0,0,1),4,4) # binary unit matrix
mu_0 <- c(0.5*log((1+T.rho_1)/(1-T.rho_1)),0.5*log((1+T.rho_2)/(1-T.rho_2)),0.5*log((1+T.rho_3)/(1-T.rho_3)),0.5*log((1+T.rho_4)/(1-T.rho_4))) 
# mu_0 is experts' prior information. It is the true item-to-domain correlation after Fisher's transformation.
# mu_0 <- mu_0-.1

#Define all needed matrices
sim.ystar=sim.data=matrix(NA,N,P)
rhoSIM.mean=cfaSIM.mean=matrix(nrow=S,ncol=P)
rhoSIM.se=cfaSIM.se=matrix(nrow=S,ncol=P)
MOBID.corr=cfa.corr=matrix(NA,S,1)
CFAfail=CFAbadest=matrix(NA,S,1)

for (s in 1:S){
  set.seed(s+1000)
  mean <- rep(0,M)                                     # Specify the means of the variables
  sigma <- matrix(c(1, .8, .8, .8,
                    .8, 1, .8, .8,
                    .8, .8, 1, .8,
                    .8, .8, .8, 1), ncol = M)          # Specify the covariance matrix of the variables
  ft <- mvrnorm(n = N, mu = mean, Sigma = sigma)
  for (j in 1:P1){
    sim.ystar[,j]=T.rho_1[j]*ft[,1]+rnorm(N,0,sqrt(1-T.rho_1[j]^2))
  }
  for (j in (P1+1):(P1+P2)){
    sim.ystar[,j]=T.rho_2[j-P1]*ft[,2]+rnorm(N,0,sqrt(1-T.rho_2[j-P1]^2))
  }
  for (j in (P1+P2+1):(P1+P2+P3)){
    sim.ystar[,j]=T.rho_3[j-P1-P2]*ft[,3]+rnorm(N,0,sqrt(1-T.rho_3[j-P1-P2]^2))
  }
  for (j in (P1+P2+P3+1):P){
    sim.ystar[,j]=T.rho_4[j-P1-P2-P3]*ft[,4]+rnorm(N,0,sqrt(1-T.rho_4[j-P1-P2-P3]^2))
  }
  
  #Transform continuous y* back to ordinal y
  if(C == 2){
    for (j in 1:P) {
      sim.data[,j]=cut(sim.ystar[,j],breaks=qnorm(seq(0,1,0.5)),
                       include.lowest=TRUE,labels=FALSE,ordered_result=TRUE)
    }
    
  } else {
    for (j in 1:P) {
      sim.data[,j]=cut(sim.ystar[,j],breaks=qnorm(seq(0,C)/C),
                       include.lowest=TRUE,labels=FALSE,ordered_result=TRUE)
    }
  }
  
  # MCMC for parameter estimation
  Y = sim.data
  stan_MIRT_data <- list(N=N, C=C, M=M, P1=P1, P2=P2, P3=P3, P4=P4, P=P, K=K, eta=eta,
                         Y=Y, Zero=Zero, W=W, mu_0=mu_0)
  fit_MIRT <- stan(file="MIRT.stan",data=stan_MIRT_data,cores = 4, chains = 2,warmup=250,iter=500,refresh=100)
  fit_ss <- extract(fit_MIRT, permuted = TRUE)
  rhosim <- rep(NA,P)
  for (i in 1:20) {rhosim[i] = mean(fit_ss$rho[,i])}
  rhoSIM.mean[s,] <- rhosim
  
  # Traditional CFA model
  Y_cfa = data.frame(sim.data)
  colnames(Y_cfa) = item
  Y_cfa[,item]=lapply(Y_cfa[,item], ordered)
  model=c(paste("factor1 =~", equation1), 
          paste("factor2 =~", equation2),
          paste("factor3 =~", equation3), 
          paste("factor4 =~", equation4))
  fitCat = cfa(model,data=Y_cfa,estimator="WLSMV",std.lv=TRUE)
  cfaSIM.mean[s,] = parameterEstimates(fitCat)$est[1:P]
  cfaSIM.se[s,] = parameterEstimates(fitCat)$se[1:P]
  
  if(lavInspect(fitCat, "converged")) {CFAfail=0} else {CFAfail=1}
  cfaSIM.mean[s,]=parameterEstimates(fitCat)$est[1:P]
  cfaSIM.se[s,]=parameterEstimates(fitCat)$se[1:P]
  cfafi=predict(fitCat)
  cfa.corr[s]=cor(cfafi,ft)
  
  if(any(CFAfail[s]!=0,na.rm=TRUE)) {
    cfaSIM.mean[s,]=NA
    cfaSIM.se[s,]=NA
  }
  if(any((cfaSIM.mean[s,]>1 | cfaSIM.mean[s,]< (-1)),na.rm=TRUE)){
    CFAbadest[s]=1  
  } else {CFAbadest[s]=0}
  
}

