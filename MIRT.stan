data{
  int<lower = 1> N; // Number of participants in a pilot study
  int<lower = 1> C; // The number of categories for each item
  int<lower = 1> M; // Number of subtests
  int<lower = 1> P1; //Number of items in Subtest 1
  int<lower = 1> P2; //Number of items in Subtest 2
  int<lower = 1> P3; //Number of items in Subtest 3
  int<lower = 1> P4; //Number of items in Subtest 4
  int<lower = 1> P;  //Number of total items
  int<lower = 1> K;  //Number of experts
  real<lower = 0> eta; //Variance for mu
  int<lower = 1, upper=5> Y[N,P];  //responses
  vector[M] Zero; // a vector of Zeros
  matrix[M,M] W; // scale matrix 
  vector[P] mu_0; // experts' prior information
}
parameters{
  cov_matrix[M] Sigma; // covariance matrix for multivariate normal distribution
  vector[P] mu; // item-to-domain correlation after Fisher's transformation
  row_vector[M] f[N]; // four factor scores
  real<lower = -0.999, upper = 0.999> rho[P]; // item-to-domain correlation
  ordered[C-1] thre[P]; // Each item has 4 thresholds.
}

model{
  Sigma ~ inv_wishart (4.0, W); // An inverse Wishart distribution is used as the prior distribution for Sigma
  f ~ multi_normal(Zero, Sigma); 
  // prior information
  mu ~ normal(mu_0, eta);
  for (p in 1:P){
    // We specify the prior distribution of the item-to-domain correlation through the following normal approximation
    rho[p] ~ normal((exp(2*mu[p])-1)/(exp(2*mu[p])+1),16*exp(4*mu[p])/(5*K*(exp(2*mu[p])+1)^4));
    // For each item, we specify the four thresholds through the normal distribution N(0,1).
    for (c in 1:(C-1)){
     thre[p,c] ~ normal(0,1);
   }
  }  
  // The likelihood
  for (i in 1:N){
    // We use "ordered_logistic()" function to transform the continuous variable back to the ordinal variable, using the four thresholds.
    for (p in 1:P1){
      Y[i,p] ~ ordered_logistic(f[i,1]*rho[p],thre[p]); 
    }
    for (p in (P1+1):(P1+P2)){
      Y[i,p] ~ ordered_logistic(f[i,2]*rho[p],thre[p]); 
    }
    for (p in (P1+P2+1):(P1+P2+P3)){
      Y[i,p] ~ ordered_logistic(f[i,3]*rho[p],thre[p]); 
    }
    for (p in (P1+P2+P3+1):P){
      Y[i,p] ~ ordered_logistic(f[i,4]*rho[p],thre[p]); 
    }
  }
}

