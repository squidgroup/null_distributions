data {
//Phenotype dataset
  // Number of clusters (an integer)
   int<lower=0> n; // number of observations for phenotypes
   int<lower=0> n_ind; //number of individual year combinations
   
  // Clusters indentifiers (an integer)
   int<lower=0> individual[n];  //  Ringnr year identites for z
   
  // Continuous outcome
   real x[n];  // phenotypic observations
   real y[n];  // phenotypic observations
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real beta_0; //intercept
     real beta_1; //intercept
  
   // Random effects
   matrix[2,n_ind]         zI; //matrix flok year year blups (intercepts and slopes)
   vector<lower=0>[2]      sigma_I; // sd flok year intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;   // residual var y
 }
 
 transformed parameters{
    matrix[2,n_ind] I; //  Unsclaed blups intercept and slope for w island year
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
}
 
model {
// Create vector of predicted values
   real u_y[n]; // predicted values for individual year phenotype
   
   for (i in 1:n) 
     u_y[i]  = beta_0  +  I[1, individual[i]]  + (beta_1+ I[2, individual[i]])*x[i];

    to_vector([beta_0, beta_1]) ~ normal(0, 1); 

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 0.5);
    L ~ lkj_corr_cholesky(4);
   
 // Likelihood function
     y~normal(u_y, sigma_e);
}   

generated quantities {
real<lower=0> sigma2_e;
real<lower=0> sigma2_I1;
real<lower=0> sigma2_I2;
matrix[2, 2] Omega_I;
real cov_I;

sigma2_e=sigma_e^2;
sigma2_I1=sigma_I[1]^2;
sigma2_I2=sigma_I[2]^2;
Omega_I = L * L';
cov_I = Omega_I[1,2]*sqrt(sigma2_I2*sigma2_I1);
}

