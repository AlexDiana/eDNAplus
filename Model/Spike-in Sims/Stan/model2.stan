data {
  // sampling unit-level occupancy covariates
  int<lower = 1> n;
  int<lower = 1> M;
  int<lower = 1> K;
  int<lower = 1> S;
  matrix[n, 1] X_z;
  
  // data  
  // matrix[sumM, S] v;
  int<lower=0> y[n * M, K, S];
  
  // priors
  real<lower = 0> tau[S];
  real<lower = 0> sigma[S];
  
}

parameters {
  matrix[ncov_z + 1, S] beta_z;
  matrix[n, S] logz;
  matrix[n * M, S] v;
}

transformed parameters {
  // vector[N] logit_theta = beta_theta * delta_p;
  matrix[n, S] Xzbeta = X_z * beta_z;
  matrix[sumM, S] Xwbeta = X_w * beta_w;
  matrix[sumM, S] Xwbetatheta = X_wtheta * beta_theta;
  matrix[sumM, S] logit_theta;
  for(i in 1:n){
    for(j in 1:S){
      for(m in 1:M_site[i]){
        logit_theta[summ[i] + m, j] = Xwbetatheta[summ[i] + m, j] + 
          beta_theta_1[j] * logz[i,j];  
      }
    }
  }
}

model {
  matrix[sumM, S] log_theta = log_inv_logit(logit_theta);
  matrix[sumM, S] log1m_theta = log1m_inv_logit(logit_theta);
  
  real target_detected;
  real target_notdetected;
  
  
  for(i in 1:n){
    for(j in 1:S){
      logz[i,j] ~ normal(Xzbeta[i,j], tau[j]);
    }
  }
  
  for(j in 1:S){
    for(ncov in 1:ncov_z){
      beta_z[ncov,j] ~ normal(0, 1);  
    }
    for(ncov in 1:ncov_w){
      beta_w[ncov,j] ~ normal(0, 1);  
    }
  }
  
  
  for (i in 1:n) {
    for(j in 1:S){
      for(m in 1:M_site[i]){
        
        target_detected = 0;
        target_notdetected = 0;
        
        
        for(k in 1:K[summ[m] + m]){
          
          target_detected += 
            neg_binomial_2_lpmf(
              y[summ[m] + m, k, j] | exp(lambda[j] + v[summ[m] + m, j]), 
              r_nb[j]);
          
          target_notdetected += 
            neg_binomial_2_lpmf(
              y[summ[m] + m, k, j] | n_0, p_0);
          
        }
        
        target += log_sum_exp(
          log_theta[summ[m] + m,j] + 
            normal_lpdf(v[summ[m] + m] | logz[i, j] + Xwbeta[summ[m] + m,j], sigma[j]) + 
            target_detected,
          log1m_theta[summ[m] + m,j] + 
            normal_lpdf(v[summ[m] + m] | -10, sigma_gamma) + 
            target_notdetected);
        
        // target += log_sum_exp(
          //   log_theta[summ[m] + m,j] + 
            //     normal_lpdf(v[summ[m] + m] | logz[i, j] + Xwbeta[summ[m] + m,j], sigma[j]),
          //   log1m_theta[summ[m] + m,j] + 
            //     normal_lpdf(v[summ[m] + m] | -10, sigma_gamma));
        
      }
    }
  }
  
}
