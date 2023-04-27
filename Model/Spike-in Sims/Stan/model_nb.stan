data {
  // sampling unit-level occupancy covariates
  int<lower = 1> n;
  int<lower = 1> M;
  int<lower = 1> K;
  int<lower = 1> S;
  int<lower = 0> S_star;
  
  matrix[n * M, S_star] v_star;
  // data  
  // matrix[sumM, S] v;
  int<lower = 0> y[n * M, K, S + S_star];
  
  // priors
  real<lower = 0> sigma_u;
  real mean_lambda;
  real<lower = 0> sigma_lambda;
  real<lower = 0> sigma[S];
  real<lower = 0> r[S + S_star];
  
}

parameters {
  matrix[n, S] logz;
  matrix[n * M, S] v;
  matrix[n * M, K] u;
  vector[S + S_star] lambda;
  // real<lower = 0> sigma[S];
  // real<lower = 0> r[S + S_star];
}


model {
  
  for(i in 1:n){
    for(m in 1:M){
      for(k in 1:K){
        u[(i - 1)*M + m,k] ~ normal(0, sigma_u);
      }
    }
  }
  
  for(j in 1:S){
    lambda[j] ~ normal(mean_lambda, sigma_lambda);
  }

  if(S_star > 0){
    for(j in 1:S_star){
      lambda[S + j] ~ normal(mean_lambda, sigma_lambda);
    }
  }
  
  for(i in 1:n){
    for(m in 1:M){
      for(j in 1:S){
        v[(i - 1)*M + m,j] ~ normal(logz[i,j], sigma);
      } 
    }
  }
  
  for (i in 1:n) {
    
    for(m in 1:M){
      
      for(k in 1:K){
        
        for(j in 1:S){      
          
          y[(i - 1)*M + m, k, j] ~ neg_binomial_2(
            exp(lambda[j] + v[(i - 1)*M + m, j] + u[(i - 1)*M + m, k]), r[j]);
          // y[(i - 1)*M + m, k, j] ~ normal(lambda[j] + v[(i - 1)*M + m, j] + u[(i - 1)*M + m, k], 
                                             //                                 sigma_y);
          
        }
        
        if(S_star > 0){
          
          for(j in 1:S_star){      
            
            y[(i - 1)*M + m, k, S + j] ~ neg_binomial_2(
              exp(lambda[S + j] + v_star[(i - 1)*M + m, j] + u[(i - 1)*M + m, k]), r[S + j]);
            
          }
          
        }
        
      }
    }
  }
  
  // for(i in 1:n){
    //   for(j in 1:S){
      //     logz[i,j] ~ normal(0, tau);
      //   }
    // }
  
}
