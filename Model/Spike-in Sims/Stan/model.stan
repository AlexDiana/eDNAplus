data {
  // sampling unit-level occupancy covariates
  int<lower = 1> n;
  int<lower = 1> M;
  int<lower = 1> K;
  int<lower = 1> S;
  int<lower = 1> S_star;
  
  matrix[n * M, S_star] v_star;
  // data  
  // matrix[sumM, S] v;
  real y[n * M, K, S + S_star];
  
  // priors
  real<lower = 0> tau;
  real<lower = 0> sigma;
  real<lower = 0> sigma_y;
  real<lower = 0> sigma_u;
  real<lower = 0> sigma_lambda;
  
}

parameters {
  matrix[n, S] logz;
  matrix[n * M, S] v;
  matrix[n * M, K] u;
  vector[S + S_star] lambda;
}


model {
  
  for(i in 1:n){
    for(m in 1:M){
      for(k in 1:K){
        u[(i - 1)*M + m,k] ~ normal(0, sigma_u);
      }
    }
  }
  
  for(j in 1:(S+S_star)){
      lambda[j] ~ normal(0, sigma_lambda);
    }
  
  for(i in 1:n){
    for(j in 1:S){
      logz[i,j] ~ normal(0, tau);
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
            
            y[(i - 1)*M + m, k, j] ~ normal(lambda[j] + v[(i - 1)*M + m, j] + u[(i - 1)*M + m, k], 
            sigma_y);
            
          }
          
          for(j in 1:S_star){      
            
            y[(i - 1)*M + m, k, S + j] ~ normal(lambda[S + j] + v_star[(i - 1)*M + m, j] + u[(i - 1)*M + m, k], 
            sigma_y);
            
          }
      }
    }
  }
  
}
