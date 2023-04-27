// #include <RcppArmadilloExtensions/sample.h>
  
  #include<RcppArmadillo.h>
  
  using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int sample_cpp(arma::vec x, arma::vec probs){
  
  int K = probs.size();
  
  double u = R::runif(0, 1);
  
  double sum = 0;
  int index = -1;
  for(int l = 0; l < K; l++){
    
    index += 1;
    sum += probs[l];
    
    if(u < sum){
      return(x[index]);  
    } 
    
  }
  
  return 0;
}


double logprob_c0(double y, double pi0, double r, double pi){
  
  if(y == 0){
    return(log(pi0));    
  } else {
    return(log(1 - pi0) + R::dnbinom(y - 1, r, 1 - pi, 1));
  }
  
  return(0);
}



// [[Rcpp::export]]
double compute_logprob_y_delta0_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                    int currentK, double n0, double p0, double pi0,
                                    double lambda, double lambdatilde){
  
  double sum = 0;
  for(int k = 0; k < currentK; k++){
    if(c_imk_current[k] == 2){
      sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
    } else {
      // sum += R::dnbinom(y_counts[k], n0, p0, 1);
      sum += logprob_c0(y_counts[k], pi0, n0, p0);
    }
  }
  
  return(sum);
}

// // [[Rcpp::export]]
// double compute_logprob_y_delta1_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                       //                                     double lambda, 
                                       //                                     int currentK, double lambda0, double lambdatilde,
                                       //                                     arma::vec u_im, double v_im){
  //   
    //   double sum = 0;
    //   for(int k = 0; k < currentK; k++){
      //     if(c_imk_current[k] == 1){
        //       sum += R::dpois(y_counts[k], exp(lambda + u_im[k] + v_im), 1);
        //     } else if(c_imk_current[k] == 2){
          //       sum += R::dpois(y_counts[k], lambdatilde, 1);
          //     } else {
            //       sum += R::dpois(y_counts[k], lambda0, 1);
            //     }
      //   }
    //   
      //   return(sum);
    // }

// [[Rcpp::export]]
double compute_logprob_y_delta1_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                    int currentK, double n0, double p0, 
                                    double v_im, double lambdatilde,
                                    double lambda, arma::vec u_im){
  
  double sum = 0;
  for(int k = 0; k < currentK; k++){
    if(c_imk_current[k] == 1){
      sum += R::dpois(y_counts[k], exp(lambda + v_im + u_im[k]), 1);
    } else if(c_imk_current[k] == 2){
      sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
    } else {
      sum += R::dnbinom(y_counts[k], n0, p0, 1);
    }
  }
  
  return(sum);
}

// [[Rcpp::export]]
arma::vec DecToBin_cpp(int l, int currentK){
  
  if(currentK == 0){
    arma::vec toReturn = arma::zeros(1);
    return(toReturn);
  }
  
  arma::vec toReturn = arma::zeros(currentK);
  if(currentK == 4){
    
    if(l == 1){
      toReturn[0] = 1;
    } else if(l == 2){
      toReturn[1] = 1;
    } else if(l == 3){
      toReturn[0] = 1;
      toReturn[1] = 1;
    } else if(l == 4){
      toReturn[2] = 1;
    } else if(l == 5){
      toReturn[0] = 1;
      toReturn[2] = 1;
    } else if(l == 6){
      toReturn[1] = 1;
      toReturn[2] = 1;
    } else if(l == 7){
      toReturn[0] = 1;
      toReturn[1] = 1;
      toReturn[2] = 1;
    } else if(l == 8){
      toReturn[3] = 1;
    } else if(l == 9){
      toReturn[0] = 1;
      toReturn[3] = 1;
    } else if(l == 10){
      toReturn[1] = 1;
      toReturn[3] = 1;
    } else if(l == 11){
      toReturn[0] = 1;
      toReturn[1] = 1;
      toReturn[3] = 1;
    } else if(l == 12){
      toReturn[2] = 1;
      toReturn[3] = 1;
    } else if(l == 13){
      toReturn[0] = 1;
      toReturn[2] = 1;
      toReturn[3] = 1;
    } else if(l == 14){
      toReturn[1] = 1;
      toReturn[2] = 1;
      toReturn[3] = 1;
    } else if(l == 15){
      toReturn[0] = 1;
      toReturn[1] = 1;
      toReturn[2] = 1;
      toReturn[3] = 1;
    } 
    
    
  } else if(currentK == 2){
    
    if(l == 1){
      toReturn[0] = 1;
    } else if(l == 2){
      toReturn[1] = 1;
    } else if(l == 3){
      toReturn[0] = 1;
      toReturn[1] = 1;
    } 
    
  } else if(currentK == 1){
    
    if(l == 1){
      toReturn[0] = 1;
    }
    
  } else if(currentK == 3){
    
    if(l == 1){
      toReturn[0] = 1;
    } else if(l == 2){
      toReturn[1] = 1;
    } else if(l == 3){
      toReturn[0] = 1;
      toReturn[1] = 1;
    } else if(l == 4){
      toReturn[2] = 1;
    } else if(l == 5){
      toReturn[0] = 1;
      toReturn[2] = 1;
    } else if(l == 6){
      toReturn[1] = 1;
      toReturn[2] = 1;
    } else if(l == 7){
      toReturn[0] = 1;
      toReturn[1] = 1;
      toReturn[2] = 1;
    }
    
  }
  
  return(toReturn);
}

// [[Rcpp::export]]
double compute_logprob_y_delta1_rnb_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                        int currentK, double n0, double p0, double pi0,
                                        double r_nb,
                                        double v_im, double lambdatilde,
                                        double lambda, arma::vec u_im){
  
  double sum = 0;
  for(int k = 0; k < currentK; k++){
    if(c_imk_current[k] == 1){
      double pi = exp(lambda + v_im + u_im[k]) / (exp(lambda + v_im + u_im[k]) + r_nb);
      sum += R::dnbinom(y_counts[k], r_nb, 1 - pi, 1);
      // sum += logprob_c0(y_counts[k], pi0, r_nb, pi);
    } else if(c_imk_current[k] == 2){
      sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
    } else {
      // sum += R::dnbinom(y_counts[k], n0, p0, 1);
      sum += logprob_c0(y_counts[k], pi0, n0, p0);
    }
  }
  
  return(sum);
}

// [[Rcpp::export]]
List update_delta_c_d_rjmcmc(arma::cube y, arma::mat v, 
                             arma::vec lambda,
                             arma::vec r_nb,
                             arma::vec M_site, arma::vec K,
                             arma::vec lambdatilde, 
                             double mu0, double n0, double pi0,
                             arma::mat u,
                             arma::mat logz, 
                             arma::mat X_w,
                             arma::mat beta_w,
                             arma::vec sigma,
                             arma::vec mu,
                             double sigma_gamma,
                             double v_sd,
                             arma::vec p11, 
                             arma::vec p10, 
                             arma::mat theta11, 
                             arma::vec theta10,
                             arma::vec emptySites){
  
  int S = y.n_slices;
  int n = M_site.size();
  
  int index_m = 0;
  int index_mj = 0;
  
  arma::mat Xw_beta = X_w * beta_w;
  
  arma::cube c_imk = arma::cube(v.n_rows, max(K), S);
  arma::mat delta = arma::mat(v.n_rows, S);
  arma::mat gamma = arma::mat(v.n_rows, S);
  
  double p0 = n0 / (n0 + mu0);
  
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        
        int currentK = K[index_m];
        
        arma::vec y_counts = arma::zeros(currentK);
        arma::vec u_im = arma::zeros(currentK);
        bool allLessThanTwo = true;
        for(int k = 0; k < currentK; k++){
          y_counts[k] = y(index_m, k, j);
          u_im[k] = u(index_m, k);
          if(y_counts[k] > 2){
            allLessThanTwo = false;
          }
        }
        
        arma::vec log_allProbs;
        arma::mat mat_delta_c_d;
        
        double v_star;
        
        if(!allLessThanTwo) {
          
          if(delta(index_m, j) == 1 | gamma(index_m, j) == 1){ // value existing
            
            v_star = v(index_m, j);
            
            if(emptySites[i] == 0){ // not an empty tube
              
              log_allProbs = arma::zeros(3 * pow(2, currentK));
              mat_delta_c_d = arma::zeros(3 * pow(2, currentK),
                                          2 + currentK);
              
              // delta = 0, gamma = 0
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta0_cpp(y_counts,
                                                                 c_imk_current,
                                                                 currentK,
                                                                 n0, p0, pi0,
                                                                 lambda[j],
                                                                 lambdatilde[j]);
                
                double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                
                double prob_d = 0;
                for(int l3 = 0; l3 < currentK; l3++){
                  prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                }
                
                log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                // log_allProbs[l] = 0;
                
                mat_delta_c_d(l, 0) = 0;
                mat_delta_c_d(l, 1) = 0;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                }
                
              }
              
              // delta = 1
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts,
                                                                     c_imk_current,
                                                                     currentK, n0, p0, pi0,
                                                                     r_nb[j],
                                                                     v_star,
                                                                     lambdatilde[j],
                                                                     lambda[j],
                                                                     u_im);
                
                double prob_delta = log(theta11(index_m,j));
                
                double log_prior_v = R::dnorm(v_star,
                                              logz(i, j) + Xw_beta(index_m, j),//r(index_m) * alpha[j],
                                              sigma[j], 1);
                
                double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                for(int k = 0; k < currentK; k++){
                  prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                }
                
                
                log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + log_prior_v + prob_c;
                // log_allProbs[pow(2, currentK) + l] = 0;
                
                mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                }
                
              }
              
              // gamma = 1
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts,
                                                                     c_imk_current,
                                                                     currentK, n0, p0, pi0,
                                                                     r_nb[j],
                                                                     v_star,
                                                                     lambdatilde[j],
                                                                     lambda[j],
                                                                     u_im);
                
                double log_prior_v = R::dnorm(v_star,
                                              mu[j],
                                              sigma_gamma, 1);
                
                double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                
                double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                for(int k = 0; k < currentK; k++){
                  prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                }
                
                log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + log_prior_v +
                  prob_c;
                // log_allProbs[2 * pow(2, currentK) + l] = 0;
                
                mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                }
                
                
              }
              
            }
            
          } else {
            
            double mean_v = log(mean(y_counts / exp(lambda[j] + u_im)));
            if(mean_v < -10){
              mean_v = - 10;
            }
            v_star = R::rnorm(mean_v, v_sd);
            
            if(emptySites[i] == 0){ // not an empty tube
              
              log_allProbs = arma::zeros(3 * pow(2, currentK));
              mat_delta_c_d = arma::zeros(3 * pow(2, currentK),
                                          2 + currentK);
              
              // delta = 0, gamma = 0
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta0_cpp(y_counts,
                                                                 c_imk_current,
                                                                 currentK,
                                                                 n0, p0, pi0,
                                                                 lambda[j],
                                                                 lambdatilde[j]);
                
                double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                
                double prob_d = 0;
                for(int l3 = 0; l3 < currentK; l3++){
                  prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                }
                
                log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                // log_allProbs[l] = 0;//log_prob_y + prob_delta + prob_d;
                
                mat_delta_c_d(l, 0) = 0;
                mat_delta_c_d(l, 1) = 0;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                }
                
              }
              
              // delta = 1
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts,
                                                                     c_imk_current,
                                                                     currentK, n0, p0, pi0,
                                                                     r_nb[j],
                                                                     v_star,
                                                                     lambdatilde[j],
                                                                     lambda[j],
                                                                     u_im);
                
                double log_prior_v = R::dnorm(v_star,
                                              logz(i, j) + Xw_beta(index_m, j),//r(index_m) * alpha[j],
                                              sigma[j], 1);
                
                double log_proposal_v = R::dnorm(v_star, mean_v, v_sd, 1);
                
                double prob_delta = log(theta11(index_m,j));
                
                double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                for(int k = 0; k < currentK; k++){
                  prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                }
                
                log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c +
                  log_prior_v - log_proposal_v;
                // log_allProbs[pow(2, currentK) + l] = 0;
                
                mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                }
                
                
              }
              
              // gamma = 1
              for(int l = 0; l < pow(2, currentK); l++){
                
                arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                
                double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts,
                                                                     c_imk_current,
                                                                     currentK, n0, p0, pi0,
                                                                     r_nb[j],
                                                                     v_star,
                                                                     lambdatilde[j],
                                                                     lambda[j],
                                                                     u_im);
                
                double log_prior_v = R::dnorm(v_star,
                                              mu[j],
                                              sigma_gamma, 1);
                
                double log_proposal_v = R::dnorm(v_star, mean_v, v_sd, 1);
                
                double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                
                double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                for(int k = 0; k < currentK; k++){
                  prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                }
                
                log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c +
                  log_prior_v - log_proposal_v;
                // log_allProbs[2 * pow(2, currentK) + l] = 0;
                
                mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                for(int k = 0; k < currentK; k++){
                  mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                }
                
                
              }
              
            }
            
          }
          
          log_allProbs = log_allProbs - max(log_allProbs);
          arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
          
          // Sample new customer assigment
          // NumericVector toSample(allProbs.size());
          arma::vec toSample(allProbs.size());
          for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
          int sampledComb_idx = sample_cpp(toSample, allProbs);//RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
          
          delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
          gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
          for(int k = 0; k < currentK; k++){
            c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
          }
          
          v(index_m, j) = v_star;
          
        } else {
          
          delta(index_m, j) = 0;
          gamma(index_m, j) = 0;
          for(int k = 0; k < currentK; k++){
            c_imk(index_m, k, j) = 0;
          }
          
        }
        
        
        index_mj += 1;
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                       _["c_imk"] = c_imk,
                       _["gamma"] = gamma,
                       _["v"] = v);
}

///// WITH PROPOSALS

// [[Rcpp::export]]
int convertDeltaIndexes(double delta, 
                        double gamma, 
                        arma::vec c,
                        int K){
  
  int index;
  
  int c_index = 0;
  
  for(int k = 0; k < K; k++){
    c_index += c[k] * pow(2,k);
  }
  
  if(delta == 0 & gamma == 0){
    
    index = c_index;
    
  } else if(delta == 1 & gamma == 0){
    
    index = c_index + pow(2, K);
    
  } else {
    
    index = c_index + 2 * pow(2, K);
    
  }
  
  return(index);
}

// [[Rcpp::export]]
arma::vec convertIndexToDeltaGammaC(int index, 
                                    int K){
  
  arma::vec indexes = arma::zeros(2 + K);
  int twoK = pow(2, K);
  int index_c = index % twoK;
  arma::vec c_vector = DecToBin_cpp(index_c, K);
  
  if(index < pow(2,K)){
    
    indexes[0] = 0;
    indexes[1] = 0;
    
  } else if(index < 2 * pow(2,K)){
    
    indexes[0] = 1;
    indexes[1] = 0;
    
  } else {
    
    indexes[0] = 0;
    indexes[1] = 1;
    
  }
  
  for(int k = 0; k < K; k++){
    
    indexes[2 + k] = c_vector[k];
    
  }
  
  return(indexes);
}

double computeProb(int delta, 
                   int gamma,
                   arma::vec c_imk,
                   arma::vec y_counts,
                   double n0, 
                   double p0,
                   double pi0,
                   double lambda_j,
                   double lambdatilde_j,
                   double theta11,
                   double theta10,
                   int currentK,
                   double p11,
                   double p10,
                   double r_nb,
                   double v_star,
                   arma::vec u_im,
                   double prior_v,
                   double sigma_j,
                   double mu,
                   double sigma_gamma){
  
  double prob;
  
  if(delta == 0 & gamma == 0){
    
    // arma::vec c_imk_new = 2 * DecToBin_cpp(l, currentK);
    arma::vec c_imk_new = 2 * c_imk;//DecToBin_cpp(l, currentK);
    
    double log_prob_y_new = compute_logprob_y_delta0_cpp(y_counts, 
                                                         c_imk_new, 
                                                         currentK, 
                                                         n0, p0, pi0,
                                                         lambda_j, 
                                                         lambdatilde_j);
    
    double prob_delta_new = log(1 - theta11) + log(1 - theta10);
    double prob_delta_current = log(1 - theta11) + log(1 - theta10);
    
    double prob_d_new = 0;
    for(int l3 = 0; l3 < currentK; l3++){
      prob_d_new += R::dbinom(c_imk_new[l3]/ 2.0, 1, p10, 1);
    }
    
    prob = log_prob_y_new + prob_delta_new + prob_d_new;
    
  } else if(delta == 1 & gamma == 0){
    
    double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts, 
                                                         c_imk, 
                                                         currentK, n0, p0, pi0, 
                                                         r_nb,
                                                         v_star, 
                                                         lambdatilde_j, 
                                                         lambda_j,
                                                         u_im);
    
    double prob_delta = log(theta11);
    
    double log_prior_v = R::dnorm(v_star, 
                                  prior_v,
                                  sigma_j, 1);
    
    double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
    for(int k = 0; k < currentK; k++){
      prob_c += R::dbinom(c_imk[k], 1, p11, 1);
    }
    
    prob = log_prob_y + prob_delta + log_prior_v + prob_c;
    
  } else {
    
    double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts, 
                                                         c_imk, 
                                                         currentK, n0, p0, pi0,
                                                         r_nb,
                                                         v_star, 
                                                         lambdatilde_j, 
                                                         lambda_j,
                                                         u_im);
    
    double log_prior_v = R::dnorm(v_star, 
                                  mu,
                                  sigma_gamma, 1);
    
    double prob_delta = log(1 - theta11) + log(theta10);
    
    double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
    for(int k = 0; k < currentK; k++){
      prob_c += R::dbinom(c_imk[k], 1, p11, 1);
    }
    
    prob = log_prob_y + prob_delta + log_prior_v + 
      prob_c; 
  }
  
  return(prob);
}

// [[Rcpp::export]]
List update_delta_c_d_proposals(arma::cube c_imk,
                                arma::mat delta,
                                arma::mat gamma,
                                arma::cube A,
                                arma::cube y, 
                                arma::mat v, 
                                arma::vec lambda,
                                arma::vec r_nb,
                                arma::vec M_site, arma::vec K,
                                arma::vec lambdatilde, 
                                double mu0, double n0, double pi0,
                                arma::mat u,
                                arma::mat logz, 
                                arma::vec r,
                                arma::vec alpha,
                                arma::vec sigma,
                                arma::vec mu,
                                double sigma_gamma,
                                double v_sd,
                                arma::vec p11, 
                                arma::vec p10, 
                                arma::mat theta11, 
                                arma::vec theta10,
                                arma::vec emptySites){
  
  int S = y.n_slices;
  int n = M_site.size();
  
  int index_m = 0;
  
  // arma::cube c_imk = arma::cube(v.n_rows, max(K), S);
  // arma::mat delta = arma::mat(v.n_rows, S);
  // arma::mat gamma = arma::mat(v.n_rows, S);
  
  double p0 = n0 / (n0 + mu0);
  
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        
        int currentK = K[index_m];
        
        arma::vec y_counts = arma::zeros(currentK);
        arma::vec u_im = arma::zeros(currentK);
        for(int k = 0; k < currentK; k++){
          y_counts[k] = y(index_m, k, j);
          u_im[k] = u(index_m, k);
        }
        
        // current value
        arma::vec c_currents = arma::zeros(currentK);
        for(int k = 0; k < currentK; k++){
          c_currents[k] = c_imk(index_m, k, j);
        }
        int delta_current = delta(index_m, j);
        int gamma_current = gamma(index_m, j);
        
        int currentIndex = convertDeltaIndexes(delta_current,
                                               gamma_current,
                                               c_currents,
                                               currentK);
        
        // propose new value 
        arma::vec proposalProbs = arma::zeros(3 * pow(2, currentK));
        for(int l = 0; l < proposalProbs.size(); l++){
          
          proposalProbs[l] = A(index_m, j, l);
          
        }
        
        NumericVector toSample(proposalProbs.size());
        for(int l5 = 0; l5 < proposalProbs.size(); l5++) toSample(l5) = l5;
        int newIndex = 0;//RcppArmadillo::sample(toSample, 1, 1, proposalProbs)[0];
        
        arma::vec newIndexes = convertIndexToDeltaGammaC(newIndex, currentK);
        int delta_new = newIndexes[0];
        int gamma_new = newIndexes[1];
        arma::vec c_new = arma::zeros(currentK);
        for(int k = 0; k < currentK; k++){
          c_new[k] = newIndexes[2 + k];
        }
        
        if(newIndex != currentIndex){
          
          double v_star;
          
          if(delta_current == 1 | gamma_current == 1){ // value existing
            
            v_star = v(index_m, j);
            
          } else {
            
            double mean_v = log(mean(y_counts / exp(lambda[j] + u_im)));
            if(mean_v < -10){
              mean_v = - 10;
            }
            v_star = R::rnorm(mean_v, v_sd);
            
          }
          
          double log_prob_current = computeProb(delta_current, 
                                                gamma_current,
                                                c_currents,
                                                y_counts,
                                                n0, 
                                                p0,
                                                pi0,
                                                lambda[j],
                                                lambdatilde[j],
                                                theta11(index_m,j),
                                                theta10[j],
                                                currentK,
                                                p11[j],
                                                p10[j],
                                                r_nb[j],
                                                v_star,
                                                u_im,
                                                logz(i, j) + r(index_m) * alpha[j],
                                                sigma[j],
                                                mu[j],
                                                sigma_gamma);
          
          double log_prob_new = computeProb(delta_new, 
                                            gamma_new,
                                            c_new,
                                            y_counts,
                                            n0, 
                                            p0,
                                            pi0,
                                            lambda[j],
                                            lambdatilde[j],
                                            theta11(index_m,j),
                                            theta10[j],
                                            currentK,
                                            p11[j],
                                            p10[j],
                                            r_nb[j],
                                            v_star,
                                            u_im,
                                            logz(i, j) + r(index_m) * alpha[j],
                                            sigma[j],
                                            mu[j],
                                            sigma_gamma);
          
          double prop_current_to_new = A(index_m, j, currentIndex);
          double prop_new_to_current = A(index_m, j, newIndex);
          
          double proposal_ratio = prop_current_to_new / prop_new_to_current;
          
          double log_posterior_ratio =  log_prob_new - log_prob_current;
          
          if(R::runif(0, 1) < exp(log_posterior_ratio) * proposal_ratio){
            
            delta(index_m, j) = delta_new;
            gamma(index_m, j) = gamma_new;
            for(int k = 0; k < currentK; k++){
              c_imk(index_m, k, j) = c_new[k];
            }
            
            v(index_m, j) = v_star;
            
          }
          
          
        }
        
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                       _["c_imk"] = c_imk,
                       _["gamma"] = gamma,
                       _["v"] = v);
}

/// PROBABILITIES

// [[Rcpp::export]]
arma::vec update_theta10_cpp(arma::vec theta_10, arma::mat delta, arma::mat gamma,
                             arma::vec M_site, arma::vec K, arma::vec emptySites,
                             double a0, double b0){
  
  int n = M_site.size();
  int S = theta_10.size();
  
  for(int j = 0; j < S; j++){
    
    int numPresents = 0;
    int numSuccesses = 0;
    
    int l = 0;
    
    for (int i = 0; i < n; i++) {
      
      for (int m = 0; m < M_site[i]; m++) {
        
        int currentK = K[l];
        
        if(emptySites[i] == 0){
          
          if(delta(l, j) == 0){
            
            numPresents++;
            
            if(gamma(l, j) == 1){
              
              numSuccesses++;
              
            }
            
            
          }
          
        }
        
        l++;
      }
      
    }
    
    theta_10[j] = R::rbeta(a0 + numSuccesses, b0 + numPresents - numSuccesses);
    
  }
  
  return theta_10;
}

// [[Rcpp::export]]
arma::vec update_p_11_cpp(arma::vec p_11, arma::mat delta, arma::mat gamma,
                          arma::cube c_imk, arma::vec M_site, arma::vec K,
                          double a0, double b0){
  
  int n = M_site.size();
  int S = c_imk.n_slices;
  
  for(int j = 0; j < S; j++){
    
    int numPresents = 0;
    int numSuccesses = 0;
    
    int l = 0;
    
    for (int i = 0; i < n; i++) {
      
      for (int m = 0; m < M_site[i]; m++) {
        
        int currentK = K[l];
        
        if(delta(l,j) == 1 | gamma(l,j) == 1){
          
          for (int k = 0; k < currentK; k++) {
            
            if(c_imk(l,k,j) == 1){
              
              numSuccesses++;
              
            }
            
            
            numPresents++;
            
          }
          
        }
        
        l++;
      }
      
    }
    
    p_11[j] = R::rbeta(a0 + numSuccesses, b0 + numPresents - numSuccesses);
    
  }
  
  return p_11;
}

// [[Rcpp::export]]
arma::vec update_p_10_cpp(arma::vec p_10, arma::mat delta, arma::mat gamma, 
                          arma::cube c_imk, arma::vec M_site, arma::vec K,
                          double a_p1, double b_p1){
  
  int n = M_site.size();
  int S = c_imk.n_slices;
  
  for(int j = 0; j < S; j++){
    
    int numPresents = 0;
    int numSuccesses = 0;
    
    int l = 0;
    
    for (int i = 0; i < n; i++) {
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(delta(l, j) == 0 & gamma(l, j) == 0){
          
          int currentK = K[l];
          
          for (int k = 0; k < currentK; k++) {
            
            if(c_imk(l,k,j) == 2){
              
              numSuccesses++;
              
            }
            
            numPresents++;
            
          }
          
        }
        
        l++;
      }
      
    }
    
    p_10[j] = R::rbeta(a_p1 + numSuccesses, b_p1 + numPresents - numSuccesses);
    
  }
  
  return p_10;
}

///// WITH POISSON-GAMMA VARIABLES

// // [[Rcpp::export]]
// double compute_logprob_y_delta1_pg(arma::vec y_counts, arma::vec c_imk_current, 
                                      //                                    arma::vec lambda_k,
                                      //                                    int currentK, double n0, double p0, 
                                      //                                    double lambdatilde,
                                      //                                    double lambda){
  //   
    //   double sum = 0;
    //   for(int k = 0; k < currentK; k++){
      //     if(c_imk_current[k] == 1){
        //       sum += R::dpois(y_counts[k], lambda_k[k], 1);
        //       // double pi = exp(lambda + v_im + u_im[k]) / (exp(lambda + v_im + u_im[k]) + r_nb);
        //       // sum += R::dnbinom(y_counts[k], r_nb, 1 - pi, 1);
        //     } else if(c_imk_current[k] == 2){
          //       sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
          //     } else {
            //       sum += R::dnbinom(y_counts[k], n0, p0, 1);
            //     }
      //   }
    //   
      //   return(sum);
    // }
// 
  // // [[Rcpp::export]]
// double compute_logprob_y_delta0_pg(arma::vec y_counts, arma::vec c_imk_current, 
                                      //                                     int currentK, double n0, double p0, 
                                      //                                     double lambda, double lambdatilde){
  //   
    //   double sum = 0;
    //   for(int k = 0; k < currentK; k++){
      //     if(c_imk_current[k] == 2){
        //       sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
        //     } else {
          //       sum += R::dnbinom(y_counts[k], n0, p0, 1);
          //     }
      //   }
    //   
      //   return(sum);
    // }
// 
  // // [[Rcpp::export]]
// List update_delta_c_d_cpp_pg(arma::cube y, 
                                //                              arma::cube lambda_ijk,
                                //                              arma::vec lambda,
                                //                              arma::vec r_nb,
                                //                              arma::vec M_site, arma::vec K,
                                //                              arma::vec lambdatilde, 
                                //                              double mu0, double n0, 
                                //                              arma::mat u, arma::vec p11, 
                                //                              arma::vec p10, arma::mat theta11, arma::vec theta10,
                                //                              arma::vec emptySites){
  //   
    //   int S = y.n_slices;
    //   int n = M_site.size();
    //   
      //   int index_m = 0;
      //   int index_mj = 0;
      //   
        //   arma::cube c_imk = arma::cube(lambda_ijk.n_rows, max(K), S);
        //   arma::mat delta = arma::mat(lambda_ijk.n_rows, S);
        //   arma::mat gamma = arma::mat(lambda_ijk.n_rows, S);
        //   
          //   double p0 = n0 / (n0 + mu0);
          //   
            //   for(int i = 0; i < n; i++){
              //     for(int m = 0; m < M_site[i]; m++){
                //       for(int j = 0; j < S; j++){
                  //         
                    //         int currentK = K[index_m];
                    //         
                      //         arma::vec y_counts = arma::zeros(currentK);
                      //         arma::vec u_im = arma::zeros(currentK);
                      //         for(int k = 0; k < currentK; k++){
                        //           y_counts[k] = y(index_m, k, j);
                        //           u_im[k] = u(index_m, k);
                        //         }
                      //         
                        //         arma::vec lambda_k = arma::zeros(currentK);
                        //         for(int k = 0 ; k < currentK; k++){
                          //           lambda_k[k] = arma::as_scalar(lambda_ijk(index_m, k, j));
                          //         }
                        //         
                          //         // arma::conv_to<arma::vec>::from(lambda_ijk.subcube(arma::span(index_m), 
                                                                                          //         //                                                                        arma::span(), 
                                                                                          //         //                                                                        arma::span(j)));
                        //         
                          //         arma::vec log_allProbs;
                        //         arma::mat mat_delta_c_d;
                        //         
                          //         if(emptySites[i] == 0){ // not an empty tube
                            //           
                              //           log_allProbs = arma::zeros(3 * pow(2, currentK));
                              //           mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                                                       //                                       2 + currentK);
                              //           
                                //           // delta = 0, gamma = 0
                                //           for(int l = 0; l < pow(2, currentK); l++){
                                  //             
                                    //             arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                                    //             
                                      //             double log_prob_y = compute_logprob_y_delta0_pg(y_counts, 
                                                                                                     //                                                              c_imk_current, 
                                                                                                     //                                                              currentK, 
                                                                                                     //                                                              n0, p0, 
                                                                                                     //                                                              lambda[j], 
                                                                                                     //                                                                    lambdatilde[j]);
                                      //             
                                        //             double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                                        //             
                                          //             double prob_d = 0;
                                          //             for(int l3 = 0; l3 < currentK; l3++){
                                            //               prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                                            //             }
                                          //             
                                            //             log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                                            //             
                                              //             mat_delta_c_d(l, 0) = 0;
                                              //             mat_delta_c_d(l, 1) = 0;
                                              //             for(int k = 0; k < currentK; k++){
                                                //               mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                                                //             }
                                              //             
                                                //           }
                                //           
                                  //           // delta = 1
                                  //           for(int l = 0; l < pow(2, currentK); l++){
                                    //             
                                      //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                      //             
                                        //             double log_prob_y = compute_logprob_y_delta1_pg(y_counts, 
                                                                                                       //                                                             c_imk_current, 
                                                                                                       //                                                             lambda_k,
                                                                                                       //                                                               currentK, n0, p0, 
                                                                                                       //                                                               lambdatilde[j], 
                                                                                                       //                                                                          lambda[j]);
                                        //             
                                          //             double prob_delta = log(theta11(index_m,j));
                                          //             
                                            //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                            //             for(int k = 0; k < currentK; k++){
                                              //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                              //             }
                                            //             
                                              //             mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                                              //             mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                                              //             for(int k = 0; k < currentK; k++){
                                                //               mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                //             }
                                              //             
                                                //             log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                //             
                                                  //           }
                                  //           
                                    //           // gamma = 1
                                    //           for(int l = 0; l < pow(2, currentK); l++){
                                      //             
                                        //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                        //             
                                          //             double log_prob_y = compute_logprob_y_delta1_pg(y_counts, 
                                                                                                         //                                                             c_imk_current, 
                                                                                                         //                                                             lambda_k,
                                                                                                         //                                                             currentK, n0, p0, 
                                                                                                         //                                                             lambdatilde[j], 
                                                                                                         //                                                                        lambda[j]);
                                          //             
                                            //             double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                                            //             
                                              //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                              //             for(int k = 0; k < currentK; k++){
                                                //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                                //             }
                                              //             
                                                //             mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                                                //             mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                                                //             for(int k = 0; k < currentK; k++){
                                                  //               mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                  //             }
                                                //             
                                                  //             log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                  //             
                                                    //           }
                                    //           
                                      //         } 
                        //         
                          //         log_allProbs = log_allProbs - max(log_allProbs);
                          //         arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
                          //         
                            //         // Sample new customer assigment
                          //         NumericVector toSample(allProbs.size());
                          //         for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
                          //         int sampledComb_idx = 0;//RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
                          //         
                            //         delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
                            //         gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
                            //         for(int k = 0; k < currentK; k++){
                              //           c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
                              //         }
                            //         
                              //         // // copy output
                            //         // for(int l3 = 0; l3 < allProbs.size(); l3++){
                              //         //   all_probs(index_m, j, l3) = allProbs(l3);
                              //         //   for(int l2 = 0; l2 < (2 + currentK); l2++){
                                //         //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
                                //         //   }
                              //         // }
                            //         
                              //         index_mj += 1;
                              //       }
                //       
                  //       index_m += 1;
                  //     }
              //   }
          //   
            //   return List::create(_["delta"] = delta,
                                      //                       _["c_imk"] = c_imk,
                                      //                       _["gamma"] = gamma);
          // }

// OLD STUFF

// // [[Rcpp::export]]
// List update_delta_c_d_cpp_nb(arma::cube y, arma::mat v, 
                                //                              arma::vec lambda,
                                //                              arma::vec r_nb,
                                //                              arma::vec M_site, arma::vec K,
                                //                              arma::vec lambdatilde, 
                                //                              double mu0, double n0, 
                                //                              arma::mat u, arma::vec p11, 
                                //                              arma::vec p10, arma::mat theta11, arma::vec theta10,
                                //                              arma::vec emptySites){
  //   
    //   int S = y.n_slices;
    //   int n = M_site.size();
    //   
      //   int index_m = 0;
      //   int index_mj = 0;
      //   
        //   arma::cube c_imk = arma::cube(v.n_rows, max(K), S);
        //   arma::mat delta = arma::mat(v.n_rows, S);
        //   arma::mat gamma = arma::mat(v.n_rows, S);
        //   
          //   double p0 = n0 / (n0 + mu0);
          //   
            //   for(int i = 0; i < n; i++){
              //     for(int m = 0; m < M_site[i]; m++){
                //       for(int j = 0; j < S; j++){
                  //         
                    //         int currentK = K[index_m];
                    //         
                      //         arma::vec y_counts = arma::zeros(currentK);
                      //         arma::vec u_im = arma::zeros(currentK);
                      //         for(int k = 0; k < currentK; k++){
                        //           y_counts[k] = y(index_m, k, j);
                        //           u_im[k] = u(index_m, k);
                        //         }
                      //         
                        //         double v_im = v(index_m, j);
                        //         
                          //         arma::vec log_allProbs;
                        //         arma::mat mat_delta_c_d;
                        //         
                          //         if(emptySites[i] == 0){ // not an empty tube
                            //           
                              //           log_allProbs = arma::zeros(3 * pow(2, currentK));
                              //           mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                                                       //                                       2 + currentK);
                              //           
                                //           // delta = 0, gamma = 0
                                //           for(int l = 0; l < pow(2, currentK); l++){
                                  //             
                                    //             arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                                    //             
                                      //             double log_prob_y = compute_logprob_y_delta0_cpp(y_counts, 
                                                                                                      //                                                              c_imk_current, 
                                                                                                      //                                                              currentK, 
                                                                                                      //                                                              n0, p0, 
                                                                                                      //                                                              lambda[j], 
                                                                                                      //                                                                    lambdatilde[j]);
                                      //             
                                        //             double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                                        //             
                                          //             double prob_d = 0;
                                          //             for(int l3 = 0; l3 < currentK; l3++){
                                            //               prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                                            //             }
                                          //             
                                            //             log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                                            //             
                                              //             mat_delta_c_d(l, 0) = 0;
                                              //             mat_delta_c_d(l, 1) = 0;
                                              //             for(int k = 0; k < currentK; k++){
                                                //               mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                                                //             }
                                              //             
                                                //           }
                                //           
                                  //           // delta = 1
                                  //           for(int l = 0; l < pow(2, currentK); l++){
                                    //             
                                      //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                      //             
                                        //             double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts, 
                                                                                                            //                                                                  c_imk_current, 
                                                                                                            //                                                                  currentK, n0, p0, 
                                                                                                            //                                                                  r_nb[j],
                                                                                                            //                                                                      v_im, 
                                                                                                            //                                                                      lambdatilde[j], 
                                                                                                            //                                                                                 lambda[j],
                                                                                                            //                                                                                       u_im);
                                        //             
                                          //             double prob_delta = log(theta11(index_m,j));
                                          //             
                                            //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                            //             for(int k = 0; k < currentK; k++){
                                              //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                              //             }
                                            //             
                                              //             mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                                              //             mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                                              //             for(int k = 0; k < currentK; k++){
                                                //               mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                //             }
                                              //             
                                                //             log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                //             
                                                  //           }
                                  //           
                                    //           // gamma = 1
                                    //           for(int l = 0; l < pow(2, currentK); l++){
                                      //             
                                        //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                        //             
                                          //             double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts, 
                                                                                                              //                                                                  c_imk_current, 
                                                                                                              //                                                                  currentK, n0, p0, 
                                                                                                              //                                                                  r_nb[j],
                                                                                                              //                                                                      v_im, 
                                                                                                              //                                                                      lambdatilde[j], 
                                                                                                              //                                                                                 lambda[j],
                                                                                                              //                                                                                       u_im);
                                          //             
                                            //             double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                                            //             
                                              //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                              //             for(int k = 0; k < currentK; k++){
                                                //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                                //             }
                                              //             
                                                //             mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                                                //             mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                                                //             for(int k = 0; k < currentK; k++){
                                                  //               mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                  //             }
                                                //             
                                                  //             log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                  //             
                                                    //           }
                                    //           
                                      //         } 
                        //         
                          //         log_allProbs = log_allProbs - max(log_allProbs);
                          //         arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
                          //         
                            //         // Sample new customer assigment
                          //         NumericVector toSample(allProbs.size());
                          //         for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
                          //         int sampledComb_idx = 0;//RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
                          //         
                            //         delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
                            //         gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
                            //         for(int k = 0; k < currentK; k++){
                              //           c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
                              //         }
                            //         
                              //         // // copy output
                            //         // for(int l3 = 0; l3 < allProbs.size(); l3++){
                              //         //   all_probs(index_m, j, l3) = allProbs(l3);
                              //         //   for(int l2 = 0; l2 < (2 + currentK); l2++){
                                //         //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
                                //         //   }
                              //         // }
                            //         
                              //         index_mj += 1;
                              //       }
                //       
                  //       index_m += 1;
                  //     }
              //   }
          //   
            //   return List::create(_["delta"] = delta,
                                      //                       _["c_imk"] = c_imk,
                                      //                       _["gamma"] = gamma);
          // }
// 
  // 
  // // [[Rcpp::export]]
// List update_delta_c_d_cpp_existingprobs(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                                           //                                         arma::vec M_site, arma::vec K,
                                           //                                         arma::cube &allProbs, arma::cube &allCombs, 
                                           //                                         arma::vec emptySites){
  //   
    //   int S = c_imk.n_slices;
    //   int n = M_site.size();
    //   
      //   int index_m = 0;
      //   int index_mj = 0;
      //   
        //   for(int i = 0; i < n; i++){
          //     for(int m = 0; m < M_site[i]; m++){
            //       for(int j = 0; j < S; j++){
              //         
                //         int currentK = K[index_m];
                //         
                  //         arma::vec currentAllProbs = arma::zeros(3 * pow(2, currentK));
                  //         arma::mat mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                                                   //                                               2 + currentK);
                  //         
                    //         for(int l = 0; l < currentAllProbs.size(); l++){
                      //           currentAllProbs[l] = allProbs(index_m, j, l);
                      //           for(int l2 = 0; l2 < (2 + currentK); l2++){
                        //             mat_delta_c_d(l, l2) = allCombs(index_mj, l, l2);
                        //           }
                      //         }
                  //         
                    //         // Sample new customer assigment
                  //         NumericVector toSample(currentAllProbs.size());
                  //         for(int l5 = 0; l5 < currentAllProbs.size(); l5++) toSample(l5) = l5;
                  //         int sampledComb_idx = 0;//RcppArmadillo::sample(toSample, 1, 1, currentAllProbs)[0];
                  //         
                    //         delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
                    //         gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
                    //         for(int k = 0; k < currentK; k++){
                      //           c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
                      //         }
                    //         
                      //         index_mj += 1;
                      //       }
            //       
              //       index_m += 1;
              //     }
          //   }
      //   
        //   return List::create(_["delta"] = delta,
                                  //                       _["c_imk"] = c_imk,
                                  //                       _["gamma"] = gamma);
      // }
// 
  // 
  // List update_delta_c_d_cpp_old(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                                   //                               arma::cube y, arma::mat v_bar, arma::vec lambda,
                                   //                               arma::mat X_z, arma::mat beta_z,
                                   //                               arma::mat r, arma::mat alpha,
                                   //                               arma::vec M_site, arma::vec K,
                                   //                               arma::vec lambdatilde, 
                                   //                               double mu0, double n0, 
                                   //                               arma::mat u, arma::vec p11, 
                                   //                               arma::vec p10, arma::mat theta11, arma::vec theta10,
                                   //                               arma::vec emptySites){
    //   
      //   int S = y.n_slices;
      //   int n = M_site.size();
      //   // arma::cube all_probs = arma::cube(sum(M_site), S, 3 * pow(2, max(K)));
      //   // arma::cube all_combinations = arma::cube(sum(M_site) * S, 3 * pow(2, max(K)), 2 + max(K));
      //   
        //   arma::mat Xzbeta = X_z * beta_z;
        //   arma::mat ralpha = r * alpha;
        //   
          //   int index_m = 0;
          //   int index_mj = 0;
          //   
            //   double p0 = n0 / (n0 + mu0);
            //   
              //   for(int i = 0; i < n; i++){
                //     for(int m = 0; m < M_site[i]; m++){
                  //       for(int j = 0; j < S; j++){
                    //         
                      //         int currentK = K[index_m];
                      //         
                        //         arma::vec y_counts = arma::zeros(currentK);
                        //         arma::vec u_im = arma::zeros(currentK);
                        //         for(int k = 0; k < currentK; k++){
                          //           y_counts[k] = y(index_m, k, j);
                          //           u_im[k] = u(index_m, k);
                          //         }
                        //         
                          //         double v_im = v_bar(index_m, j) + 
                            //           Xzbeta(i, j) +
                            //           ralpha(index_m, j);
                          //         
                            //         arma::vec log_allProbs;
                          //         arma::mat mat_delta_c_d;
                          //         
                            //         if(emptySites[i] == 0){ // not an empty tube
                              //           
                                //           log_allProbs = arma::zeros(3 * pow(2, currentK));
                                //           mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                                                         //                                       2 + currentK);
                                //           
                                  //           // delta = 0, gamma = 0
                                  //           for(int l = 0; l < pow(2, currentK); l++){
                                    //             
                                      //             arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                                      //             
                                        //             double log_prob_y = compute_logprob_y_delta0_cpp(y_counts, 
                                                                                                        //                                                              c_imk_current, 
                                                                                                        //                                                              currentK, 
                                                                                                        //                                                              n0, p0, 
                                                                                                        //                                                              lambda[j], 
                                                                                                        //                                                                    lambdatilde[j]);
                                        //             
                                          //             double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                                          //             
                                            //             double prob_d = 0;
                                            //             for(int l3 = 0; l3 < currentK; l3++){
                                              //               prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                                              //             }
                                            //             
                                              //             log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                                              //             
                                                //             mat_delta_c_d(l, 0) = 0;
                                                //             mat_delta_c_d(l, 1) = 0;
                                                //             for(int k = 0; k < currentK; k++){
                                                  //               mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                                                  //             }
                                                //             
                                                  //           }
                                  //           
                                    //           // delta = 1
                                    //           for(int l = 0; l < pow(2, currentK); l++){
                                      //             
                                        //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                        //             
                                          //             double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                                                                          //                                                              c_imk_current, 
                                                                                                          //                                                              currentK, n0, p0, 
                                                                                                          //                                                              v_im, 
                                                                                                          //                                                              lambdatilde[j], 
                                                                                                          //                                                                         lambda[j],
                                                                                                          //                                                                               u_im);
                                          //             
                                            //             double prob_delta = log(theta11(index_m,j));
                                            //             
                                              //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                              //             for(int k = 0; k < currentK; k++){
                                                //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                                //             }
                                              //             
                                                //             mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                                                //             mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                                                //             for(int k = 0; k < currentK; k++){
                                                  //               mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                  //             }
                                                //             
                                                  //             log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                  //             
                                                    //           }
                                    //           
                                      //           // gamma = 1
                                      //           for(int l = 0; l < pow(2, currentK); l++){
                                        //             
                                          //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                          //             
                                            //             double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                                                                            //                                                              c_imk_current, 
                                                                                                            //                                                              currentK, n0, p0, 
                                                                                                            //                                                              v_im, 
                                                                                                            //                                                              lambdatilde[j], 
                                                                                                            //                                                                         lambda[j],
                                                                                                            //                                                                               u_im);
                                            //             
                                              //             double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                                              //             
                                                //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                                //             for(int k = 0; k < currentK; k++){
                                                  //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                                  //             }
                                                //             
                                                  //             mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                                                  //             mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                                                  //             for(int k = 0; k < currentK; k++){
                                                    //               mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                    //             }
                                                  //             
                                                    //             log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                    //             
                                                      //           }
                                      //           
                                        //         } 
                          //         
                            //         log_allProbs = log_allProbs - max(log_allProbs);
                            //         arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
                            //         
                              //         // Sample new customer assigment
                            //         NumericVector toSample(allProbs.size());
                            //         for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
                            //         int sampledComb_idx = sample_cpp(toSample, allProbs);//RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
                            //         
                              //         delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
                              //         gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
                              //         for(int k = 0; k < currentK; k++){
                                //           c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
                                //         }
                              //         
                                //         // // copy output
                              //         // for(int l3 = 0; l3 < allProbs.size(); l3++){
                                //         //   all_probs(index_m, j, l3) = allProbs(l3);
                                //         //   for(int l2 = 0; l2 < (2 + currentK); l2++){
                                  //         //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
                                  //         //   }
                                //         // }
                              //         
                                //         index_mj += 1;
                                //       }
                  //       
                    //       index_m += 1;
                    //     }
                //   }
            //   
              //   return List::create(_["delta"] = delta,
                                        //                       _["c_imk"] = c_imk,
                                        //                       _["gamma"] = gamma);//,
            //   // _["all_probs"] = all_probs,
            //   // _["all_combinations"] = all_combinations);
// }
// 
  // // [[Rcpp::export]]
// List update_delta_c_d_cpp(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                             //                           arma::cube y, arma::mat v, 
                             //                           arma::vec lambda,
                             //                           arma::vec M_site, arma::vec K,
                             //                           arma::vec lambdatilde, 
                             //                           double mu0, double n0, 
                             //                           arma::mat u, arma::vec p11, 
                             //                           arma::vec p10, arma::mat theta11, arma::vec theta10,
                             //                           arma::vec emptySites){
  //   
    //   int S = y.n_slices;
    //   int n = M_site.size();
    //   
      //   int index_m = 0;
      //   int index_mj = 0;
      //   
        //   double p0 = n0 / (n0 + mu0);
        //   
          //   for(int i = 0; i < n; i++){
            //     for(int m = 0; m < M_site[i]; m++){
              //       for(int j = 0; j < S; j++){
                //         
                  //         int currentK = K[index_m];
                  //         
                    //         arma::vec y_counts = arma::zeros(currentK);
                    //         arma::vec u_im = arma::zeros(currentK);
                    //         for(int k = 0; k < currentK; k++){
                      //           y_counts[k] = y(index_m, k, j);
                      //           u_im[k] = u(index_m, k);
                      //         }
                    //         
                      //         double v_im = v(index_m, j);
                      //         
                        //         arma::vec log_allProbs;
                      //         arma::mat mat_delta_c_d;
                      //         
                        //         if(emptySites[i] == 0){ // not an empty tube
                          //           
                            //           log_allProbs = arma::zeros(3 * pow(2, currentK));
                            //           mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                                                     //                                       2 + currentK);
                            //           
                              //           // delta = 0, gamma = 0
                              //           for(int l = 0; l < pow(2, currentK); l++){
                                //             
                                  //             arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
                                  //             
                                    //             double log_prob_y = compute_logprob_y_delta0_cpp(y_counts, 
                                                                                                    //                                                              c_imk_current, 
                                                                                                    //                                                              currentK, 
                                                                                                    //                                                              n0, p0, 
                                                                                                    //                                                              lambda[j], 
                                                                                                    //                                                                    lambdatilde[j]);
                                    //             
                                      //             double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
                                      //             
                                        //             double prob_d = 0;
                                        //             for(int l3 = 0; l3 < currentK; l3++){
                                          //               prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
                                          //             }
                                        //             
                                          //             log_allProbs[l] = log_prob_y + prob_delta + prob_d;
                                          //             
                                            //             mat_delta_c_d(l, 0) = 0;
                                            //             mat_delta_c_d(l, 1) = 0;
                                            //             for(int k = 0; k < currentK; k++){
                                              //               mat_delta_c_d(l, 2 + k) = c_imk_current[k];
                                              //             }
                                            //             
                                              //           }
                              //           
                                //           // delta = 1
                                //           for(int l = 0; l < pow(2, currentK); l++){
                                  //             
                                    //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                    //             
                                      //             double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                                                                      //                                                              c_imk_current, 
                                                                                                      //                                                              currentK, n0, p0, 
                                                                                                      //                                                              v_im, 
                                                                                                      //                                                              lambdatilde[j], 
                                                                                                      //                                                                         lambda[j],
                                                                                                      //                                                                               u_im);
                                      //             
                                        //             double prob_delta = log(theta11(index_m,j));
                                        //             
                                          //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                          //             for(int k = 0; k < currentK; k++){
                                            //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                            //             }
                                          //             
                                            //             mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
                                            //             mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
                                            //             for(int k = 0; k < currentK; k++){
                                              //               mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                              //             }
                                            //             
                                              //             log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                              //             
                                                //           }
                                //           
                                  //           // gamma = 1
                                  //           for(int l = 0; l < pow(2, currentK); l++){
                                    //             
                                      //             arma::vec c_imk_current = DecToBin_cpp(l, currentK);
                                      //             
                                        //             double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                                                                        //                                                              c_imk_current, 
                                                                                                        //                                                              currentK, n0, p0, 
                                                                                                        //                                                              v_im, 
                                                                                                        //                                                              lambdatilde[j], 
                                                                                                        //                                                                         lambda[j],
                                                                                                        //                                                                               u_im);
                                        //             
                                          //             double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
                                          //             
                                            //             double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
                                            //             for(int k = 0; k < currentK; k++){
                                              //               prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
                                              //             }
                                            //             
                                              //             mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
                                              //             mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
                                              //             for(int k = 0; k < currentK; k++){
                                                //               mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
                                                //             }
                                              //             
                                                //             log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
                                                //             
                                                  //           }
                                  //           
                                    //         } 
                      //         
                        //         log_allProbs = log_allProbs - max(log_allProbs);
                        //         arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
                        //         
                          //         // Sample new customer assigment
                        //         NumericVector toSample(allProbs.size());
                        //         for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
                        //         int sampledComb_idx = 0;//RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
                        //         
                          //         delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
                          //         gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
                          //         for(int k = 0; k < currentK; k++){
                            //           c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
                            //         }
                          //         
                            //         // // copy output
                          //         // for(int l3 = 0; l3 < allProbs.size(); l3++){
                            //         //   all_probs(index_m, j, l3) = allProbs(l3);
                            //         //   for(int l2 = 0; l2 < (2 + currentK); l2++){
                              //         //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
                              //         //   }
                            //         // }
                          //         
                            //         index_mj += 1;
                            //       }
              //       
                //       index_m += 1;
                //     }
            //   }
        //   
          //   return List::create(_["delta"] = delta,
                                    //                       _["c_imk"] = c_imk,
                                    //                       _["gamma"] = gamma);
        // }

