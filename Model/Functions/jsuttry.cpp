

// [[Rcpp::export]]
double logposterior_logz_cpp(arma::vec v_samples, double sigma,
                             double logz_current, double logz_star,
                             double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < v_samples.size(); l++){
    
    loglikelihood += (R::dnorm(v_samples[l], logz_star, sigma, 1) - 
      R::dnorm(v_samples[l], logz_current, sigma, 1));
    
  }
  
  double logprior = R::dnorm(logz_star, prior_mean, sqrt(prior_var), 1) - 
    R::dnorm(logz_current, prior_mean, sqrt(prior_var), 1);
  // Rcout << loglikelihood << " - " << logprior << std::endl;
  return(loglikelihood + logprior);
  
  // sum(dpois(PCR_counts, lambda = lambda_j * u_im * wstar, log = T)) -
  //   sum(dpois(PCR_counts, lambda = lambda_j * u_im * wcurrent, log = T)) + 
  //   dnorm(w_star, logz, sigmasq, log = T) - 
  //   dnorm(wcurrent, logz, sigmasq, log = T)
}

// [[Rcpp::export]]
double logposterior_logz_logistic_cpp(arma::vec y_all, double x_all,
                                      arma::vec v_all,
                                      double logz_current, 
                                      double logz_star){
  
  double logzstar_x = x_all * exp(logz_star);
  double logzcurrent_x = x_all * exp(logz_current);
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < y_all.size(); l++){
    
    double p_current = 1 / (1 + exp(- logzcurrent_x - v_all[l]));
    double p_star = 1 / (1 + exp(- logzstar_x - v_all[l]));
    
    loglikelihood += (R::dbinom(y_all[l], 1, p_star, 1) - 
      R::dbinom(y_all[l], 1, p_current, 1));
    
  }
  
  return(loglikelihood);
  
}

// [[Rcpp::export]]
double h_f_cpp(double l,
               double x_all,
               arma::vec y_all,
               arma::vec v_all,
               arma::vec v_samples,
               double tauj,
               double sigma_j,
               double prior_mean){
  
  double x_all_l = exp(l) * x_all;
  
  double term1 = - 1 / pow(sigma_j,2) * v_samples.size();// + 
  
  // double term2_asscalar = 
  
  double term2_1 = sum(y_all * x_all_l);
  double term2_2 = sum(x_all * x_all * exp(2 * (x_all_l + v_all + l)) /
                       ((1 + exp(v_all + x_all_l)) % (1 + exp(v_all + x_all_l))) );
  double term2_3 = sum(x_all * (x_all_l + 1) * exp(x_all_l + v_all + l) / (1 + exp(v_all + x_all_l)));
  
  double term2 = term2_1 + term2_2 - term2_3;
  
  // double term2 = -sum(x_all * x_all * exp(v_all + x_all_l) /
  //                     ((1 + exp(v_all + x_all_l)) % (1 + exp(v_all + x_all_l)) ) );
  double term3 = - 1 / pow(tauj,2);
  // -sum(x_all * exp(v_all + x_all_l)) + sum(y_all * x_all) -
  // sum(y_all * (v_all + x_all_l)) - sum(log(1 + exp(x_all_l + v_all))) - 
  // Rcout << term2_1 << " - " << term2_2 << " - " << term2_3 << std::endl;
  return(term1 + term2 + term3);
  
}


// [[Rcpp::export]]
double logf_cpp(double l,
                double x_all,
                double sigmaj,
                arma::vec v_samples,
                arma::vec v,
                arma::vec y,
                double tauj,
                double prior_mean){
  
  double x_all_l = exp(l) * x_all;
  
  double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
    sum(- (v_samples - l));
  
  // double loglikelihood_v =  sum(y * x_all) - 
  // sum(x_all * exp(v + x_all_l) / (1 + exp(v + x_all_l)));
  
  double loglikelihood_v = 0;
  for(int i = 0; i < y.size(); i++){
    // loglikelihood_v += y[i] * x_all - 
    // x_all * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
    loglikelihood_v += y[i] * x_all - 
      x_all * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
  }
  // sum(y * x_all) - 
  // sum(x_all * exp(v + x_all_l) / (1 + exp(v + x_all_l)));
  
  double logprior = - 1 / (pow(tauj,2)) * (l - prior_mean); 
  
  return( loglikelihood + loglikelihood_v + logprior);
}


// [[Rcpp::export]]
List update_logz_cpp(arma::mat logz, arma::vec beta0,
                     arma::mat X_z, arma::mat beta_z,
                     arma::vec mu, arma::mat v,
                     arma::vec lambda, arma::mat beta_theta,
                     arma::vec r,
                     arma::vec alpha,
                     arma::vec tau, arma::mat delta,
                     arma::mat gamma,
                     arma::vec sigma, arma::vec M_site,
                     arma::vec emptySites, double sigma_prop){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  int n = M_site.size();
  int S = beta_bar.size();
  
  arma::mat Xz_beta = X_z * beta_z;
  
  // update parameters
  
  int sum_m = 0;
  for(int i = 0; i < n; i++){
    if(!(emptySites[i] == 1)){
      for(int j = 0; j < S; j++){
        
        double logzbar_current = logz_bar(i, j);
        // double logzbar_star = R::rnorm(logz_bar(i, j), sigma_prop);
        
        arma::vec v_samples = arma::zeros(M_site[i]);
        
        int l2 = 0;
        for(int m = 0; m < M_site[i]; m++){
          if(delta(sum_m + m, j) == 1){
            v_samples[l2] = v_bar(m + sum_m, j) - r[m + sum_m] * alpha[j];
            l2 += 1;
          }
        }
        
        arma::vec v_samples2 = arma::zeros(l2);
        for(int l = 0; l < l2; l++){
          v_samples2[l] = v_samples[l];
        }
        
        arma::vec y_all = arma::zeros(M_site[i]);
        arma::vec v_all = arma::zeros(M_site[i]);
        double x_all = beta_theta(j, 1) / exp(lambda[j]);
        for(int m = 0; m < M_site[i]; m++){
          
          y_all[m] = delta(sum_m + m, j);
          v_all[m] = beta_theta(j, 0) +
            r[m + sum_m] * beta_theta(j, 2);
          
        }
        
        double prior_mean = beta_bar[j] + Xz_beta(i, j);
        double prior_var = tau[j] * tau[j];
        
        double logz_star = findzero_cpp(logzbar_current - 50,
                                        logzbar_current + 50,
                                        .01,
                                        x_all,
                                        y_all, v_all, v_samples2, tau[j],
                                        sigma[j], prior_mean);
        // Rcout << "i = " << i << " - j = " << j << " - zero = " << logz_star << std::endl;
        double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
                                            x_all,
                                            y_all, v_all, v_samples2, tau[j],
                                            sigma[j], prior_mean)));
        
        
        double logz_new = R::rnorm(logz_star, sd_star);
        
        // findzero_cpp(logz_bar[i,j] - 4,
                        //              logz_bar[i,j] + 4,
                        //              tol = .01,
                        //              x_all[1], y_all, v_all, v_delta, tau[j],
                        //                                                  sigma[j], prior_mean)
        
        
        // double posterior_var = 1 / (1 / prior_var + l2 / lik_var);
        // double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
        
        double loglikelihood = logposterior_logz_cpp(v_samples2, sigma[j],
                                                     logzbar_current, logz_new,
                                                     prior_mean, prior_var);
        
        double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
                                                                            x_all,
                                                                            v_all,
                                                                            logzbar_current,
                                                                            logz_new);
        
        double logproposal_ratio = R::dnorm(logzbar_current, logz_star, sd_star, 1) - 
          R::dnorm(logz_new, logz_star, sd_star, 1);
        
        double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
        
        if(R::runif(0,1) < exp(logposterior)){
          
          logz_bar(i, j) = logz_new;
          
        }
        
      }
    }
    sum_m += M_site[i];
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  
  return List::create(_["logz"] = logz2,
                       _["v"] = v2);
  
}