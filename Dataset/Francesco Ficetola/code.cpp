#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double dmvnorm_cpp_fast(arma::vec data, arma::vec m, arma::mat rooti, bool returnLog){
  
  int xdim = data.size();
  // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  double constants = -(xdim/2) * log2pi;
  arma::vec z = rooti * ( data - m) ;  
  
  if(returnLog){
    return (constants - 0.5 * arma::sum(z%z) + rootisum);
  } else {
    return exp(constants - 0.5 * arma::sum(z%z) + rootisum);     
  }
  
}

// [[Rcpp::export]]
double dmvnorm_cpp(arma::vec data, arma::vec m, arma::mat Sigma, bool returnLog){
  
  int xdim = data.size();
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  double constants = -(xdim/2) * log2pi;
  arma::vec z = rooti * ( data - m) ;  
  
  if(returnLog){
    return (constants - 0.5 * arma::sum(z%z) + rootisum);
  } else {
    return exp(constants - 0.5 * arma::sum(z%z) + rootisum);     
  }
  
}

// [[Rcpp::export]]
arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + cholsigma * Y;
}


// [[Rcpp::export]]
double logbeta_fun(arma::vec alpha){
  
  return(sum(lgamma(alpha))- lgamma(sum(alpha)));
  
}

// [[Rcpp::export]]
double logsingleterm1(arma::vec gamma_i){
  
  int S = gamma_i.size();
  
  double loglikelihood = 0;
  
  loglikelihood += lgamma(sum(gamma_i));
  // loglikelihood += log(gamma(sum(gamma_i)));
  
  loglikelihood -= sum(lgamma(gamma_i));
  
  return(loglikelihood);
  
}

// [[Rcpp::export]]
double logsingleterm2(int i, int m, arma::vec gamma_i, arma::vec eta_i_jm){
  
  int S = gamma_i.size();
  
  double loglikelihood = 0;
  
  arma::vec gamma_eta = gamma_i + eta_i_jm;
  arma::vec sumgammai = arma::zeros(1);
  sumgammai[0] = sum(gamma_eta);
  
  loglikelihood -= lgamma(sumgammai)[0];
  
  arma::vec log_gamma_eta = lgamma(gamma_eta);
  
  loglikelihood += sum(log_gamma_eta);
  
  return(loglikelihood);
  
}

// [[Rcpp::export]]
double loglikelihood_y(arma::mat beta, arma::mat X, arma::cube& eta_jim, arma::vec M_y){
  
  int Y = eta_jim.n_rows;
  
  arma::mat gamma = exp(X * beta);
  
  double loglikelihood = 0;
  
  for(int i = 0; i < Y; i++){
    
    arma::vec gamma_i = arma::conv_to<arma::vec>::from(gamma.row(i));
    loglikelihood += M_y[i] * logsingleterm1(gamma_i);
    
    for(int m = 0; m < M_y[i]; m++){
      
      arma::vec eta_i_jm = eta_jim.subcube(arma::span(i), arma::span(m), arma::span());
      
      loglikelihood += logsingleterm2(i, m, gamma_i, eta_i_jm);
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double lgammadiff(double x, double a){
  return(a * log(x) + (lgamma(x + a) - (lgamma(x) + a * log(x))));
}

// [[Rcpp::export]]
double loglikelihood_y_optim(arma::mat beta, arma::mat X, arma::cube& eta_jim, arma::vec M_y){
  
  int Y = eta_jim.n_rows;
  int S = beta.n_cols;
  
  arma::mat gamma = exp(X * beta);
  
  double loglikelihood = 0;
  
  for(int i = 0; i < Y; i++){
    
    arma::vec gamma_i = arma::conv_to<arma::vec>::from(gamma.row(i));
    
    for(int m = 0; m < M_y[i]; m++){
      
      arma::vec eta_i_jm = eta_jim.subcube(arma::span(i), arma::span(m), arma::span());
      
      for(int j = 0; j < S; j++){
        
        loglikelihood += (lgamma(gamma_i[j] + eta_i_jm[j]) - lgamma(gamma_i[j]));
        
      }
      
      // int iterToUse = 100000;
      // if(sum(eta_i_jm) > iterToUse){
      //   
      //   for(int l = 1; l <= iterToUse; l++){
      //     loglikelihood -= log(sum(gamma_i) + sum(eta_i_jm) - l);
      //   }
      //   loglikelihood += (lgamma(sum(gamma_i)) - lgamma(sum(gamma_i + eta_i_jm) - iterToUse));
      //   // Rcout << (lgamma(sum(gamma_i)) - lgamma(sum(gamma_i + eta_i_jm) - iterToUse)) << std::endl;
      // } else {
        // loglikelihood += (lgamma(sum(gamma_i)) - lgamma(sum(gamma_i + eta_i_jm)));
      loglikelihood -= lgammadiff(sum(gamma_i), sum(eta_i_jm));
      // }
      
      
    }
  }
    
  return(loglikelihood);
}

// [[Rcpp::export]]
double dlogsingleterm1(int j, arma::vec gamma_i){
  
  int S = gamma_i.size();
  
  double results = 0;
  
  results += R::digamma(sum(gamma_i));
  
  results -= R::digamma(gamma_i[j]);
  
  return(results);
  
}

// [[Rcpp::export]]
double dlogsingleterm2(int j, arma::vec gamma_i, arma::vec eta_i_jm){
  
  int S = gamma_i.size();
  
  double results = 0;
  
  results += R::digamma(gamma_i[j] + eta_i_jm[j]);
  
  results -= R::digamma(sum(gamma_i + eta_i_jm));
  
  return(results);
  
}

// [[Rcpp::export]]
double grad_loglikelihood_y(int p, int j, arma::mat beta, arma::mat X, arma::cube& eta_jim, arma::vec M_y){
  
  int Y = eta_jim.n_rows;
  
  arma::mat gamma = exp(X * beta);
  
  double loglikelihood = 0;
  
  for(int i = 0; i < Y; i++){
    
    arma::vec gamma_i = arma::conv_to<arma::vec>::from(gamma.row(i));
    loglikelihood += M_y[i] * dlogsingleterm1(j, gamma_i) * gamma_i[j] * X(i,p);
    
    for(int m = 0; m < M_y[i]; m++){
      
      arma::vec eta_i_jm = eta_jim.subcube(arma::span(i), arma::span(m), arma::span());
      
      loglikelihood += dlogsingleterm2(j, gamma_i, eta_i_jm) * gamma_i[j] * X(i,p);
      
    }
    
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
double loglikelihood_y_old(arma::mat beta, arma::mat X, arma::cube& y_jim, arma::vec M_y){
  
  int Y = y_jim.n_rows;
  
  arma::mat gamma = exp(X * beta);
  
  double loglikelihood = 0;
  
  for(int i = 0; i < Y; i++){
    
    arma::vec gamma_i = arma::conv_to<arma::vec>::from(gamma.row(i));
    loglikelihood -= M_y[i] * logbeta_fun(gamma_i);
    
    for(int m = 0; m < M_y[i]; m++){
      
      arma::vec y_jim_im = y_jim.subcube(arma::span(i), arma::span(m), arma::span());
      arma::vec gamma_tilde_ik = y_jim_im + gamma_i;
      
      loglikelihood += logbeta_fun(gamma_tilde_ik);
      
    }
    
  }
  
  return(loglikelihood);
}

int delta_kron(int i, int j){
  
  if(i == j){
    return(1);
  } 

  return(0);
  
}

// [[Rcpp::export]]
double hess_logsingleterm1(int j, int j2, arma::vec gamma_i){
  
  double result = gamma_i[j] * ( ( R::trigamma(sum(gamma_i)) - R::trigamma(gamma_i[j]) * delta_kron(j, j2) ) * gamma_i[j2] + 
    (R::digamma(sum(gamma_i)) - R::digamma(gamma_i[j])) * delta_kron(j, j2) );
  
  return(result);
  
}

// [[Rcpp::export]]
double hess_logsingleterm2(int j, int j2, arma::vec gamma_i, arma::vec eta_jim){
  
  double result = gamma_i[j] * ( ( R::trigamma(gamma_i[j] + eta_jim[j]) * delta_kron(j, j2) - R::trigamma(sum(gamma_i + eta_jim))  ) * gamma_i[j2] + 
    (R::digamma(gamma_i[j] + eta_jim[j]) - R::digamma(sum(gamma_i + eta_jim))) * delta_kron(j, j2) );
  
  return(result);
  
}

// [[Rcpp::export]]
double hessian_loglikelihood_y(int p, int j, int p2, int j2, arma::mat beta, arma::mat X, arma::cube& eta_jim, arma::vec M_y){
  
  int Y = eta_jim.n_rows;
  
  arma::mat gamma = exp(X * beta);
  
  double loglikelihood = 0;
  
  for(int i = 0; i < Y; i++){
    
    arma::vec gamma_i = arma::conv_to<arma::vec>::from(gamma.row(i));
    loglikelihood += M_y[i] * X(i,p) * X(i,p2) * hess_logsingleterm1(j, j2, gamma_i);
    
    for(int m = 0; m < M_y[i]; m++){
      
      arma::vec eta_i_jm = eta_jim.subcube(arma::span(i), arma::span(m), arma::span());
      
      loglikelihood += X(i,p) * X(i,p2) * hess_logsingleterm2(j, j2, gamma_i, eta_i_jm);
      
    }
    
  }
  
  return(loglikelihood);
}



