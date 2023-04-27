#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


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
arma::cube updateDeltaGammaC_iter(arma::mat delta,
                                  arma::mat gamma,
                                  arma::cube c_imk,
                                  arma::vec K,
                                  int n,
                                  arma::vec M_site,
                                  int S){
  
  arma::cube deltagammac_output_iter = arma::zeros(sum(M_site), 3 * pow(2, max(K)), S);
  
  int index_m = 0;
  for(int i = 0; i < n; i++){
    
    for(int m = 0; m < M_site[i]; m++){
      
      int currentK = K[index_m];
      
      for(int j = 0; j < S; j++){
        
        arma::vec c_imk_current = arma::zeros(currentK);
        for(int k = 0; k < currentK; k++){
          c_imk_current[k] = c_imk(index_m , k, j);
        }
        int index_deltagammac = convertDeltaIndexes(delta(index_m, j), 
                                                     gamma(index_m, j), 
                                                     c_imk_current, 
                                                     currentK);
        
        deltagammac_output_iter(index_m, index_deltagammac, j) += 1;
        
      }
      
      index_m += 1;
    }
    
  }
  
  return(deltagammac_output_iter);
}


//   // [[Rcpp::export]]
// arma::cube updateA(arma::cube A, arma::){
//   
//   
// }
// pi_hat <- apply(deltagammac_output_iter[,,,1:(iter - 1)], c(1,2,3), mean)
//   
//   for (i in 1:n) {
//     for (m in 1:M_site[i]) {
//       currentK <- K[m + sum(M_site[seq_len(i-1)])]
//       for (j in 1:S) {
//         currentProbs <- pi_hat[m + sum(M_site[seq_len(i-1)]),1:(3 * 2^currentK),j]
//         newProbs <- .95 * currentProbs + .05 * rep(1 / (3 * 2^currentK), 3 * 2^currentK)
//         A[m + sum(M_site[seq_len(i-1)]), j, 1:(3 * 2^currentK)] <- newProbs
//       }
//     }
//   }
  
  