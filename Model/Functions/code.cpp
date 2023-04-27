#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;

const double log2pi = std::log(2.0 * M_PI);

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

double aterm(int n, double x, double t) {
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

double exprnd(double mu) {
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

double randinvg(double mu) {
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}

double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// [[Rcpp::export]]
double rpg(int n, double z){
  
  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }
  
  return(x);
}

///


double jj_m1(double b, double z)  {
  z = fabs(z);
  double m1 = 0.0;
  if (z > 1e-12)
    m1 = b * tanh(z) / z;
  else
    m1 = b * (1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6));
  return m1;
}

double jj_m2(double b, double z) {
  z = fabs(z);
  double m2 = 0.0;
  if (z > 1e-12)
    m2 = (b+1) * b * pow(tanh(z)/z,2) + b * ((tanh(z)-z)/pow(z,3));
  else
    m2 = (b+1) * b * pow(1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6), 2) +
      b * ((-1.0/3) + (2.0/15) * pow(z,2) - (17.0/315) * pow(z,4));
  return m2;
}

double pg_m1(double b, double z) {
  return jj_m1(b, 0.5 * z) * 0.25;
}

double pg_m2(double b, double z) {
  return jj_m2(b, 0.5 * z) * 0.0625;
}

// [[Rcpp::export]]
double rpg_fast(double b, double z){
  
  double m = pg_m1(b, z);
  double v = pg_m2(b,z) - m*m;
  double x = x = R::rnorm(m, sqrt(v));
  
  return(x);
}

// sample from normal distribution

arma::vec mvrnormArma(arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + chol(Sigma) * Y;
}

// compute product between a matrix and a  diagonal matrix (summarised in a vector) 

// [[Rcpp::export]]
arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){
  
  arma::mat result(X.n_rows, D.size());
  for(int i = 0; i < result.n_rows; i++){
    for(int j = 0; j < result.n_cols; j++){
      result(i, j) = X(i,j) * D(j);
    }  
  }
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b, arma::vec& Omega, arma::vec& k){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  arma::vec result = mvrnormArma(invXtOmegaXpB * XtkpBb, invXtOmegaXpB);
  
  // arma::mat L = arma::trans(arma::chol(tXOmegaX + arma::inv(B))); 
  // arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  // arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  // 
  // arma::vec result = mvrnormArma(alpha, arma::trans(arma::inv(arma::trimatl(L))));
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sample_Omega_cpp(arma::mat& X, arma::vec& beta, arma::vec& n){
  
  int nsize = n.size();
  arma::vec Omega_vec(nsize);
  
  for(int i = 0; i < nsize; i++){
    
    arma::vec b = X.row(i) * beta;
    Omega_vec[i] = rpg(n[i], b[0]);
    
  }
  
  return(Omega_vec);
}

// [[Rcpp::export]]
arma::vec sample_betaPG(arma::vec beta, arma::mat X, arma::vec b,
                        arma::mat B, arma::vec n, arma::vec k){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp(X, B, b, Omega, k);
  
  return(beta);
}


// // [[Rcpp::export]]
// double compute_logprob_y_delta0_cpp(arma::vec y_counts, arma::vec c_imk_current, 
//                                     int currentK, double lambda0, double lambdatilde){
//   
//   double sum = 0;
//   for(int k = 0; k < currentK; k++){
//     if(c_imk_current[k] == 2){
//       sum += R::dpois(y_counts[k], lambdatilde, 1);
//     } else {
//       sum += R::dpois(y_counts[k], lambda0, 1);
//     }
//   }
//   
//   return(sum);
// }

// [[Rcpp::export]]
double compute_logprob_y_delta0_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                    int currentK, double n0, double p0, 
                                    double lambda, double lambdatilde){
  
  double sum = 0;
  for(int k = 0; k < currentK; k++){
    if(c_imk_current[k] == 2){
      sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
    } else {
      sum += R::dnbinom(y_counts[k], n0, p0, 1);
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

List update_delta_c_d_cpp_old(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                              arma::cube y, arma::mat v_bar, arma::vec lambda,
                              arma::mat X_z, arma::mat beta_z,
                              arma::mat r, arma::mat alpha,
                              arma::vec M_site, arma::vec K,
                              arma::vec lambdatilde, 
                              double mu0, double n0, 
                              arma::mat u, arma::vec p11, 
                              arma::vec p10, arma::mat theta11, arma::vec theta10,
                              arma::vec emptySites){
  
  int S = y.n_slices;
  int n = M_site.size();
  // arma::cube all_probs = arma::cube(sum(M_site), S, 3 * pow(2, max(K)));
  // arma::cube all_combinations = arma::cube(sum(M_site) * S, 3 * pow(2, max(K)), 2 + max(K));
  
  arma::mat Xzbeta = X_z * beta_z;
  arma::mat ralpha = r * alpha;
  
  int index_m = 0;
  int index_mj = 0;
  
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
        
        double v_im = v_bar(index_m, j) + 
          Xzbeta(i, j) +
          ralpha(index_m, j);
        
        arma::vec log_allProbs;
        arma::mat mat_delta_c_d;
        
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
                                                             n0, p0, 
                                                             lambda[j], 
                                                                   lambdatilde[j]);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
            
            double prob_d = 0;
            for(int l3 = 0; l3 < currentK; l3++){
              prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
            }
            
            log_allProbs[l] = log_prob_y + prob_delta + prob_d;
            
            mat_delta_c_d(l, 0) = 0;
            mat_delta_c_d(l, 1) = 0;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(l, 2 + k) = c_imk_current[k];
            }
            
          }
          
          // delta = 1
          for(int l = 0; l < pow(2, currentK); l++){
            
            arma::vec c_imk_current = DecToBin_cpp(l, currentK);
            
            double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                             c_imk_current, 
                                                             currentK, n0, p0, 
                                                             v_im, 
                                                             lambdatilde[j], 
                                                                        lambda[j],
                                                                              u_im);
            
            double prob_delta = log(theta11(index_m,j));
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
            mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
          // gamma = 1
          for(int l = 0; l < pow(2, currentK); l++){
            
            arma::vec c_imk_current = DecToBin_cpp(l, currentK);
            
            double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                             c_imk_current, 
                                                             currentK, n0, p0, 
                                                             v_im, 
                                                             lambdatilde[j], 
                                                                        lambda[j],
                                                                              u_im);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
            mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
        } 
        
        log_allProbs = log_allProbs - max(log_allProbs);
        arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
        
        // Sample new customer assigment
        NumericVector toSample(allProbs.size());
        for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
        int sampledComb_idx = RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
        
        delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
        gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
        for(int k = 0; k < currentK; k++){
          c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
        }
        
        // // copy output
        // for(int l3 = 0; l3 < allProbs.size(); l3++){
        //   all_probs(index_m, j, l3) = allProbs(l3);
        //   for(int l2 = 0; l2 < (2 + currentK); l2++){
        //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
        //   }
        // }
        
        index_mj += 1;
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                      _["c_imk"] = c_imk,
                      _["gamma"] = gamma);//,
  // _["all_probs"] = all_probs,
  // _["all_combinations"] = all_combinations);
}

// [[Rcpp::export]]
List update_delta_c_d_cpp(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                          arma::cube y, arma::mat v, 
                          arma::vec lambda,
                          arma::vec M_site, arma::vec K,
                          arma::vec lambdatilde, 
                          double mu0, double n0, 
                          arma::mat u, arma::vec p11, 
                          arma::vec p10, arma::mat theta11, arma::vec theta10,
                          arma::vec emptySites){
  
  int S = y.n_slices;
  int n = M_site.size();
  
  int index_m = 0;
  int index_mj = 0;
  
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
        
        double v_im = v(index_m, j);
        
        arma::vec log_allProbs;
        arma::mat mat_delta_c_d;
        
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
                                                             n0, p0, 
                                                             lambda[j], 
                                                                   lambdatilde[j]);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
            
            double prob_d = 0;
            for(int l3 = 0; l3 < currentK; l3++){
              prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
            }
            
            log_allProbs[l] = log_prob_y + prob_delta + prob_d;
            
            mat_delta_c_d(l, 0) = 0;
            mat_delta_c_d(l, 1) = 0;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(l, 2 + k) = c_imk_current[k];
            }
            
          }
          
          // delta = 1
          for(int l = 0; l < pow(2, currentK); l++){
            
            arma::vec c_imk_current = DecToBin_cpp(l, currentK);
            
            double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                             c_imk_current, 
                                                             currentK, n0, p0, 
                                                             v_im, 
                                                             lambdatilde[j], 
                                                                        lambda[j],
                                                                              u_im);
            
            double prob_delta = log(theta11(index_m,j));
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
            mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
          // gamma = 1
          for(int l = 0; l < pow(2, currentK); l++){
            
            arma::vec c_imk_current = DecToBin_cpp(l, currentK);
            
            double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
                                                             c_imk_current, 
                                                             currentK, n0, p0, 
                                                             v_im, 
                                                             lambdatilde[j], 
                                                                        lambda[j],
                                                                              u_im);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
            mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
        } 
        
        log_allProbs = log_allProbs - max(log_allProbs);
        arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
        
        // Sample new customer assigment
        NumericVector toSample(allProbs.size());
        for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
        int sampledComb_idx = RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
        
        delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
        gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
        for(int k = 0; k < currentK; k++){
          c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
        }
        
        // // copy output
        // for(int l3 = 0; l3 < allProbs.size(); l3++){
        //   all_probs(index_m, j, l3) = allProbs(l3);
        //   for(int l2 = 0; l2 < (2 + currentK); l2++){
        //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
        //   }
        // }
        
        index_mj += 1;
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                      _["c_imk"] = c_imk,
                      _["gamma"] = gamma);
}


// [[Rcpp::export]]
double compute_logprob_y_delta1_rnb_cpp(arma::vec y_counts, arma::vec c_imk_current, 
                                        int currentK, double n0, double p0, 
                                        double r_nb,
                                        double v_im, double lambdatilde,
                                        double lambda, arma::vec u_im){
  
  double sum = 0;
  for(int k = 0; k < currentK; k++){
    if(c_imk_current[k] == 1){
      double pi = exp(lambda + v_im + u_im[k]) / (exp(lambda + v_im + u_im[k]) + r_nb);
      sum += R::dnbinom(y_counts[k], r_nb, 1 - pi, 1);
    } else if(c_imk_current[k] == 2){
      sum += R::dpois(y_counts[k], exp(lambda) * lambdatilde, 1);
    } else {
      sum += R::dnbinom(y_counts[k], n0, p0, 1);
    }
  }
  
  return(sum);
}

// [[Rcpp::export]]
List update_delta_c_d_cpp_nb(arma::cube y, arma::mat v, 
                             arma::vec lambda,
                             arma::vec r_nb,
                             arma::vec M_site, arma::vec K,
                             arma::vec lambdatilde, 
                             double mu0, double n0, 
                             arma::mat u, arma::vec p11, 
                             arma::vec p10, arma::mat theta11, arma::vec theta10,
                             arma::vec emptySites){
  
  int S = y.n_slices;
  int n = M_site.size();
  
  int index_m = 0;
  int index_mj = 0;
  
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
        for(int k = 0; k < currentK; k++){
          y_counts[k] = y(index_m, k, j);
          u_im[k] = u(index_m, k);
        }
        
        double v_im = v(index_m, j);
        
        arma::vec log_allProbs;
        arma::mat mat_delta_c_d;
        
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
                                                             n0, p0, 
                                                             lambda[j], 
                                                                   lambdatilde[j]);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
            
            double prob_d = 0;
            for(int l3 = 0; l3 < currentK; l3++){
              prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
            }
            
            log_allProbs[l] = log_prob_y + prob_delta + prob_d;
            
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
                                                                 currentK, n0, p0, 
                                                                 r_nb[j],
                                                                     v_im, 
                                                                     lambdatilde[j], 
                                                                                lambda[j],
                                                                                      u_im);
            
            double prob_delta = log(theta11(index_m,j));
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
            mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
          // gamma = 1
          for(int l = 0; l < pow(2, currentK); l++){
            
            arma::vec c_imk_current = DecToBin_cpp(l, currentK);
            
            double log_prob_y = compute_logprob_y_delta1_rnb_cpp(y_counts, 
                                                                 c_imk_current, 
                                                                 currentK, n0, p0, 
                                                                 r_nb[j],
                                                                     v_im, 
                                                                     lambdatilde[j], 
                                                                                lambda[j],
                                                                                      u_im);
            
            double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
            
            double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
            for(int k = 0; k < currentK; k++){
              prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
            }
            
            mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
            mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
            for(int k = 0; k < currentK; k++){
              mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
            }
            
            log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
            
          }
          
        } 
        
        log_allProbs = log_allProbs - max(log_allProbs);
        arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
        
        // Sample new customer assigment
        NumericVector toSample(allProbs.size());
        for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
        int sampledComb_idx = RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
        
        delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
        gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
        for(int k = 0; k < currentK; k++){
          c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
        }
        
        // // copy output
        // for(int l3 = 0; l3 < allProbs.size(); l3++){
        //   all_probs(index_m, j, l3) = allProbs(l3);
        //   for(int l2 = 0; l2 < (2 + currentK); l2++){
        //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
        //   }
        // }
        
        index_mj += 1;
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                      _["c_imk"] = c_imk,
                      _["gamma"] = gamma);
}

// // [[Rcpp::export]]
// List update_delta_c_d_cpp_single(int i, int m, int j, int index_m, 
//                                  arma::cube y, arma::mat v, 
//                                  arma::vec M_site, arma::vec K,
//                                  arma::vec lambda, arma::vec lambdatilde, 
//                                  double mu0, double n0, 
//                                  arma::mat u, arma::vec p11, 
//                                  arma::vec p10, arma::mat theta11, arma::vec theta10,
//                                  arma::vec emptySites){
//   
//   int S = y.n_slices;
//   int n = M_site.size();
//   // arma::cube all_probs = arma::cube(sum(M_site), S, 3 * pow(2, max(K)));
//   // arma::cube all_combinations = arma::cube(sum(M_site) * S, 3 * pow(2, max(K)), 2 + max(K));
//   
//   // int index_mj = 0;
//   
//   double p0 = n0 / (n0 + mu0);
//   
//   
//   int currentK = K[index_m];
//   
//   arma::vec y_counts = arma::zeros(currentK);
//   arma::vec u_im = arma::zeros(currentK);
//   for(int k = 0; k < currentK; k++){
//     y_counts[k] = y(index_m, k, j);
//     u_im[k] = u(index_m, k);
//   }
//   
//   double v_im = v(index_m, j);
//   
//   arma::vec log_allProbs;
//   arma::mat mat_delta_c_d;
//   
//   if(emptySites[i] == 0){ // not an empty tube
//     
//     log_allProbs = arma::zeros(3 * pow(2, currentK));
//     mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
//                                 2 + currentK);
//     
//     Rcout << "here " << std::endl;
//     
//     // delta = 0, gamma = 0
//     for(int l = 0; l < pow(2, currentK); l++){
//       
//       arma::vec c_imk_current = 2 * DecToBin_cpp(l, currentK);
//       
//       Rcout << "here 1" << std::endl;
//       
//       double log_prob_y = compute_logprob_y_delta0_cpp(y_counts, 
//                                                        c_imk_current, 
//                                                        currentK, 
//                                                        n0, p0, 
//                                                        lambda[j], lambdatilde[j]);
//       
//       double prob_delta = log(1 - theta11(index_m,j)) + log(1 - theta10[j]);
//       
//       double prob_d = 0;
//       for(int l3 = 0; l3 < currentK; l3++){
//         prob_d += R::dbinom(c_imk_current[l3]/ 2.0, 1, p10[j], 1);
//       }
//       
//       log_allProbs[l] = log_prob_y + prob_delta + prob_d;
//       
//       Rcout << "here 2" << std::endl;
//       
//       mat_delta_c_d(l, 0) = 0;
//       mat_delta_c_d(l, 1) = 0;
//       for(int k = 0; k < currentK; k++){
//         mat_delta_c_d(l, 2 + k) = c_imk_current[k];
//       }
//       
//     }
//     
//     // delta = 1
//     for(int l = 0; l < pow(2, currentK); l++){
//       
//       Rcout << "here 3" << std::endl;
//       
//       arma::vec c_imk_current = DecToBin_cpp(l, currentK);
//       
//       double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
//                                                        c_imk_current, 
//                                                        currentK, n0, p0, 
//                                                        lambda[j], lambdatilde[j], 
//                                                                              u_im, v_im);
//       
//       double prob_delta = log(theta11(index_m,j));
//       
//       double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
//       for(int k = 0; k < currentK; k++){
//         prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
//       }
//       
//       Rcout << "here 4" << std::endl;
//       
//       mat_delta_c_d(pow(2, currentK) + l, 0) = 1;
//       mat_delta_c_d(pow(2, currentK) + l, 1) = 0;
//       for(int k = 0; k < currentK; k++){
//         mat_delta_c_d(pow(2, currentK) + l, 2 + k) = c_imk_current[k];
//       }
//       
//       log_allProbs[pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
//       
//     }
//     
//     // gamma = 1
//     for(int l = 0; l < pow(2, currentK); l++){
//       
//       Rcout << "here 5" << std::endl;
//       
//       arma::vec c_imk_current = DecToBin_cpp(l, currentK);
//       
//       double log_prob_y = compute_logprob_y_delta1_cpp(y_counts, 
//                                                        c_imk_current, 
//                                                        currentK, n0, p0, 
//                                                        lambda[j], lambdatilde[j], 
//                                                                              u_im, v_im);
//       
//       double prob_delta = log(1 - theta11(index_m,j)) + log(theta10[j]);
//       
//       
//       Rcout << "here 6" << std::endl;
//       
//       double prob_c = 0;// <- sum(dbinom(c_imk_current, 1, p_11[i], log = T))
//       for(int k = 0; k < currentK; k++){
//         prob_c += R::dbinom(c_imk_current[k], 1, p11[j], 1);
//       }
//       
//       mat_delta_c_d(2 * pow(2, currentK) + l, 0) = 0;
//       mat_delta_c_d(2 * pow(2, currentK) + l, 1) = 1;
//       for(int k = 0; k < currentK; k++){
//         mat_delta_c_d(2 * pow(2, currentK) + l, 2 + k) = c_imk_current[k];
//       }
//       
//       log_allProbs[2 * pow(2, currentK) + l] = log_prob_y + prob_delta + prob_c;
//       
//     }
//     
//   } 
//   
//   log_allProbs = log_allProbs - max(log_allProbs);
//   arma::vec allProbs = exp(log_allProbs) / sum(exp(log_allProbs));
//   
//   // Sample new customer assigment
//   NumericVector toSample(allProbs.size());
//   for(int l5 = 0; l5 < allProbs.size(); l5++) toSample(l5) = l5;
//   int sampledComb_idx = RcppArmadillo::sample(toSample, 1, 1, allProbs)[0];
//   
//   
//   // // copy output
//   // for(int l3 = 0; l3 < allProbs.size(); l3++){
//   //   all_probs(index_m, j, l3) = allProbs(l3);
//   //   for(int l2 = 0; l2 < (2 + currentK); l2++){
//   //     all_combinations(index_mj, l3, l2) = mat_delta_c_d(l3, l2);
//   //   }
//   // }
//   
//   
//   return List::create(_["log_allProbs"] = log_allProbs,
//                       _["mat_delta_c_d"] = mat_delta_c_d);//,
//   // _["all_probs"] = all_probs,
//   // _["all_combinations"] = all_combinations);
// }



// [[Rcpp::export]]
List update_delta_c_d_cpp_existingprobs(arma::cube c_imk, arma::mat delta, arma::mat gamma, 
                                        arma::vec M_site, arma::vec K,
                                        arma::cube &allProbs, arma::cube &allCombs, 
                                        arma::vec emptySites){
  
  int S = c_imk.n_slices;
  int n = M_site.size();
  
  int index_m = 0;
  int index_mj = 0;
  
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        
        int currentK = K[index_m];
        
        arma::vec currentAllProbs = arma::zeros(3 * pow(2, currentK));
        arma::mat mat_delta_c_d = arma::zeros(3 * pow(2, currentK), 
                                              2 + currentK);
        
        for(int l = 0; l < currentAllProbs.size(); l++){
          currentAllProbs[l] = allProbs(index_m, j, l);
          for(int l2 = 0; l2 < (2 + currentK); l2++){
            mat_delta_c_d(l, l2) = allCombs(index_mj, l, l2);
          }
        }
        
        // Sample new customer assigment
        NumericVector toSample(currentAllProbs.size());
        for(int l5 = 0; l5 < currentAllProbs.size(); l5++) toSample(l5) = l5;
        int sampledComb_idx = RcppArmadillo::sample(toSample, 1, 1, currentAllProbs)[0];
        
        delta(index_m, j) = mat_delta_c_d(sampledComb_idx, 0);
        gamma(index_m, j) = mat_delta_c_d(sampledComb_idx, 1);
        for(int k = 0; k < currentK; k++){
          c_imk(index_m, k, j) = mat_delta_c_d(sampledComb_idx, 2 + k);
        }
        
        index_mj += 1;
      }
      
      index_m += 1;
    }
  }
  
  return List::create(_["delta"] = delta,
                      _["c_imk"] = c_imk,
                      _["gamma"] = gamma);
}

arma::vec logisticXb(arma::mat X, arma::vec beta){
  
  arma::vec Xbeta = X * beta;
  
  arma::vec plogistic = 1 / (1 + exp(-Xbeta));
  
  return(plogistic);
}

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

// [[Rcpp::export]]
double loglikbetatheta11_cpp(arma::vec y, arma::mat X, arma::vec beta){
  
  arma::vec plogistic = logisticXb(X, beta);
  
  double loglikelihood = 0;
  for(int l = 0; l < X.n_rows; l++){
    loglikelihood += R::dbinom(y[l], 1, plogistic[l], 1);
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
List update_betatheta11_cpp(arma::mat logz, arma::mat X_z, arma::mat beta_z,
                            arma::mat beta_theta11, 
                            arma::mat theta11, arma::mat delta, 
                            arma::vec r, arma::vec M_site, arma::vec emptySites,
                            arma::vec b_theta11, arma::mat B_theta11){
  
  int S = logz.n_cols;
  int n = M_site.size();
  
  arma::mat Xzbeta = X_z * beta_z;
  
  for(int j = 0; j < S; j++){
    
    arma::mat X_long = arma::zeros(sum(M_site), 3);
    arma::vec y_long = arma::zeros(sum(M_site));
    
    arma::vec beta_theta11_current = arma::conv_to<arma::vec>::from(beta_theta11.row(j));
    
    int l2 = 0;
    int l = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(emptySites[i] == 0){
          y_long[l2] = delta(l, j);
          X_long(l2, 0) = 1;
          X_long(l2, 1) = logz(i, j);
          X_long(l2, 2) = r[l];
          
          l2++;  
        }
        
        l++;
      }
      
    }
    
    arma::vec n2 = arma::ones(l2);
    arma::vec y2 = arma::zeros(l2);
    arma::mat X2 = arma::zeros(l2, 3);
    for(int l3 = 0; l3 < l2; l3++){
      y2[l3] = y_long[l3];
      X2.row(l3) = X_long.row(l3);
    }
    arma::vec k2 = y2 - .5;
    
    arma::vec beta_j = sample_betaPG(beta_theta11_current, X2, b_theta11,
                                     B_theta11, n2, k2);
    beta_theta11.row(j) = arma::conv_to<arma::rowvec>::from(beta_j);
    
    l = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(emptySites[i] == 0){
          
          double Xbeta = beta_theta11(j, 0) +
            beta_theta11(j, 1) * logz(i, j) + 
            beta_theta11(j, 2) * r[l];
          
          theta11(l, j) = 1 / (1 + exp(-Xbeta));
          
        }
        
        
        l++;
      }
      
    }
    
  }
  
  return List::create(_["beta_theta11"] = beta_theta11,
                      _["theta11"] = theta11);
}

// [[Rcpp::export]]
double rinvgamma_cpp(double a, double b){
  return 1 / R::rgamma(a, 1 / b);
}


// [[Rcpp::export]]
arma::mat update_betaz_cpp(arma::mat beta_z, arma::mat l_noempty, arma::vec tau,
                           arma::mat X_beta_noempty, double sigma_beta,
                           arma::vec emptySites, arma::mat B_z){
  
  // X_beta_noempty <- X_beta[emptySites == 0,]
  // l_noempty <- logz[emptySites == 0,]
  int S = tau.size();
  
  int ncov_z = X_beta_noempty.n_cols;
  arma::mat tXX = arma::trans(X_beta_noempty) * X_beta_noempty;
  for(int j = 0; j < S; j++){
    
    arma::mat Lambda_beta = (tXX / (tau[j] * tau[j])) + B_z / (sigma_beta * sigma_beta);
    arma::vec mu_beta = arma::zeros(ncov_z);
    
    for(int i = 0; i < X_beta_noempty.n_rows;i++){
      
      arma::rowvec Xbetal = X_beta_noempty.row(i) * l_noempty(i,j) / (tau[j] * tau[j]);
      arma::vec Xbetal_vec = arma::conv_to<arma::vec>::from(Xbetal);
      mu_beta += Xbetal_vec;
      
    }
    
    // Rcout << mu_beta << std::endl;
    
    // mu_beta <- apply(t(sapply(1:nrow(X_beta_noempty), function(i){
    //   X_beta_noempty[i,] * l_noempty[i,j] / tau[j]^2
    // })), 2, sum)
    
    arma::vec beta_zj = mvrnormArma(arma::inv(Lambda_beta) * mu_beta, arma::inv(Lambda_beta));
    beta_z.col(j) = arma::conv_to<arma::colvec>::from(beta_zj);
    // beta_z[,j] <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
    
  }
  
  return beta_z;
}

// [[Rcpp::export]]
arma::vec update_alpha_cpp(arma::vec alpha, arma::mat v, arma::mat delta, 
                           arma::mat logz, arma::vec r, arma::vec sigma, 
                           double sigma_beta, arma::vec M_site){
  
  int S = alpha.size();
  int n = logz.n_rows;
  
  for(int j = 0; j < S; j++){
    
    arma::vec r_long = arma::zeros(sum(M_site));
    arma::vec y_long = arma::zeros(sum(M_site));
    
    int l = 0;
    int r_idx = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(delta(l, j) == 1){
          y_long[r_idx] = v(l, j) - logz(i,j);
          r_long[r_idx] = r[l];
          
          r_idx++;  
        }
        
        l++;
      }
      
    }
    
    arma::vec r_long2 = arma::zeros(r_idx);
    arma::vec y_long2 = arma::zeros(r_idx);
    for(int l = 0; l < r_idx; l++){
      r_long2[l] = r_long[l];
      y_long2[l] = y_long[l];
    }
    
    arma::vec r_longtrlong = arma::trans(r_long) * r_long;
    double Lambda_beta = (r_longtrlong[0] / pow(sigma[j],2)) + (1.0 / pow(sigma_beta,2));
    double mu_beta = sum(y_long % r_long) / pow(sigma[j],2);
    
    alpha[j] = R::rnorm(mu_beta / Lambda_beta, 1 / Lambda_beta);
    
  }
  
  return(alpha);
}

// [[Rcpp::export]]
arma::vec update_lambda0_NB_cpp(arma::cube &y, arma::cube &c_imk, 
                                double mu0, double n0, 
                                double sd_mu0, double sd_n0,
                                arma::vec M_site, arma::vec K){
  
  int n = M_site.size();
  int S = y.n_slices;
  
  arma::vec nonPCRcounts(sum(M_site) * max(K) * S);
  int index_m = 0;
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[index_m]; k++){
        for(int j = 0; j < S; j++){
          if(c_imk(index_m,k,j) == 0){
            nonPCRcounts[l] = y(index_m,k,j);
            l++;
          }
        }
      }
      index_m++;
    }
  }
  
  NumericVector nonPCRcounts2 = NumericVector(l);
  for(int i = 0; i < l; i++){
    nonPCRcounts2[i] = nonPCRcounts[i];
  }
  
  // nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
  // nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  // propose new sets of parameters
  double mu0_star = R::rnorm(mu0, sd_mu0);
  double n0_star = R::rnorm(n0, sd_n0);
  
  double p0_star = n0_star / (n0_star + mu0_star);
  double p0 = n0 / (n0 + mu0);
  
  if(n0_star > 0){
    
    arma::vec lik_star_all = dnbinom(nonPCRcounts2, n0_star, p0_star, 1);
    double lik_star = sum(lik_star_all);
    arma::vec lik_current_all = dnbinom(nonPCRcounts2, n0, p0, 1);
    double lik_current = sum(lik_current_all);
    
    if(R::runif(0,1) < exp(lik_star - lik_current)){
      mu0 = mu0_star;
      n0 = n0_star;
    }
    
  }
  
  arma::vec toReturn(2);
  toReturn[0] = mu0;
  toReturn[1] = n0;
  
  return(toReturn);
}

// INTERWEAVING 

// [[Rcpp::export]]
List convertSPtoCP_cpp(arma::vec lambda, arma::mat beta_z,
                       arma::vec beta0,
                       arma::vec mu, 
                       arma::mat logz, 
                       arma::mat v, 
                       arma::mat delta, 
                       arma::mat gamma, 
                       arma::mat beta_theta,
                       arma::vec M_site){
  
  int n = M_site.size();
  int S = lambda.size();
  
  arma::vec beta_bar = lambda + beta0;
  arma::vec mu_bar = lambda + mu;
  
  arma::mat logz_bar = arma::zeros(n, S);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < S; j++){
      logz_bar(i,j) = logz(i,j) + beta_bar[j] - beta0[j];
    }  
  }
  
  arma::mat v_bar = arma::zeros(sum(M_site), S);
  
  int l = 0;
  for (int i = 0; i < n; i++) {
    for (int m = 0; m < M_site[i]; m++) {
      for (int j = 0; j < S; j++) {
        if(delta(l, j) == 1){
          v_bar(l, j) = v(l, j)  + logz_bar(i,j) -
            logz(i,j);
        } else if (gamma(l,j) == 1){
          v_bar(l,j) = v(l,j) + mu_bar[j] - mu[j];
        } 
        
      }
      l += 1;
    }
  }
  
  arma::vec beta_theta0_bar = arma::zeros(S);
  for(int j = 0; j < S; j++){
    beta_theta0_bar[j] = beta_theta(j, 0) - lambda[j] * beta_theta(j, 1);
  }
  
  return List::create(_["beta_bar"] = beta_bar,
                      _["beta_theta0_bar"] = beta_theta0_bar,
                      _["mu_bar"] = mu_bar,
                      _["logz_bar"] = logz_bar,
                      _["v_bar"] = v_bar);
  
}

// [[Rcpp::export]]
List convertCPtoSP_cpp(arma::mat beta0_bar,
                       arma::vec lambda, arma::vec mu_bar, 
                       arma::mat logz_bar, arma::mat v_bar, 
                       arma::mat delta, arma::mat gamma, 
                       arma::vec beta_theta0_bar,
                       arma::mat beta_theta,
                       arma::vec M_site){
  
  int n = M_site.size();
  int S = lambda.size();
  
  arma::vec beta0 = beta0_bar - lambda;
  // for(int j = 0; j < S; j++){
  // beta_z(0, j) = beta_bar[j] - lambda[j];
  // }
  arma::vec mu = mu_bar - lambda;
  
  // for(int j = 0; j < S; j++){
  //   arma::vec beta_bar = lambda + arma::conv_to<arma::vec>::from(beta_z.row(0));
  // }
  
  arma::mat logz = arma::zeros(n, S);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < S; j++){
      logz(i,j) = logz_bar(i,j) + beta0[j] - beta0_bar[j];
    }  
  }
  
  arma::mat v = arma::zeros(sum(M_site), S);
  
  int l = 0;
  for (int i = 0; i < n; i++) {
    for (int m = 0; m < M_site[i]; m++) {
      for (int j = 0; j < S; j++) {
        if(delta(l, j) == 1){
          v(l, j) = v_bar(l, j) + logz(i,j) - 
            logz_bar(i,j);
        } else if (gamma(l,j) == 1){
          v(l,j) = v_bar(l,j) + mu[j] - mu_bar[j];
        } 
        
      }
      l += 1;
    }
  }
  
  for(int j = 0; j < S; j++){
    // beta_theta0_bar[j] = beta_theta(j, 0) - lambda[j] * beta_theta(j, 1);
    beta_theta(j, 0) = beta_theta0_bar[j] + lambda[j] * beta_theta(j, 1);
  }
  
  return List::create(_["beta0"] = beta0,
                      _["mu"] = mu,
                      _["logz"] = logz,
                      _["v"] = v,
                      _["beta_theta"] = beta_theta);
  
}

// [[Rcpp::export]]
List convertSPtoNP_cpp(arma::mat logz, arma::mat beta_z, 
                       arma::mat v, arma::vec mu, 
                       arma::mat delta, arma::mat gamma, 
                       arma::vec M_site){
  
  int n = M_site.size();
  int S = mu.size();
  
  arma::mat logz_tilde = arma::zeros(n, S);
  for(int j = 0; j < S; j++){
    for(int i = 0; i < n; i++){
      logz_tilde(i,j) = logz(i,j) - beta_z(0,j); 
    }
  }
  
  arma::mat v_tilde = arma::zeros(sum(M_site), S);
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for (int m = 0; m < M_site[i]; m++) {
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1){
          v_tilde(l, j) = v(l, j) - logz(i,j); 
        } else if (gamma(l, j) == 1){
          v_tilde(l, j) = v(l, j) - mu[j]; 
        } 
      }
      l += 1;
    }
  }
  
  return List::create(_["logz_tilde"] = logz_tilde,
                      _["v_tilde"] = v_tilde);
  
}

// [[Rcpp::export]]
List convertNPtoSP_cpp(arma::mat beta_z, arma::mat logz_tilde, 
                       arma::mat v_tilde, arma::vec mu, 
                       arma::mat delta, arma::mat gamma,
                       arma::vec M_site){
  
  int n = M_site.size();
  int S = mu.size();
  
  arma::mat logz = arma::zeros(n, S);
  for(int j = 0; j < S; j++){
    for(int i = 0; i < n; i++){
      logz(i,j) = logz_tilde(i,j) + beta_z(0,j); 
    }
  }
  
  arma::mat v = arma::zeros(sum(M_site), S);
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for (int m = 0; m < M_site[i]; m++) {
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1){
          v(l, j) = v_tilde(l, j) + logz(i,j); 
        } else if (gamma(l, j) == 1){
          v(l, j) = v_tilde(l, j) + mu[j]; 
        } 
      }
      l += 1;
    }
  }
  
  return List::create(_["logz"] = logz,
                      _["v"] = v);
  
}


double logposterior_lambda_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                               double lambda_current, double lambda_star){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(lambda_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(lambda_current + PCR_v[l]), 1));
    
  }
  
  return(loglikelihood);
  
}

// // [[Rcpp::export]]
// List update_lambda_iw_cpp(arma::mat beta_z, arma::vec mu, arma::vec lambda,
//                           arma::mat logz, arma::mat v, arma::mat u, arma::cube y,
//                           arma::cube c_imk, arma::mat delta, arma::mat gamma,
//                           double sigma_beta, double sigma_mu, arma::vec M_site,
//                           arma::vec K){
//   
//   int n = logz.n_rows;
//   int S = mu.size();
//   
//   // UPDATE LAMBDA CP ------
//   
//   List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site);
//   arma::mat beta_bar = list_CP_cpp["beta_bar"];
//   arma::vec mu_bar = list_CP_cpp["mu_bar"];
//   arma::mat logz_bar = list_CP_cpp["logz_bar"];
//   arma::mat v_bar = list_CP_cpp["v_bar"];
//   
//   // update paramters
//   
//   for(int j = 0; j < S; j++){
//     
//     double sigma_inv = 1 / pow(sigma_mu,2) + 1 / pow(sigma_beta,2);
//     
//     double mu_inv = mu_bar[j] / pow(sigma_mu,2) + beta_bar[j] / pow(sigma_beta,2);
//     
//     double posterior_var = 1 /  sigma_inv;
//     double posterior_mean = mu_inv * posterior_var;
//     
//     lambda[j] = R::rnorm(posterior_mean, sqrt(posterior_var));
//     
//   }
//   
//   List list_SP_cpp = convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
//                                        v_bar, delta, gamma, M_site);
//   arma::mat beta_z2 = list_SP_cpp["beta_z"];
//   arma::vec mu2 = list_SP_cpp["mu"];
//   arma::mat logz2 = list_SP_cpp["logz"];
//   arma::mat v2 = list_SP_cpp["v"];
//   
//   // UPDATE LAMBDA NP ------
//   
//   List list_NP = convertSPtoNP_cpp(logz2, beta_z2, v2, mu2,
//                                    delta, gamma, M_site);
//   arma::mat logz_tilde = list_NP["logz_tilde"];
//   arma::mat v_tilde = list_NP["v_tilde"];
//   
//   // update parameters
//   
//   for(int j = 0; j < S; j++){
//     
//     double sum_v = 0;
//     
//     arma::vec PCR_counts = arma::zeros(sum(M_site) * max(K));
//     arma::vec PCR_v = arma::zeros(sum(M_site) * max(K));
//     
//     int l = 0;
//     int l2 = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         for (int k = 0; k < K[l]; k++) {
//           
//           if(delta(l,j) == 1 &
//              c_imk(l,k,j) == 1){
//             PCR_counts[l2] = y(l,k,j);
//             PCR_v[l2] =  beta_z2(0, j) + logz_tilde(i,j) +
//               v_tilde(l,j) + u(l,k);
//             l2 += 1;
//           } else if (gamma(l,j) == 1 &
//             c_imk(l,k,j) == 1){
//             PCR_counts[l2] = y(l,k,j);
//             PCR_v[l2] = mu2[j] + v_tilde(l,j) + u(l,k);
//             l2 += 1;
//           }
//           
//         }
//         l += 1;
//       }
//     }
//     
//     arma::vec PCR_counts2 = arma::zeros(l2);
//     arma::vec PCR_v2 = arma::zeros(l2);
//     for(int l = 0; l < l2; l++){
//       PCR_counts2[l] = PCR_counts[l];
//       PCR_v2[l] = PCR_v[l];
//     }
//     
//     double mean_likelihood = log(sum(PCR_counts2)) - log(sum(exp(PCR_v2)));
//     double sigma_likelihood = 1.0 / sum(PCR_counts2);
//     
//     double mean_prior = 0;
//     double prior_var = 1000;
//     
//     double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
//     double posterior_mean = ((mean_prior / prior_var) + 
//                              (mean_likelihood / sigma_likelihood)) * posterior_var;
//     double lambda_star = R::rnorm(posterior_mean, sqrt(posterior_var));
//     
//     double lambda_current = lambda[j];
//     
//     double logposterior = logposterior_lambda_cpp(PCR_counts2, PCR_v2,
//                                                   lambda_current, lambda_star);
//     
//     if(R::runif(0,1) < logposterior){
//       lambda[j] = lambda_star;
//     }
//     
//   }
//   
//   list_SP_cpp = convertNPtoSP_cpp(beta_z, logz_tilde, v_tilde,
//                                   mu, delta, gamma, M_site);
//   arma::mat logz3 = list_SP_cpp["logz"];
//   arma::mat v3 = list_SP_cpp["v"];
//   
//   return List::create(_["lambda"] = lambda,
//                       _["mu"] = mu2,
//                       _["beta_z"] = beta_z2,
//                       _["logz"] = logz3,
//                       _["v"] = v3);
//   
// }

// [[Rcpp::export]]
arma::vec sampleLambdaNP(arma::vec lambda, arma::cube c_imk, arma::cube y,
                         arma::mat delta, arma::mat gamma, 
                         arma::mat beta_z, arma::mat logz_tilde, 
                         arma::mat v_tilde, arma::mat u, 
                         arma::vec lambdatilde, arma::vec mu,
                         arma::vec M_site, arma::vec K){
  
  int S = y.n_slices;
  int n = M_site.n_rows;
  
  for (int j = 0; j < S; j++) {
    
    double sum_v = 0;
    double sum_y = 0;
    
    double l = 0;
    for(int i = 0; i < n; i++){
      for(int m = 0; m < M_site[i]; m++){
        for (int k = 0; k < K[l]; k++) {
          
          if(delta(l, j) == 1 &
             c_imk(l, k, j) == 1){
            
            sum_y += y(l, k, j);
            
            sum_v += exp(beta_z(0, j) +
              logz_tilde(i, j) +
              v_tilde(l, j) +
              u(l, k));
            
          } else if (gamma(l, j) == 1 &
            c_imk(l, k, j) == 1){
            sum_y += y(l, k, j);
            sum_v += exp(v_tilde(l, j) +
              mu[j] + u(l, k));
          } else if(c_imk(l, k, j) == 2){
            
            sum_y += y(l, k, j);
            sum_v += lambdatilde[j];
            
          }
          
        }
        
        l++;
      }
    }
    
    lambda[j] = log(R::rgamma(sum_y + 1, 1 / (sum_v + .001)));
    
  }
  
  return(lambda);
}

double logposterior_v_bar_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                              double v_current, double v_star,
                              double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(v_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(v_current + PCR_v[l]), 1));
    
  }
  
  double logprior = R::dnorm(v_star, prior_mean, sqrt(prior_var), 1) - 
    R::dnorm(v_current, prior_mean, sqrt(prior_var), 1);
  
  return(loglikelihood + logprior);
  
  // sum(dpois(PCR_counts, lambda = lambda_j * u_im * wstar, log = T)) -
  //   sum(dpois(PCR_counts, lambda = lambda_j * u_im * wcurrent, log = T)) + 
  //   dnorm(w_star, logz, sigmasq, log = T) - 
  //   dnorm(wcurrent, logz, sigmasq, log = T)
}

// [[Rcpp::export]]
arma::mat update_v_cpp(arma::mat v, arma::mat logz, 
                       arma::vec lambda,  
                       arma::mat X_z, 
                       arma::mat beta_theta,
                       arma::mat u, arma::mat beta_z, 
                       arma::vec beta0,
                       arma::vec mu, arma::cube y,
                       arma::cube c_imk, arma::mat delta, arma::mat gamma,
                       arma::vec sigma, arma::vec sigma_gamma, arma::vec M_site,
                       arma::vec r, arma::vec alpha,
                       arma::vec K, arma::vec emptySites){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  int n = logz_bar.n_rows;
  int S = mu_bar.size();
  
  // update paramters
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      if(!(emptySites[i] == 1)){
        for(int j = 0; j < S; j++){
          if(delta(l, j) == 1 | 
             gamma(l, j) == 1){
            
            arma::vec PCR_counts = arma::zeros(K[l]);
            arma::vec PCR_v = arma::zeros(K[l]);
            
            int l2 = 0;
            for(int k = 0; k < K[l]; k++){
              if(c_imk(l,k,j) == 1){
                
                PCR_counts[l2] = y(l, k, j);
                PCR_v[l2] = u(l,k);
                
                l2 += 1;
              }
            }
            
            arma::vec PCR_counts2 = arma::zeros(l2);
            arma::vec PCR_v2 = arma::zeros(l2);
            for(int k = 0; k < l2; k++){
              PCR_counts2[k] = PCR_counts[k];
              PCR_v2[k] = PCR_v[k];
            }
            
            double prior_mean; 
            double prior_var; 
            
            if(delta(l,j) == 1){
              
              prior_mean = logz_bar(i, j) + alpha[j] * r[l];
              prior_var = sigma[j] * sigma[j];
              
            } else {
              
              prior_mean = mu_bar[j];
              prior_var = sigma_gamma[j] * sigma_gamma[j];
              
            }
            
            double mean_likelihood;
            double sigma_likelihood;
            if(l2 > 0){
              
              mean_likelihood = log(sum(PCR_counts2)) - log(sum(exp(PCR_v2)));
              sigma_likelihood = 1.0 / sum(PCR_counts2);
              
            } else {
              
              mean_likelihood = 0;
              sigma_likelihood = exp(50);
              
            }
            
            double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
            double posterior_mean = ((prior_mean / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var;
            
            double v_current = v_bar(l, j);
            double v_star = R::rnorm(posterior_mean, sqrt(posterior_var));
            
            double logposterior_ratio = logposterior_v_bar_cpp(PCR_counts2, PCR_v2,
                                                               v_current, v_star,
                                                               prior_mean, prior_var);
            
            if(R::runif(0,1) < exp(logposterior_ratio)){
              v_bar(l, j) = v_star;
            }
            
          }
        }
      }
      l += 1;
    }
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  return v2;
  
}

double logposterior_z_iw_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                             double logz_current, double logz_star,
                             double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(logz_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(logz_current + PCR_v[l]), 1));
    
  }
  
  double logprior = R::dnorm(logz_star, prior_mean, sqrt(prior_var), 1) - 
    R::dnorm(logz_current, prior_mean, sqrt(prior_var), 1);
  
  return(loglikelihood + logprior);
  
  // sum(dpois(PCR_counts, lambda = lambda_j * u_im * wstar, log = T)) -
  //   sum(dpois(PCR_counts, lambda = lambda_j * u_im * wcurrent, log = T)) + 
  //   dnorm(w_star, logz, sigmasq, log = T) - 
  //   dnorm(wcurrent, logz, sigmasq, log = T)
}

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
  
  return(loglikelihood + logprior);
  
  // sum(dpois(PCR_counts, lambda = lambda_j * u_im * wstar, log = T)) -
  //   sum(dpois(PCR_counts, lambda = lambda_j * u_im * wcurrent, log = T)) + 
  //   dnorm(w_star, logz, sigmasq, log = T) - 
  //   dnorm(wcurrent, logz, sigmasq, log = T)
}

double logposterior_logz_logistic_cpp(arma::vec y_all, double x_all,
                                      arma::vec v_all,
                                      double logz_current, 
                                      double logz_star){
  
  double logzstar_x = x_all * logz_star;
  double logzcurrent_x = x_all * logz_current;
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < y_all.size(); l++){
    
    double p_current = 1 / (1 + exp(- logzcurrent_x - v_all[l]));
    double p_star = 1 / (1 + exp(- logzstar_x - v_all[l]));
    
    loglikelihood += (R::dbinom(y_all[l], 1, p_star, 1) - 
      R::dbinom(y_all[l], 1, p_current, 1));
    
  }
  
  return(loglikelihood);
  
}

double logf_cpp(double l,
                double x_all,
                double sigmaj,
                arma::vec v_samples,
                arma::vec v,
                arma::vec y,
                double tauj,
                double prior_mean){
  
  double x_all_l = l * x_all;
  
  double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
    sum(- (v_samples - l));
  
  double loglikelihood_v = - sum(x_all * exp(v + x_all_l) / (1 + exp(v + x_all_l))) + 
    sum(y * x_all);
  
  double logprior = - 1 / (pow(tauj,2)) * (l - prior_mean); 
  
  return( loglikelihood + loglikelihood_v + logprior);
}

double findzero_cpp(double a,
                    double b,
                    double tol,
                    double x_all,
                    arma::vec y_all,
                    arma::vec v_all,
                    arma::vec v_samples,
                    double tauj,
                    double sigma_j,
                    double prior_mean){
  
  double c = (a + b) / 2;
  
  double fc = logf_cpp(c,
                       x_all,
                       sigma_j,
                       v_samples,
                       v_all,
                       y_all,
                       tauj,
                       prior_mean);
  
  int nsteps = 0;
  
  while( (b - a) / 2 > tol  & nsteps < 50){
    
    double fa = logf_cpp(a,
                         x_all,
                         sigma_j,
                         v_samples,
                         v_all,
                         y_all,
                         tauj,
                         prior_mean);
    
    if((fc < 0 & fa < 0) | (fc > 0 & fa > 0)){
      a = c;
    } else {
      b = c;
    }  
    
    c = (a + b) / 2;
    
    fc = logf_cpp(c,
                  x_all,
                  sigma_j,
                  v_samples,
                  v_all,
                  y_all,
                  tauj,
                  prior_mean);
    
    nsteps++;
  }
  
  return (a + b) / 2;
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
  
  double x_all_l = l * x_all;
  
  double term1 = - 1 / pow(sigma_j,2) * v_samples.size();// + 
  
  // double term2_asscalar = 
  
  double term2 = -sum(x_all * x_all * exp(v_all + x_all_l) / 
                      ((1 + exp(v_all + x_all_l)) % (1 + exp(v_all + x_all_l)) ) );
  double term3 = - 1 / pow(tauj,2);
  // -sum(x_all * exp(v_all + x_all_l)) + sum(y_all * x_all) -
  // sum(y_all * (v_all + x_all_l)) - sum(log(1 + exp(x_all_l + v_all))) - 
  return(term1 + term2 + term3);
  
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
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
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
        double x_all = beta_theta(j, 1);
        for(int m = 0; m < M_site[i]; m++){
          
          y_all[m] = delta(sum_m + m, j);
          v_all[m] = beta_theta0_bar[j] +
            r[m + sum_m] * beta_theta(j, 2);
          
        }
        
        double prior_mean = beta_bar[j] + Xz_beta(i, j);
        double prior_var = tau[j] * tau[j];
        
        double logz_star = findzero_cpp(logzbar_current - 4,
                                        logzbar_current + 4,
                                        .01,
                                        x_all,
                                        y_all, v_all, v_samples2, tau[j],
                                                                     sigma[j], prior_mean);
        
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
        
        
        double logposterior = loglikelihood + logposterior_ratio_logistic;
        
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
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  
  return List::create(_["logz"] = logz2,
                      _["v"] = v2);
  
}

// [[Rcpp::export]]
arma::mat update_beta0_bar_cpp(arma::vec beta0_bar,
                               arma::mat logz_bar, arma::vec tau,
                               double sigma_beta, arma::vec emptySites){
  
  int n = logz_bar.n_rows;
  int S = beta0_bar.size();
  
  // update parameters
  for(int j = 0; j < S; j++){
    
    
    double sum_logzbar = 0;
    int numSamples = 0;
    
    for(int i = 0; i < n; i++){
      
      if(emptySites[i] != 1){
        
        sum_logzbar += logz_bar(i, j);
        
        numSamples++;
        
      }
      
    }
    
    double lik_mean;
    double lik_var;
    if(numSamples > 0){
      
      lik_mean = sum_logzbar;
      lik_var = tau[j] * tau[j];
      
    } else {
      
      lik_mean = 0;
      lik_var = exp(50);
      
    }
    
    double prior_mean = 0;
    double prior_var = sigma_beta * sigma_beta;
    
    double posterior_var = 1 / (1 / prior_var + numSamples / lik_var);
    double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
    
    beta0_bar[j] = R::rnorm(posterior_mean, sqrt(posterior_var));
    
    
  }
  
  return beta0_bar;
  
}

double logposterior_u_iw_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                             double u_current, double u_star,
                             double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(u_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(u_current + PCR_v[l]), 1));
    
  }
  
  double logprior = R::dnorm(u_star, prior_mean, sqrt(prior_var), 1) - 
    R::dnorm(u_current, prior_mean, sqrt(prior_var), 1);
  
  return(loglikelihood + logprior);
  
}

// [[Rcpp::export]]
List update_u_cpp(arma::mat v, 
                  double zeta, 
                  arma::mat u,
                  arma::vec lambda,
                  arma::vec beta0,
                  arma::mat beta_z, 
                  arma::mat logz,
                  arma::vec mu, 
                  arma::cube y,
                  arma::vec r, 
                  arma::vec alpha,
                  arma::cube c_imk, arma::mat delta, arma::mat gamma,
                  double sigma_u, arma::mat beta_theta,
                  arma::vec sigma, 
                  double sigma_gamma, arma::vec M_site,
                  arma::vec K, arma::vec emptySites){
  
  int n = M_site.size();
  int S = alpha.size();
  
  // UPDATE U CP ------
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  // define parameters in PX
  arma::mat u_hat = u + zeta;
  arma::mat v_hat = v_bar - zeta;
  
  // update zeta
  
  double sum_sigma = 0;
  double sum_mu = 0;
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++)  {
        sum_sigma += (1.0 / pow(sigma_u, 2));
        sum_mu += u_hat(l, k) / pow(sigma_u, 2);
      }
      l += 1;
    }
  }
  
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1){
          sum_mu += (logz_bar(i, j) + r[l] * alpha[j] - v_hat(l, j)) / pow(sigma[j], 2);
          sum_sigma += 1 / pow(sigma[j], 2);
        } else if(gamma(l, j) == 1){
          sum_mu += (mu_bar[j] - v_hat(l, j)) / pow(sigma_gamma, 2);
          sum_sigma += 1 / pow(sigma_gamma, 2);
        } else {
          sum_mu += (- v_hat(l, j));
          sum_sigma += 1;
        }
      }
      l += 1;
    }
  }
  
  double posterior_var = 1.0 / sum_sigma;
  double posterior_mean = sum_mu * posterior_var;
  
  zeta = R::rnorm(posterior_mean, sqrt(posterior_var));
  
  // reupdate parameters
  u = u_hat - zeta;
  v_bar = v_hat + zeta;
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  // UPDATE U CP --------------------------------------------
  
  // update parameters
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        
        arma::vec PCR_counts = arma::zeros(S);
        arma::vec PCR_v = arma::zeros(S);
        
        int l2 = 0;
        for(int j = 0; j < S; j++){
          
          if(c_imk(l, k, j) == 1){
            
            PCR_counts[l2] = y(l, k, j);
            PCR_v[l2] = v_bar(l, j);
            
            l2 += 1;
            
          }
          
        }
        
        arma::vec PCR_counts2 = arma::zeros(l2);
        arma::vec PCR_v2 = arma::zeros(l2);
        for(int j = 0; j < l2; j++){
          PCR_counts2[j] = PCR_counts[j];
          PCR_v2[j] = PCR_v[j];
        }
        
        double mean_prior = 0;
        double prior_var = sigma_u * sigma_u;
        
        double mean_likelihood;
        double sigma_likelihood;
        if(l2 > 0){
          
          mean_likelihood = log(sum(PCR_counts2)) - log(sum(exp(PCR_v2)));
          sigma_likelihood = 1.0 / sum(PCR_counts2);
          
        } else {
          
          mean_likelihood = 0;
          sigma_likelihood = exp(50);
          
        }
        
        double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
        double posterior_mean = ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var;
        
        double u_current = u(l, k);
        double u_star = R::rnorm(posterior_mean, sqrt(posterior_var));
        
        double logposterior_ratio = logposterior_u_iw_cpp(PCR_counts2, PCR_v2,
                                                          u_current, u_star,
                                                          mean_prior, prior_var);
        
        if(R::runif(0,1) < exp(logposterior_ratio)){
          u(l, k) = u_star;
        }
        
      }
      l += 1;
    }
    
  }
  
  return List::create(_["u"] = u,
                      _["zeta"] = zeta,
                      _["v"] = v2);
  
}

double logposterior_mu_iw_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                              double mu_current, double mu_star,
                              double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(mu_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(mu_current + PCR_v[l]), 1));
    
  }
  
  double logprior = R::dnorm(mu_star, prior_mean, sqrt(prior_var), 1) - 
    R::dnorm(mu_current, prior_mean, sqrt(prior_var), 1);
  
  return(loglikelihood + logprior);
  
}

// [[Rcpp::export]]
List update_mu_cpp(arma::vec mu, 
                   arma::vec lambda, 
                   arma::mat delta,
                   arma::mat gamma,
                   arma::vec sigma, 
                   double sigma_gamma, 
                   arma::vec beta0,
                   arma::mat beta_z,
                   arma::mat logz,
                   arma::mat v,
                   arma::mat beta_theta,
                   arma::vec M_site,
                   double sigma_mu, arma::vec emptySites){
  
  int n = M_site.size();
  int S = mu.size();
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  // UPDATE PARAMETERS
  
  for(int j = 0; j < S; j++){
    
    int samples_vbar = 0;
    arma::vec v_bar_mu = arma::zeros(sum(M_site));
    
    int l = 0;
    for(int i = 0; i < n; i++){
      for(int m = 0; m < M_site[i]; m++){
        if(gamma(l, j) == 1){
          v_bar_mu[samples_vbar] = v_bar(l, j);
          samples_vbar++;
        }
        l++;
      }
    }
    
    double lik_mean;
    double lik_var;
    if(samples_vbar > 0){
      lik_mean = sum(v_bar_mu);
      lik_var = sigma_gamma * sigma_gamma; 
    } else {
      lik_mean = 0;
      lik_var = exp(50);
    }
    
    double prior_mean = lambda[j];
    double prior_var = sigma_mu * sigma_mu;
    
    double posterior_var = 1 / (1 / prior_var + samples_vbar / lik_var);
    double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
    
    mu_bar[j] = R::rnorm(posterior_mean, sqrt(posterior_var));
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  return List::create(_["v"] = v2,
                      _["mu"] = mu2);
  
}

double logposterior_betaz_cpp(arma::vec PCR_counts, arma::vec PCR_v,
                              arma::colvec beta_current, arma::colvec beta_star,
                              arma::mat X_z, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
    
    double X_zbeta_star = arma::as_scalar(X_z.row(l) * beta_star);
    double X_zbeta_current = arma::as_scalar(X_z.row(l) * beta_current);
    
    loglikelihood += (R::dpois(PCR_counts[l], exp(X_zbeta_star + PCR_v[l]), 1) - 
      R::dpois(PCR_counts[l], exp(X_zbeta_current + PCR_v[l]), 1));
    
  }
  
  double logprior = 0;
  for(int k = 0; k < beta_current.size(); k++){
    logprior += (R::dnorm(beta_star[k], 0, sqrt(prior_var), 1) - 
      R::dnorm(beta_current[k], 0, sqrt(prior_var), 1));
  }
  
  return(loglikelihood + logprior);
  
}

double logposterior_betaz_logistic_cpp(arma::vec delta_l, arma::mat delta_x,
                                       arma::vec delta_c,
                                       arma::colvec beta_current, arma::colvec beta_star){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < delta_l.size(); l++){
    
    double X_zbeta_star = arma::as_scalar(delta_x.row(l) * beta_star);
    double X_zbeta_current = arma::as_scalar(delta_x.row(l) * beta_current);
    
    double p_current = 1 / (1 + exp(- X_zbeta_current - delta_c[l]));
    double p_star = 1 / (1 + exp(- X_zbeta_star - delta_c[l]));
    
    loglikelihood += (R::dbinom(delta_l[l], 1, p_star, 1) - 
      R::dbinom(delta_l[l], 1, p_current, 1));
    
  }
  
  return(loglikelihood);
  
}

// [[Rcpp::export]]
arma::mat update_beta_cpp(arma::mat v_bar, double zeta, arma::mat u, 
                          arma::mat beta_z, arma::vec mu_bar, arma::cube y,
                          arma::vec r, arma::mat X_z, arma::vec alpha,
                          arma::mat logz_bar,
                          arma::cube c_imk, arma::mat delta, arma::mat gamma,
                          arma::mat Sigma_prop, arma::vec M_site,
                          double sigma_beta,
                          arma::mat beta_theta, 
                          arma::vec K, arma::vec emptySites){
  
  int n = M_site.size();
  int S = alpha.size();
  int ncov = X_z.n_cols;
  
  // update parameters
  for(int j = 0; j < S; j++){
    
    arma::colvec beta_current = beta_z.col(j);
    
    arma::vec PCR_counts = arma::zeros(sum(M_site) * max(K));
    arma::vec PCR_v = arma::zeros(sum(M_site) * max(K));
    arma::mat X_l = arma::zeros(sum(M_site) * max(K), ncov);
    
    int l = 0;
    
    int index_m = 0;
    for(int i = 0; i < n; i++){
      for(int m = 0; m < M_site[i]; m++){
        for(int k = 0; k < K[index_m]; k++){
          
          if(c_imk(index_m, k, j) == 1){
            
            X_l.row(l) = X_z.row(i);
            PCR_counts[l] = y(index_m, k, j);
            PCR_v[l] = v_bar(index_m, j) + 
              alpha[j] * r[index_m] + u(index_m, k);
            
            l += 1;
            
          }
          
        }
        index_m += 1;
      }
      
    }
    
    arma::vec PCR_counts2 = arma::zeros(l);
    arma::vec PCR_v2 = arma::zeros(l);
    arma::mat X_l2 = arma::zeros(l, ncov);
    for(int l2 = 0; l2 < l; l2++){
      PCR_counts2[l2] = PCR_counts[l2];
      PCR_v2[l2] = PCR_v[l2];
      X_l2.row(l2) = X_l.row(l2);
    }
    
    // double mean_prior = 0;
    // double prior_var = sigma_u * sigma_u;
    // 
    // double mean_likelihood;
    // double sigma_likelihood;
    // if(l2 > 0){
    //   
    //   mean_likelihood = log(sum(PCR_counts2)) - log(sum(exp(PCR_v2)));
    //   sigma_likelihood = 1.0 / sum(PCR_counts2);
    //   
    // } else {
    //   
    //   mean_likelihood = 0;
    //   sigma_likelihood = exp(50);
    //   
    // }
    // 
    // double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
    // double posterior_mean = ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var;
    
    arma::vec beta_current_vec = arma::conv_to<arma::vec>::from(beta_current);
    arma::vec beta_star_vec = mvrnormArma(beta_current_vec, Sigma_prop);
    
    arma::colvec beta_star = arma::conv_to<arma::colvec>::from(beta_star_vec);
    
    double logposterior_ratio_y = logposterior_betaz_cpp(PCR_counts2, PCR_v2,
                                                         beta_current, beta_star,
                                                         X_l2, sigma_beta);
    
    arma::vec delta_l = arma::zeros(sum(M_site));
    arma::vec delta_c = arma::zeros(sum(M_site));
    arma::mat delta_x = arma::zeros(sum(M_site), ncov);
    index_m = 0;
    l = 0;
    for(int i = 0; i < n; i++){
      for(int m = 0; m < M_site[i]; m++){
        
        delta_l[l] = delta(index_m, j);
        delta_c[l] = beta_theta(j, 0) + 
          beta_theta(j, 1) * logz_bar(i, j) +
          beta_theta(j, 2) * r[index_m];
        delta_x.row(l) = beta_theta(j, 1) * X_z.row(i);
        
        index_m++;
        l++;
      }
    }
    
    double logposterior_ratio_logistic = logposterior_betaz_logistic_cpp(delta_l,
                                                                         delta_x,
                                                                         delta_c,
                                                                         beta_current, 
                                                                         beta_star);
    
    double logposterior_ratio = logposterior_ratio_y + logposterior_ratio_logistic;
    // Rcout << exp(logposterior_ratio) << std::endl;
    if(R::runif(0,1) < exp(logposterior_ratio)){
      beta_z.col(j) = beta_star;
    }
    
  }
  
  return beta_z;
  
}

// [[Rcpp::export]]
List createMatricesNleq(int j, arma::mat X_z, arma::mat delta, 
                        arma::mat beta_theta, arma::mat logz_bar, arma::vec r,
                        arma::cube c_imk, arma::cube y, arma::vec alpha,
                        arma::mat u, arma::mat v_bar,
                        arma::vec M_site, arma::vec K){
  
  int n = M_site.size();
  int ncov = X_z.n_cols;
  
  arma::vec y_all = arma::zeros(sum(M_site));
  arma::vec v_all = arma::zeros(sum(M_site));
  arma::mat X_all = arma::zeros(sum(M_site), ncov);
  
  // data from delta 
  int index_m = 0;
  for(int i = 0; i < n; i++){
    for (int m = 0; m < M_site[i]; m++) {
      X_all.row(index_m) = beta_theta(j, 1) * X_z.row(i);
      y_all[index_m] = delta(index_m,j);
      v_all[index_m] = beta_theta(j, 0) + 
        beta_theta(j, 1) * logz_bar(i, j) +
        beta_theta(j, 2) * r[index_m];
      index_m++;
    }
  }
  
  // data from y 
  arma::mat X_all2 = arma::zeros(sum(M_site) * max(K), ncov);
  arma::vec y_all2 = arma::zeros(sum(M_site) * max(K));
  arma::vec v_all2 = arma::zeros(sum(M_site) * max(K));
  index_m = 0;
  int l = 0;
  for(int i = 0; i < n; i++){
    for (int m = 0; m < M_site[i]; m++) {
      for (int k = 0; k < K[index_m]; k++) {
        if(c_imk(index_m,k,j) == 1){
          X_all2.row(l) = X_z.row(i);
          y_all2[l] = y(index_m,k,j);
          v_all2[l] = v_bar(index_m,j) +
            alpha[j] * r[index_m] +
            u(index_m,k);
          l++;
        }
      }
      index_m++;
    }
  }
  
  arma::mat X_all3 = arma::zeros(l, ncov);
  arma::vec y_all3 = arma::zeros(l);
  arma::vec v_all3 = arma::zeros(l);
  for(int l2 = 0; l2 < l; l2++){
    X_all3.row(l2) = X_all2.row(l2);
    y_all3[l2] = y_all2[l2];
    v_all3[l2] = v_all2[l2];
  }
  
  return List::create(_["X_all"] = X_all,
                      _["y_all"] = y_all,
                      _["v_all"] = v_all,
                      _["X_all2"] = X_all3,
                      _["y_all2"] = y_all3,
                      _["v_all2"] = v_all3);
}

// [[Rcpp::export]]
arma::vec log_fp_cpp(arma::vec beta, 
                     arma::mat X_all, 
                     arma::vec v_all,
                     arma::vec y_all,
                     arma::mat X_all2, 
                     arma::vec v_all2, 
                     arma::vec y_all2){
  
  int ncov = beta.size();
  
  arma::mat Xbeta = X_all * beta;
  arma::mat Xbeta2 = X_all2 * beta;
  
  arma::vec toReturn = arma::zeros(ncov);
  
  for(int k = 0; k < ncov; k++){
    toReturn[k] = -sum(X_all.col(k) % exp(v_all + Xbeta)) + sum(y_all % X_all.col(k)) +
      (-sum(X_all2.col(k) % exp(v_all2 + Xbeta2)) + sum(y_all2 % X_all2.col(k)));
  }
  
  return(toReturn);
  
}

// [[Rcpp::export]]
arma::mat H_f_cpp(arma::vec beta, 
                  arma::mat X_all, 
                  arma::vec v_all,
                  arma::vec y_all,
                  arma::mat X_all2, 
                  arma::vec v_all2, 
                  arma::vec y_all2){
  
  int ncov = beta.size();
  
  arma::mat Xbeta = X_all * beta;
  arma::mat Xbeta2 = X_all2 * beta;
  
  arma::mat H = arma::zeros(ncov, ncov);
  for(int l1 = 0; l1 < ncov; l1++){
    for(int l2 = 0; l2 <= l1; l2++){
      double H_fl1l2 = -sum(X_all.col(l1) % X_all.col(l2) % exp(v_all + Xbeta) / 
                            ((1 + exp(v_all + Xbeta)) % (1 + exp(v_all + Xbeta))) ) + 
                            (-sum(X_all2.col(l1) % X_all2.col(l2) % exp(v_all2 + Xbeta2)));
      H(l1, l2) = H_fl1l2;
      H(l2, l1) = H(l1, l2);
    }
  }
  
  return(H);
  
}

// [[Rcpp::export]]
double logposterior_beta_cpp(arma::vec beta, 
                             arma::mat X_all, 
                             arma::vec v_all,
                             arma::vec y_all,
                             arma::mat X_all2, 
                             arma::vec v_all2, 
                             arma::vec y_all2){
  
  arma::mat Xbeta = X_all * beta;
  arma::mat Xbeta2 = X_all2 * beta;
  
  double loglikelihood = 0;
  for(int i = 0; i < y_all.size(); i++){
    loglikelihood += R::dbinom(y_all[i], 1, 1 / (1 + exp(-Xbeta[i] - v_all[i])),
                               1);
  }
  
  for(int i = 0; i < y_all2.size(); i++){
    loglikelihood += R::dpois(y_all2[i], exp(Xbeta2[i] + v_all2[i]), 1);
  }
  
  return(loglikelihood);
}

// [[Rcpp::export]]
arma::vec update_sigma_cpp(arma::vec sigma, 
                           arma::vec lambda,
                           arma::mat beta_z,
                           arma::vec beta0,
                           arma::vec mu,
                           arma::mat logz,
                           arma::mat v, 
                           arma::vec r,
                           arma::vec alpha,
                           arma::mat delta,
                           arma::mat gamma,
                           arma::mat beta_theta,
                           double a_sigma, 
                           arma::vec b_sigma, 
                           arma::vec emptySites, 
                           arma::vec M_site){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu,
                                       logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  int S = sigma.size();
  int n = M_site.size();
  
  for(int j = 0; j < S; j++){
    
    int n_samples = 0;
    double sumsq = 0;
    
    int l = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(emptySites[i] == 0){
          if(delta(l,j) == 1){
            sumsq += pow(v_bar(l,j) - (logz_bar(i,j) + r[l] * alpha[j]), 2);
            
            n_samples++;
          }
        }
        
        l++;
      }
      
    }
    
    sigma[j] = sqrt(rinvgamma_cpp(a_sigma + (n_samples / 2), 
                                  b_sigma[j] + (sumsq / 2)));
  }
  
  return sigma;
}


// [[Rcpp::export]]
List update_u_nb_cpp(arma::mat v, 
                     double zeta, 
                     arma::mat u,
                     arma::vec lambda,
                     arma::vec beta0,
                     arma::mat beta_z, 
                     arma::mat logz,
                     arma::vec mu, 
                     arma::cube y,
                     arma::vec r, 
                     arma::vec r_nb,
                     arma::vec alpha,
                     arma::cube c_imk, arma::mat delta, arma::mat gamma,
                     double sigma_u, arma::mat beta_theta,
                     arma::vec sigma, 
                     double sigma_gamma, arma::vec M_site,
                     arma::vec K, arma::vec emptySites){
  
  int n = M_site.size();
  int S = alpha.size();
  
  // UPDATE U CP ------
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  // define parameters in PX
  arma::mat u_hat = u + zeta;
  arma::mat v_hat = v_bar - zeta;
  
  // update zeta
  
  double sum_sigma = 0;
  double sum_mu = 0;
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++)  {
        sum_sigma += (1.0 / pow(sigma_u, 2));
        sum_mu += u_hat(l, k) / pow(sigma_u, 2);
      }
      l += 1;
    }
  }
  
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1){
          sum_mu += (logz_bar(i, j) + r[l] * alpha[j] - v_hat(l, j)) / pow(sigma[j], 2);
          sum_sigma += 1 / pow(sigma[j], 2);
        } else if(gamma(l, j) == 1){
          sum_mu += (mu_bar[j] - v_hat(l, j)) / pow(sigma_gamma, 2);
          sum_sigma += 1 / pow(sigma_gamma, 2);
        } else {
          sum_mu += (- v_hat(l, j));
          sum_sigma += 1;
        }
      }
      l += 1;
    }
  }
  
  double posterior_var = 1.0 / sum_sigma;
  double posterior_mean = sum_mu * posterior_var;
  
  zeta = R::rnorm(posterior_mean, sqrt(posterior_var));
  
  // reupdate parameters
  u = u_hat - zeta;
  v_bar = v_hat + zeta;
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  // UPDATE U CP --------------------------------------------
  
  // update parameters
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        
        arma::vec PCR_counts = arma::zeros(S);
        arma::vec PCR_v = arma::zeros(S);
        arma::vec PCR_rnb = arma::zeros(S);
        double u_current = u(l, k);
        
        int l2 = 0;
        for(int j = 0; j < S; j++){
          
          if(c_imk(l, k, j) == 1){
            
            PCR_counts[l2] = y(l, k, j);
            PCR_v[l2] = v_bar(l, j);
            PCR_rnb[l2] = r_nb[j];
            
            l2 += 1;
            
          }
          
        }
        
        double prior_mean = 0;
        double prior_var = sigma_u * sigma_u;
        
        double mean_likelihood = 0;
        double sum_omegas = 0;
        if(l2 > 0){
          
          arma::vec PCR_counts2 = arma::zeros(l2);
          arma::vec PCR_v2 = arma::zeros(l2);
          arma::vec PCR_rnb2 = arma::zeros(l2);
          for(int j = 0; j < l2; j++){
            PCR_counts2[j] = PCR_counts[j];
            PCR_v2[j] = PCR_v[j];
            PCR_rnb2[j] = PCR_rnb[j];
          }
          
          arma::vec omegas = arma::zeros(l2);
          for(int j = 0; j < l2; j++){
            omegas[j] = rpg_fast(PCR_counts2[j] + PCR_rnb2[j], (PCR_v2[j] + u_current) - 
              log(PCR_rnb2[j]));
          }
          
          arma::vec ytilde = (PCR_counts2 - (PCR_counts2 + PCR_rnb2) / 2.0) / omegas + 
            log(PCR_rnb2) - PCR_v2; 
          double mu_beta = arma::as_scalar(arma::trans(omegas) * ytilde);
          
          mean_likelihood = mu_beta;
          sum_omegas =  sum(omegas);
          
        }
        
        double posterior_var = 1 / ((1 / prior_var) + sum_omegas);
        double posterior_mean = ((prior_mean / prior_var) + mean_likelihood) * posterior_var;
        
        u(l, k) = R::rnorm(posterior_mean, sqrt(posterior_var));
        
      }
      l += 1;
    }
    
  }
  
  return List::create(_["u"] = u,
                      _["zeta"] = zeta,
                      _["v"] = v2);
  
}


// [[Rcpp::export]]
arma::mat update_v_nb_cpp(arma::mat v, arma::mat logz, 
                          arma::vec lambda,  
                          arma::mat X_z, 
                          arma::mat beta_theta,
                          arma::mat u, arma::mat beta_z, 
                          arma::vec beta0,
                          arma::vec r_nb,
                          arma::vec mu, arma::cube y,
                          arma::cube c_imk, arma::mat delta, arma::mat gamma,
                          arma::vec sigma,
                          double sigma_gamma, 
                          arma::vec M_site,
                          arma::vec r, arma::vec alpha,
                          arma::vec K, arma::vec emptySites){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  int n = logz_bar.n_rows;
  int S = mu_bar.size();
  
  // update paramters
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      if(!(emptySites[i] == 1)){
        for(int j = 0; j < S; j++){
          if(delta(l, j) == 1 | 
             gamma(l, j) == 1){
            
            arma::vec PCR_counts = arma::zeros(K[l]);
            arma::vec PCR_u = arma::zeros(K[l]);
            double rnb_current = r_nb[j];
            double v_current = v_bar(l, j);
            
            int l2 = 0;
            for(int k = 0; k < K[l]; k++){
              if(c_imk(l,k,j) == 1){
                
                PCR_counts[l2] = y(l, k, j);
                PCR_u[l2] = u(l,k);
                
                l2 += 1;
              }
            }
            
            double prior_mean; 
            double prior_var; 
            
            if(delta(l,j) == 1){
              
              prior_mean = logz_bar(i, j) + alpha[j] * r[l];
              prior_var = sigma[j] * sigma[j];
              
            } else {
              
              prior_mean = mu_bar[j];
              prior_var = sigma_gamma * sigma_gamma;
              
            }
            
            double mean_likelihood = 0;
            double sum_omegas = 0;
            if(l2 > 0){
              
              arma::vec PCR_counts2 = arma::zeros(l2);
              arma::vec PCR_u2 = arma::zeros(l2);
              for(int k = 0; k < l2; k++){
                PCR_counts2[k] = PCR_counts[k];
                PCR_u2[k] = PCR_u[k];
              }
              
              arma::vec omegas = arma::zeros(l2);
              for(int k = 0; k < l2; k++){
                omegas[k] = rpg_fast(PCR_counts2[k] + rnb_current, (PCR_u2[k] + v_current) - 
                  log(rnb_current));
              }
              
              // double sigma_beta = 1 / sum(omegas);
              
              arma::vec ytilde = (PCR_counts2 - (PCR_counts2 + rnb_current) / 2.0) / omegas + 
                log(rnb_current) - PCR_u2; 
              double mu_beta = arma::as_scalar(arma::trans(omegas) * ytilde);
              
              mean_likelihood = mu_beta;
              sum_omegas =  sum(omegas);
              
            } 
            
            
            double posterior_var = 1 / ((1 / prior_var) + sum_omegas);
            double posterior_mean = ((prior_mean / prior_var) + mean_likelihood) * posterior_var;
            
            v_bar(l, j) = R::rnorm(posterior_mean, sqrt(posterior_var));
            
          }
        }
      }
      l += 1;
    }
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar, 
                                       logz_bar, v_bar, 
                                       delta,
                                       gamma, 
                                       beta_theta0_bar,
                                       beta_theta,
                                       M_site);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  return v2;
  
}

// R SAMPLER 

// [[Rcpp::export]]
double d2loglik_r_cpp(double r, arma::vec y,
                      arma::vec mean_nb){
  
  double sumtrigamma = 0;
  for(int i = 0; i < y.size(); i++){
    sumtrigamma += R::trigamma(y[i] + r);
  }
  
  double out = sumtrigamma - y.size() * R::trigamma(r) + 
    y.size() / r -
    sum(1 / (mean_nb + r)) - 
    sum((mean_nb - y) / pow(mean_nb + y, 2));
  
  return out;
}

// [[Rcpp::export]]
double loglik_star(double r, 
                   arma::vec y,
                   arma::vec mu_nb){
  
  double sum = 0;
  for(int i = 0; i < y.size(); i++){
    double pi = mu_nb[i] / (mu_nb[i] + r);
    sum += R::dnbinom(y[i], r, 1 - pi, 1);
  }
  
  return sum;
} 

// changed function return type and 
// the return type of first parameter
// DO NOT EXPORT THIS FUNCTION VIA RCPP ATTRIBUTES
double obj_fun_rcpp(double& r, 
                    arma::vec& y, arma::vec& mean_nb){
  
  double loglik1 = sum(lgamma(r + y)) - y.size() * R::lgammafn(r);
  
  double loglik2 = mean_nb.size() * r * log(r);
  
  double loglik3 = - arma::as_scalar(arma::trans(y + r) * log(mean_nb + r));
  
  return (-(loglik1 + loglik2 + loglik3));
}


// [[Rcpp::export]]
double optim_rcpp(double& init_val,
                  arma::vec& y, arma::vec& mean_nb){
  
  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optimize = stats["optimize"];
  
  // // Call the optim function from R in C++
  // Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
  //                                // Make sure this function is not exported!
  //                                Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
  //                                Rcpp::_["method"] = "BFGS",
  //                                // Pass in the other parameters as everything
  //                                // is scoped environmentally
  //                                Rcpp::_["y"] = y,
  //                                Rcpp::_["mean_nb"] = mean_nb);
  
  // Call the optim function from R in C++
  Rcpp::List opt_results = optimize(Rcpp::_["f"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                    Rcpp::_["lower"]    = init_val - 1,
                                    Rcpp::_["upper"]    = init_val + 1,
                                    // Pass in the other parameters as everything
                                    // is scoped environmentally
                                    Rcpp::_["y"] = y,
                                    Rcpp::_["mean_nb"] = mean_nb);
  
  // Extract out the estimated parameter values
  double out = opt_results[0];
  
  // Return estimated values
  return out;
}

// [[Rcpp::export]]
arma::vec update_r_nb_cpp(arma::vec r_nb,
                          arma::vec lambda,  
                          arma::mat X_z, 
                          arma::mat beta_theta,
                          arma::mat u, arma::mat beta_z, 
                          arma::vec beta0,
                          arma::vec mu,
                          arma::mat v,
                          arma::mat logz,
                          arma::cube &y,
                          arma::mat delta,
                          arma::mat gamma,
                          arma::cube &c_imk,
                          arma::vec M_site,
                          arma::vec K){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  int n = M_site.size();
  int S = r_nb.size();
  
  for (int j = 0; j < S; j++) {
    
    double sum_v = 0;
    
    arma::vec y_present0 = arma::zeros(sum(M_site) * max(K));
    arma::vec mean_uv0 = arma::zeros(sum(M_site) * max(K));
    
    int l = 0;
    int l2 = 0;
    for(int i = 0; i < n; i++){
      for(int m = 0; m < M_site[i]; m++){
        for (int k = 0; k < K[l]; k++) {
          
          if(c_imk(l,k,j) == 1){
            y_present0[l2] = y(l,k,j);
            mean_uv0[l2] = exp(v_bar(l, j) + u(l, k));
            l2 += 1;
          } 
          
        }
        l += 1;
      }
    }
    
    if(l2 > 0){
      
      arma::vec y_present = arma::zeros(l2);
      arma::vec mean_uv = arma::zeros(l2);
      for(int l = 0; l < l2; l++){
        y_present[l] = y_present0[l];
        mean_uv[l] = mean_uv0[l];
      }
      
      double r_star = optim_rcpp(r_nb[j], y_present, mean_uv);
      double sd_star = - 1 / d2loglik_r_cpp(r_star, y_present, mean_uv);
      
      double r_new = R::rnorm(r_star, sqrt(sd_star));
      
      double loglik_new = loglik_star(r_new, y_present, mean_uv);
      double loglik_current = loglik_star(r_nb[j], y_present, mean_uv);
      
      if(R::runif(0, 1) < exp(loglik_new - loglik_current)){
        r_nb[j] = r_new;
      }
      
    }
    
  }
  
  return(r_nb);
  
}



