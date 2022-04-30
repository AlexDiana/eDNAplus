#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::depends(RcppParallel)]]
// #include <RcppParallel.h>
// 
// using namespace RcppParallel;

#include <cmath>
#include <algorithm>

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

// [[Rcpp::export]]
double dt2(double x, double mean, double scale, double df){
  double tstat = (x - mean) / scale;
  return(R::dt(tstat, df, 0) / scale);
}

// [[Rcpp::export]]
double rt2(double mean, double scale, double df){
  double tstat = R::rt(df);
  return(tstat * scale + mean);
}

// sample from normal distribution

arma::vec mvrnormArma(arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return mu + arma::trans(Y * arma::chol(Sigma));
}

// [[Rcpp::export]]
arma::vec mrt2(arma::vec mean, arma::mat Sigma, double df){
  arma::vec zeroVec = arma::zeros(mean.size());
  arma::vec y = mvrnormArma(zeroVec, Sigma);
  double u = R::rchisq(df);
  arma::vec x = sqrt(df / u) * y + mean;
  return x;
}

// [[Rcpp::export]]
double dmt_cpp(arma::vec x, double nu, arma::vec mu, arma::mat Sigma, bool returnLog){
  int p = x.size();
  
  double logratio = R::lgammafn(( nu + p ) / 2) - R::lgammafn( nu / 2 );
  
  arma::vec product = (arma::trans(x - mu) * arma::inv(Sigma) * (x - mu));
  double lognum = (- ( nu + p ) / 2) * log(1 + (1 / nu) * product[0]);
  double logden = (p / 2.0) * log(M_PI * nu) + (0.5) * log(arma::det(Sigma));
  
  double loglikelihood = logratio + lognum - logden;
  
  if(returnLog){
    return loglikelihood;
  } else {
    return exp(loglikelihood);
  }
}


//


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

// LAMBERT

const double EPS = 2.2204460492503131e-16;
const double M_1_E = 1.0 / M_E;

double FritschIter(double x, double w_guess){
  double w = w_guess;
  int MaxEval = 5;
  bool CONVERGED = false;
  double k = 2.0 / 3.0;
  int i = 0;
  do {
    double z = std::log(x / w) - w;
    double w1 = w + 1.0;
    double q = 2.0 * w1 * (w1 + k * z);
    double qmz = q - z;
    double e = z / w1 * qmz / (qmz - z);
    CONVERGED = std::abs(e) <= EPS;
    w *= (1.0 + e);
    ++i;
  } while (!CONVERGED && i < MaxEval);
  return(w);
}

// [[Rcpp::export]]
double lambertW0_CS(double x) {
  double result;
  double w;
  if (x == R_PosInf) {
    result = R_PosInf;
  } else if (x < -M_1_E) {
    result = R_NaN;
  } else if (std::abs(x + M_1_E) <= EPS) {
    result = -1.0;
  } else if (x <= M_E - 0.5) {
    if (std::abs(x) <= 1e-16) {
      /* This close to 0 the W_0 branch is best estimated by its Taylor/Pade
       expansion whose first term is the value x and remaining terms are below
       machine double precision. See
       https://math.stackexchange.com/questions/1700919/how-to-derive-the-lambert-w-function-series-expansion
       */
      result = x;
    } else {
      if (std::abs(x) <= 7e-3) {
        /* Use equation (5) in Fritsch */
        w = ((1.33333333333333333 * x + 1.0) * x) /
          ((0.83333333333333333 * x + 2.33333333333333333) * x + 1.0);
      } else {
        /* Use expansion in Corliss 4.22 to create (3, 2) Pade approximant
         Numerator:-10189 / 303840 * p^3 + 40529 / 303840 * p^2 + 489 / 844 * p-1
         Denominator: -14009 / 303840 * p^2 + 355 / 844 * p + 1
         Converted to digits to reduce needed operations
         */
        double p = std::sqrt(2.0 * (M_E * x + 1.0));
        double Numer = ((-0.03353409689310163 * p + 0.1333892838335966) * p +
                        0.5793838862559242) * p - 1.0;
        double Denom = (-0.04610650342285413 * p + 0.4206161137440758) * p + 1.0;
        w = Numer / Denom;
      }
      result = FritschIter(x, w);
    }
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
    w = std::log(x);
    double L_2 = std::log(w);
    double L_3 = L_2 / w;
    double L_3_sq = L_3 * L_3;
    w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
      w + L_3_sq * L_3 / 3.0;
    result = FritschIter(x, w);
  }
  return(result);
}

// [[Rcpp::export]]
double lambertWm1_CS(double x){
  double result;
  double w;
  if (x == 0.0) {
    result = R_NegInf;
  } else if (x < -M_1_E || x > 0.0) {
    result = R_NaN;
  } else if (std::abs(x + M_1_E) <= EPS) {
    result = -1.0;
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
    w = std::log(-x);
    double L_2 = std::log(-w);
    double L_3 = L_2 / w;
    double L_3_sq = L_3 * L_3;
    w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
      w + L_3_sq * L_3 / 3.0;
    result = FritschIter(x, w);
  }
  return(result);
}

// compute product between a matrix and a  diagonal matrix (summarised in a vector) 

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
arma::mat rmtrnorm(arma::mat mu, arma::mat U, arma::mat V) {
  
  arma::mat Xmat(mu.n_rows, mu.n_cols);
  for(int i = 0; i < Xmat.n_rows; i++){
    for(int j = 0; j < Xmat.n_cols; j++){
      Xmat(i, j) = R::rnorm(0, 1);
    }  
  }
  
  arma::mat A = arma::chol(U);
  arma::mat B = arma::chol(V);
  
  arma::mat W = A * Xmat * B;
  
  arma::mat Y = W + mu;
  
  return(Y);
}

/////////////////////////
///////// BETA THETA
/////////////////////////

double logistic(double x){
  return(1 / (1 + exp(-x)));
}

double rtexp(double alpha, 
             double mu_bar){
  return(mu_bar + R::rexp(1 / alpha));
}

double tnorm_std(double mu_bar){
  
  double alpha_star = (mu_bar + sqrt(mu_bar * mu_bar + 4)) / 2;
  
  double z = 0;
  double rho_z = 0;
  do {
    z = rtexp(alpha_star, mu_bar);
    
    rho_z = exp(- (z - alpha_star) * (z - alpha_star) / 2);
    
  } while (!(R::runif(0, 1) < rho_z));
  
  return(z);
  
}

double tnorm(double mu, 
             double sigma, 
             double mu_bar){
  return (mu + sigma * tnorm_std((mu_bar - mu) / sigma));
}

arma::vec truncNormal(arma::vec mu, 
                      arma::mat Sigma, 
                      int truncIndex,
                      double truncation,
                      bool updateBetaTheta1){
  
  int p = mu.size();
  
  // sample marginal of truncated variable
  double x1 = 1;
  if(updateBetaTheta1){
    x1 = tnorm(mu[truncIndex],
               sqrt(Sigma(truncIndex, truncIndex)),
               truncation);
  } 
  
  arma::vec mean1 = arma::zeros(p - 1);
  arma::mat Sigma_11 = arma::zeros(p - 1, p - 1);
  arma::vec Sigma_12 = arma::zeros(p - 1);
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      mean1[k] = mu[k];
      Sigma_12[k] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k, k2 - 1) = Sigma(k, k2);
        }
      }
    }
    if(k > truncIndex){
      mean1[k - 1] = mu[k];
      Sigma_12[k - 1] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k - 1, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k - 1, k2 - 1) = Sigma(k, k2);
        }
      }
    }
  }
  
  // x1 = 2;
  
  arma::vec mu_bar = mean1 + 
    arma::trans(arma::trans(Sigma_12) * (1 / Sigma(truncIndex, truncIndex)) * (x1 - mu[truncIndex]));
  // Rcout << mu_bar << std::endl;
  arma::mat Sigma_bar = Sigma_11 - Sigma_12 * (1 / Sigma(truncIndex, truncIndex)) * arma::trans(Sigma_12);
  
  arma::vec x_left = mvrnormArma(mu_bar, Sigma_bar);
  
  arma::vec out = arma::zeros(p);
  
  out[truncIndex] = x1;
  
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      out[k] = x_left[k];
    }  
    if(k > truncIndex){
      out[k] = x_left[k - 1];
    } 
  }
  
  return out;
  
}


arma::vec sample_beta_cpp_trunc(arma::mat& X, arma::mat& B, 
                                arma::vec& b, arma::vec& Omega, 
                                arma::vec& k, 
                                bool updateBetaTheta1){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  arma::vec result = truncNormal(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 1, 0, updateBetaTheta1);
  // arma::vec result = truncNormal2(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 0, 0);
  
  return(result);
}

arma::vec sample_betaPG_trunc(arma::vec beta, arma::mat X, arma::vec b,
                              arma::mat B, arma::vec n, arma::vec k,
                              bool updateBetaTheta1){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp_trunc(X, B, b, Omega, k, 1);
  // beta = sample_beta_cpp_trunc(X, B, b, Omega, k, updateBetaTheta1);
  
  return(beta);
}

arma::vec truncNormal2(arma::vec mu, 
                       arma::mat Sigma, 
                       int truncIndex,
                       double truncation){
  
  int p = mu.size();
  
  // sample marginal of truncated variable
  double x1 = tnorm(mu[truncIndex],
                    sqrt(Sigma(truncIndex, truncIndex)),
                    truncation);
  
  arma::vec mean1 = arma::zeros(p - 1);
  arma::mat Sigma_11 = arma::zeros(p - 1, p - 1);
  arma::vec Sigma_12 = arma::zeros(p - 1);
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      mean1[k] = mu[k];
      Sigma_12[k] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k, k2 - 1) = Sigma(k, k2);
        }
      }
    }
    if(k > truncIndex){
      mean1[k - 1] = mu[k];
      Sigma_12[k - 1] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k - 1, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k - 1, k2 - 1) = Sigma(k, k2);
        }
      }
    }
  }
  
  // x1 = 2;
  
  arma::vec mu_bar = mean1 + 
    arma::trans(arma::trans(Sigma_12) * (1 / Sigma(truncIndex, truncIndex)) * (x1 - mu[truncIndex]));
  // Rcout << mu_bar << std::endl;
  arma::mat Sigma_bar = Sigma_11 - Sigma_12 * (1 / Sigma(truncIndex, truncIndex)) * arma::trans(Sigma_12);
  
  arma::vec x_left = mvrnormArma(mu_bar, Sigma_bar);
  
  arma::vec out = arma::zeros(p);
  
  out[truncIndex] = x1;
  
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      out[k] = x_left[k];
    }  
    if(k > truncIndex){
      out[k] = x_left[k - 1];
    } 
  }
  
  return out;
  
}

arma::vec sample_beta_cpp_trunc2(arma::mat& X, arma::mat& B, 
                                 arma::vec& b, arma::vec& Omega, 
                                 arma::vec& k){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  // arma::vec result = truncNormal(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 1, 0, updateBetaTheta1);
  arma::vec result = truncNormal2(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 0, 0);
  
  return(result);
}

arma::vec sample_betaPG_trunc2(arma::vec beta, arma::mat X, arma::vec b,
                               arma::mat B, arma::vec n, arma::vec k){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp_trunc2(X, B, b, Omega, k);
  // beta = sample_beta_cpp_trunc(X, B, b, Omega, k, updateBetaTheta1);
  
  return(beta);
}


arma::vec logisticXb(arma::mat X, arma::vec beta){
  
  arma::vec Xbeta = X * beta;
  
  arma::vec plogistic = 1 / (1 + exp(-Xbeta));
  
  return(plogistic);
}


arma::mat createXMatrix(arma::vec M_site,
                        int ncov_w, 
                        int j,
                        arma::mat& logz,
                        arma::mat& X_w,
                        bool updateBetaTheta0){
  
  int n = M_site.size();
  
  arma::mat X = arma::ones(sum(M_site), updateBetaTheta0 + 1 + ncov_w);
  
  int l = 0;
  for(int i = 0; i < n; i++){
    
    for (int m = 0; m < M_site[i]; m++) {
      
      int idx = 0;
      if(updateBetaTheta0){
        X(l, 0) = 1;
        idx = 1;
      }
      
      X(l, idx + 0) = logz(i, j);
      
      for(int cov_w = 0; cov_w < ncov_w; cov_w++){
        X(l, idx + 1 + cov_w) = X_w(l, cov_w);
      }
      
      l++;
    }
    
  }
  
  if(ncov_w > 0){
    // arma::uvec idxes = linspace<uvec>(0,sum(M_site) - 1);
    // arma::uvec idxes2 = linspace<uvec>(0,ncov_w - 1);
    // Rcout << idxes2 << std::endl;
    // Rcout << idxes << std:;endl;
    // X.cols(2, 2 + ncov_w) = X_w;
  }
  
  // arma::mat X = arma::ones(sum(M_site), 2 + ncov_w);
  // 
  // int l = 0;
  // for(int i = 0; i < n; i++){
  //   
  //   for (int m = 0; m < M_site[i]; m++) {
  //     
  //     X(l, 1) = exp(logz(i, j));
  //     
  //     for(int cov_w = 0; cov_w < ncov_w; cov_w++){
  //       X(l, 2 + cov_w) = X_w(l, cov_w);
  //     }
  //     
  //     l++;
  //   }
  //   
  // }
  // 
  // if(ncov_w > 0){
  //   // arma::uvec idxes = linspace<uvec>(0,sum(M_site) - 1);
  //   // arma::uvec idxes2 = linspace<uvec>(0,ncov_w - 1);
  //   // Rcout << idxes2 << std::endl;
  //   // Rcout << idxes << std:;endl;
  //   // X.cols(2, 2 + ncov_w) = X_w;
  // }
  
  return X;
}

// [[Rcpp::export]]
List update_betatheta11_cpp(arma::mat logz, 
                            arma::mat beta_theta11, 
                            arma::mat theta11, arma::mat delta, 
                            arma::mat X_w, arma::vec M_site, 
                            arma::vec b_theta11, arma::mat B_theta11,
                            bool updateBetaTheta0,
                            bool updateBetaTheta1){
  
  int S = logz.n_cols;
  int n = M_site.size();
  int ncov_w = X_w.n_cols;
  
  for(int j = 0; j < S; j++){
    
    // arma::mat X_long = arma::zeros(sum(M_site), 2 + ncov_w);
    // arma::vec y_long = arma::zeros(sum(M_site));
    
    arma::vec beta_theta11_current = arma::zeros(updateBetaTheta0 + 1 + ncov_w);
    int idx = 0;
    if(updateBetaTheta0){
      beta_theta11_current[0] = beta_theta11(j, 0);
      idx = 1;
    }
    
    beta_theta11_current[idx + 0] = beta_theta11(j, 1);
    for(int cov_w = 0; cov_w < ncov_w; cov_w++){
      beta_theta11_current(idx + 1 + cov_w) = beta_theta11(j, 2 + cov_w);
    }
    // arma::vec beta_theta11_current = arma::conv_to<arma::vec>::from(beta_theta11.row(j));
    
    arma::vec n2 = arma::ones(sum(M_site));
    arma::vec y2 = delta.col(j);
    
    arma::mat X2 = createXMatrix(M_site,
                                 ncov_w,
                                 j,
                                 logz,
                                 X_w,
                                 updateBetaTheta0);
    arma::vec k2 = y2 - .5;
    
    arma::vec beta_j = sample_betaPG_trunc2(beta_theta11_current, X2, b_theta11,
                                            B_theta11, n2, k2);
    // beta_theta11.row(j) = arma::conv_to<arma::rowvec>::from(beta_j);
    // int idx = 0;
    
    if(updateBetaTheta0){
      beta_theta11(j, 0) = beta_j[0];
      // idx = 1;
    }
    beta_theta11(j, 1) = beta_j[idx + 0];
    for(int cov_w = 0; cov_w < ncov_w; cov_w++){
      beta_theta11(j, cov_w + 2) = beta_j(idx + cov_w + 1);
    }
    
    int l = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        double Xbeta = beta_theta11(j, 0) + beta_theta11(j, 1) * logz(i, j);
        for(int cov_w = 0; cov_w < ncov_w; cov_w++){
          Xbeta += beta_theta11(j, 2 + cov_w) * X_w(l, cov_w);
        }
        
        theta11(l, j) = logistic(Xbeta);
        
        l++;
      }
      
    }
    
  }
  
  return List::create(_["beta_theta11"] = beta_theta11,
                      _["theta11"] = theta11);
}


//////

// [[Rcpp::export]]
double rinvgamma_cpp(double a, double b){
  return 1 / R::rgamma(a, 1 / b);
}

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
arma::mat update_betaw_cpp(arma::mat beta_w, arma::mat v, arma::mat delta, 
                           arma::mat logz, arma::mat X_w, arma::vec sigma, 
                           double sigma_beta, arma::vec M_site){
  
  int S = beta_w.n_cols;
  int n = logz.n_rows;
  int ncov_w = X_w.n_cols;
  
  arma::mat prior_betaw = (1 / (sigma_beta * sigma_beta)) * arma::eye(ncov_w, ncov_w);
  
  for(int j = 0; j < S; j++){
    
    arma::mat Xw_long = arma::zeros(sum(M_site), ncov_w);
    arma::vec y_long = arma::zeros(sum(M_site));
    
    int l = 0;
    int r_idx = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(delta(l, j) == 1){
          y_long[r_idx] = v(l, j) - logz(i,j);
          Xw_long.row(r_idx) = X_w.row(l);
          
          r_idx++;  
        }
        
        l++;
      }
      
    }
    
    arma::mat Xw_long2 = arma::zeros(r_idx, ncov_w);
    arma::vec y_long2 = arma::zeros(r_idx);
    arma::vec Xy_long2_colsums = arma::zeros(ncov_w);
    for(int l = 0; l < r_idx; l++){
      Xw_long2.row(l) = Xw_long.row(l);
      y_long2[l] = y_long[l];
      Xy_long2_colsums += arma::trans(Xw_long2.row(l)) * y_long[l];
    }
    
    arma::mat r_longtrlong = arma::trans(Xw_long2) * Xw_long2;
    arma::mat Lambda_beta = (r_longtrlong / pow(sigma[j],2)) + prior_betaw;
    // double mu_beta = sum(y_long % r_long) / pow(sigma[j],2);
    arma::vec mu_beta = Xy_long2_colsums / pow(sigma[j],2);
    
    arma::mat Lambdabeta = arma::inv(Lambda_beta);
    arma::vec mean_beta = Lambdabeta * mu_beta;
    
    arma::vec betaw_j = mvrnormArma(mean_beta, Lambdabeta);
    beta_w.col(j) = arma::conv_to<arma::colvec>::from(betaw_j);
    // alpha[j] = R::rnorm(mu_beta / Lambda_beta, 1 / Lambda_beta);
    
  }
  
  return(beta_w);
}

// INTERWEAVING 

// [[Rcpp::export]]
List convertSPtoCP_cpp(arma::vec lambda, 
                       arma::mat beta_z,
                       arma::vec beta0,
                       arma::vec mu, 
                       arma::mat logz, 
                       arma::mat v, 
                       arma::mat delta, 
                       arma::mat gamma, 
                       arma::mat beta_theta,
                       arma::vec M_site,
                       int S_star,
                       int emptyTubes){
  
  int n = M_site.size();
  int S = mu.size();
  
  arma::uvec pos = regspace<uvec>(0, S -1);
  
  arma::vec beta_bar = lambda.elem(pos) + beta0;
  arma::vec mu_bar = lambda.elem(pos) + mu;
  
  arma::mat logz_bar = arma::zeros(n, S);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < S; j++){
      logz_bar(i,j) = logz(i,j) + lambda[j];
    }  
  }
  
  arma::mat v_bar = arma::zeros(sum(M_site) + emptyTubes, S + S_star);
  
  int l = 0;
  for (int l = 0; l < v.n_rows; l++) {
    for(int j = 0; j < (S + S_star); j++){
      v_bar(l, j) = v(l, j) + lambda[j];
    }
  }
  
  arma::vec beta_theta_bar = arma::zeros(S);
  for(int j = 0; j < S; j++){
    // beta_theta_bar(j, 0) = beta_theta(j, 0);
    // beta_theta_bar(j, 1) = beta_theta(j, 1) / exp(lambda[j]);
    beta_theta_bar[j] = beta_theta(j, 0) - beta_theta(j, 1) * lambda[j];
    // beta_theta_bar(j, 1) = beta_theta(j, 1);
    // for(int k = 0; k < (beta_theta.n_cols - 2); k++){
    // beta_theta_bar(j, k + 2) = beta_theta(j, k + 2);
    // }
  }
  
  return List::create(_["beta_bar"] = beta_bar,
                      _["beta_theta_bar"] = beta_theta_bar,
                      _["mu_bar"] = mu_bar,
                      _["logz_bar"] = logz_bar,
                      _["v_bar"] = v_bar);
  
}

// [[Rcpp::export]]
List convertCPtoSP_cpp(arma::mat beta0_bar,
                       arma::vec lambda, arma::vec mu_bar, 
                       arma::mat logz_bar, arma::mat v_bar, 
                       arma::mat delta, arma::mat gamma, 
                       arma::vec beta_theta_bar,
                       arma::vec M_site,
                       int S_star,
                       int emptyTubes){
  
  int n = M_site.size();
  int S = mu_bar.size();
  
  arma::uvec pos = regspace<uvec>(0, S -1);
  
  arma::vec beta0 = beta0_bar - lambda.elem(pos);
  // for(int j = 0; j < S; j++){
  // beta_z(0, j) = beta_bar[j] - lambda[j];
  // }
  arma::vec mu = mu_bar - lambda.elem(pos);
  
  // for(int j = 0; j < S; j++){
  //   arma::vec beta_bar = lambda + arma::conv_to<arma::vec>::from(beta_z.row(0));
  // }
  
  arma::mat logz = arma::zeros(n, S);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < S; j++){
      logz(i,j) = logz_bar(i,j) - lambda[j];
    }  
  }
  
  arma::mat v = arma::zeros(sum(M_site) + emptyTubes, S + S_star);
  
  for (int l = 0; l < v.n_rows; l++) {
    for (int j = 0; j < (S + S_star); j++) {
      v(l, j) = v_bar(l, j) - lambda[j];
    }
    
  }
  
  // arma::mat beta_theta = arma::zeros(S, beta_theta_bar.n_cols);
  // for(int j = 0; j < S; j++){
  //   beta_theta(j, 0) = beta_theta_bar[j] + beta_theta(j, 1) * lambda[j];
  //   // beta_theta(j, 1) = beta_theta_bar(j, 1);
  //   // beta_theta(j, 0) = beta_theta_bar(j, 0);
  //   // beta_theta(j, 1) = beta_theta_bar(j, 1) * exp(lambda[j]);
  //   // for(int k = 0; k < (beta_theta.n_cols - 2); k++){
  //   //   beta_theta(j, k + 2) = beta_theta_bar(j, k + 2);
  //   // }
  // }
  
  return List::create(_["beta0"] = beta0,
                      _["mu"] = mu,
                      _["logz"] = logz,
                      _["v"] = v);
  
}

/////// LAMBDA

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

//////////////////////////////////////////
///////// LOGZ SAMPLER
//////////////////////////////////////////

// [[Rcpp::export]]
double logposterior_logz_cpp(arma::vec v_samples, double sigma,
                             double logz_current, double logz_star){
  // double logz_current, double logz_star,
  // double prior_mean, double prior_var){
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < v_samples.size(); l++){
    
    loglikelihood += (R::dnorm(v_samples[l], logz_star, sigma, 1) - 
      R::dnorm(v_samples[l], logz_current, sigma, 1));
    
  }
  
  // double logprior = R::dnorm(logz_star, prior_mean, sqrt(prior_var), 1) -
  //   R::dnorm(logz_current, prior_mean, sqrt(prior_var), 1);
  
  return(loglikelihood);
  // return(loglikelihood + logprior);
  
}

// [[Rcpp::export]]
double logposterior_logz_logistic_cpp(arma::vec y_all, double x_all,
                                      arma::vec v_all,
                                      double logz_current, 
                                      double logz_star){
  
  double logzstar_x = x_all * logz_star;
  double logzcurrent_x = x_all * logz_current;
  // double logzstar_x = x_all * exp(logz_star);
  // double logzcurrent_x = x_all * exp(logz_current);
  
  double loglikelihood = 0;
  for(int l = 0; (unsigned)l < y_all.size(); l++){
    
    // loglikelihood += (- log(1 + exp(- logzstar_x - v_all[l])) + 
    //   log(1 + exp(- logzcurrent_x - v_all[l]))); 
    // 
    // if(y_all[l] == 0){
    //   
    //   loglikelihood += x_all * (logz_current - logz_star); 
    //   
    // } 
    
    double p_current = 1 / (1 + exp(- logzcurrent_x - v_all[l]));
    double p_star = 1 / (1 + exp(- logzstar_x - v_all[l]));
    
    loglikelihood += (R::dbinom(y_all[l], 1, p_star, 1) -
      R::dbinom(y_all[l], 1, p_current, 1));
    
  }
  
  return(loglikelihood);
  
}

// // [[Rcpp::export]]
// double logf_cpp2(double l,
//                  double x_all,
//                  double sigmaj,
//                  arma::vec v_samples,
//                  arma::vec v,
//                  arma::vec y,
//                  double tauj,
//                  double prior_mean){
//   
//   double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
//     sum(- (v_samples - l));
//   
//   double x_all_l = exp(l) * x_all;
//   
//   double loglikelihood_v = 0;
//   for(int i = 0; i < y.size(); i++){
//     loglikelihood_v += y[i] * x_all_l - 
//       x_all_l * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
//   }
//   
//   double logprior = - (l - prior_mean) / (pow(tauj,2)); 
//   
//   return( loglikelihood + loglikelihood_v + logprior);
// }

// [[Rcpp::export]]
double logf_cpp(double l,
                double x,
                double sigmaj,
                arma::vec v_samples,
                arma::vec v,
                arma::vec y,
                double tauj,
                double prior_mean){
  
  double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
    sum(- (v_samples - l));
  
  double x_all_l = l * x;
  // double x_all_l = exp(l) * x;
  
  // double loglikelihood_v = 0;
  // for(int i = 0; i < y.size(); i++){
  //   loglikelihood_v += y[i] * x_all_l -
  //     x_all_l * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
  // }
  
  double loglikelihood_v = 0;
  for(int i = 0; i < y.size(); i++){
    if(y[i] == 1){
      loglikelihood_v += ( x * exp(-v[i] - x_all_l)) / (1 + exp(-v[i] - x_all_l));
    } else {
      loglikelihood_v += (- x) / (1 + exp(-v[i] - x_all_l));
    }
    // loglikelihood_v += x * (y[i] * exp(l - v[i] - x * exp(l)) +
    //   (y[i] - 1) * exp(l)) / (1 + exp(-v[i] - x * exp(l)));
  }
  
  double logprior = - (l - prior_mean) / (pow(tauj,2)); 
  
  return( loglikelihood + loglikelihood_v + logprior);
}

// [[Rcpp::export]]
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
               double x,
               arma::vec y_all,
               arma::vec v_all,
               arma::vec v_samples,
               double tauj,
               double sigma_j,
               double prior_mean){
  
  double term1 = - 1 / pow(sigma_j,2) * v_samples.size();
  
  // double x_all_l = exp(l) * x_all;
  // 
  // double term2_1 = sum(y_all * x_all_l);
  // 
  // double term2_2 = 0;
  // for(int i = 0; i < v_all.size(); i++){
  //   term2_2 += pow(exp(l + log(x_all) ) / (1 + 1 / exp(x_all_l + v_all[i])), 2);
  // }
  // double term2_3 = 0;
  // for(int i = 0; i < v_all.size(); i++){
  //   term2_3 += x_all * (x_all * exp(l) + 1) * exp(l) /
  //     (1 + 1 / exp(v_all[i] + x_all_l)); 
  // }
  // 
  // double term2 = term2_1 + term2_2 - term2_3;
  
  double x_l = l * x;
  
  double term2 = 0;
  
  for(int i = 0; i < v_all.size(); i++){
    term2 -= x * x * exp(- x_l - v_all[i]) / pow(1 + exp(- x_l - v_all[i]), 2);
  }
  
  // double term2 = term2_1 + term2_2 - term2_3;
  
  double term3 = - 1 / pow(tauj,2);
  
  return(term1 + term2 + term3);
  
}

// [[Rcpp::export]]
arma::mat update_logz_cpp(arma::mat logz, arma::vec beta0,
                          arma::mat X_z, arma::mat beta_z,
                          arma::vec mu, arma::mat v,
                          arma::vec lambda, arma::mat beta_theta,
                          arma::mat X_w,
                          arma::mat beta_w,
                          arma::vec tau, arma::mat delta,
                          arma::mat gamma,
                          arma::vec sigma, 
                          arma::vec M_site,
                          int S_star,
                          int emptyTubes){
  
  double df_t = 3;
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  
  int n = M_site.size();
  int S = beta0.size();
  int ncov_w_theta = X_w.n_cols;
  int ncov_z = beta_z.n_rows;
  
  arma::mat Xz_beta = X_z * beta_z;
  arma::mat Xw_beta = X_w * beta_w;
  arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
  
  if(ncov_w_theta > 0){
    Xw_beta_theta = X_w * arma::trans(beta_theta.submat(0,2,S - 1,2 + (ncov_w_theta - 1)));  
  }
  
  // update parameters
  
  int sum_m = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < S; j++){
      
      double logz_current = logz_bar(i, j);
      
      arma::vec v_samples = arma::zeros(M_site[i]);
      
      int l2 = 0;
      for(int m = 0; m < M_site[i]; m++){
        if(delta(sum_m + m, j) == 1){
          v_samples[l2] = v_bar(m + sum_m, j) - Xw_beta(m + sum_m, j);
          l2 += 1;
        }
      }
      
      arma::vec v_samples2 = arma::zeros(l2);
      // v_samples = v_samples.elem();
      for(int l = 0; l < l2; l++){
        v_samples2[l] = v_samples[l];
      }
      
      arma::vec y_all = arma::zeros(M_site[i]);
      arma::vec v_all = arma::zeros(M_site[i]);
      // double x_all = beta_theta(j, 1) / exp(lambda[j]);
      double x_all = beta_theta(j, 1);
      for(int m = 0; m < M_site[i]; m++){
        
        y_all[m] = delta(sum_m + m, j);
        v_all[m] = beta_theta_bar(j) +
          // v_all[m] = beta_theta(j, 0) +
          Xw_beta_theta(m + sum_m, j);
        
      }
      
      double prior_mean = beta_bar[j] + Xz_beta(i, j);
      double prior_var = tau[j] * tau[j];
      
      double logz_star = findzero_cpp(logz_current - 50,
                                      logz_current + 50,
                                      .01,
                                      x_all,
                                      y_all, v_all, v_samples2, tau[j],
                                                                   sigma[j], prior_mean);
      
      double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
                                          x_all,
                                          y_all, v_all, v_samples2, tau[j],
                                                                       sigma[j], prior_mean)));
      
      
      double logz_new = rt2(logz_star, sd_star, df_t);
      
      double logprior = R::dnorm(logz_new, prior_mean, sqrt(prior_var), 1) -
        R::dnorm(logz_current, prior_mean, sqrt(prior_var), 1);
      
      double loglikelihood_v = logposterior_logz_cpp(v_samples2, sigma[j],
                                                     logz_current, logz_new);
      // double loglikelihood_v = logposterior_logz_cpp(v_samples2, sigma[j],
      //                                                logz_current, logz_new,
      //                                                prior_mean, prior_var);
      
      double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
                                                                          x_all,
                                                                          v_all,
                                                                          logz_current,
                                                                          logz_new);
      
      double logproposal_ratio = log(dt2(logz_current, logz_star, sd_star, df_t)) - 
        log(dt2(logz_new, logz_star, sd_star, df_t));
      
      double logposterior = logprior + loglikelihood_v +
        logposterior_ratio_logistic + logproposal_ratio;
      
      if(R::runif(0,1) < exp(logposterior)){
        
        logz_bar(i, j) = logz_new;
        
      }
      
    }
    
    sum_m += M_site[i];
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar,
                                       logz_bar, v_bar,
                                       delta,
                                       gamma,
                                       beta_theta_bar,
                                       M_site, S_star,
                                       emptyTubes);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  
  return logz2;
  
}

// // [[Rcpp::export]]
// arma::mat update_logz_int_cpp(arma::mat logz, arma::vec beta0,
//                               arma::mat X_z, arma::mat beta_z,
//                               arma::vec mu, arma::mat v,
//                               arma::vec lambda, arma::mat beta_theta,
//                               arma::mat X_w,
//                               arma::mat beta_w,
//                               arma::vec tau, arma::mat delta,
//                               arma::mat gamma,
//                               arma::vec sigma, 
//                               double sigma_beta,
//                               arma::vec M_site,
//                               int S_star){
//   
//   double df_t = 3;
//   
//   int n = M_site.size();
//   int S = beta0.size();
//   int ncov_w_theta = X_w.n_cols;
//   int ncov_z = beta_z.n_rows;
//   
//   arma::mat Xz_beta = X_z * beta_z;
//   arma::mat Xw_beta = X_w * beta_w;
//   arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
//   
//   if(ncov_w_theta > 0){
//     Xw_beta_theta = X_w * beta_theta.submat(0,2,S - 1,2 + ncov_w_theta);  
//   }
//   
//   arma::mat Sigma_beta = arma::zeros(ncov_z, ncov_z);
//   for(int l = 0; l < ncov_z; l++){
//     Sigma_beta(l, l) = sigma_beta * sigma_beta;
//   }
//   
//   // update parameters
//   
//   int sum_m = 0;
//   for(int i = 0; i < n; i++){
//     for(int j = 0; j < S; j++){
//       
//       double logz_current = logz(i, j);
//       
//       arma::vec v_samples = arma::zeros(M_site[i]);
//       
//       int l2 = 0;
//       for(int m = 0; m < M_site[i]; m++){
//         if(delta(sum_m + m, j) == 1){
//           v_samples[l2] = v(m + sum_m, j) - Xw_beta(m + sum_m, j);
//           l2 += 1;
//         }
//       }
//       
//       arma::vec v_samples2 = arma::zeros(l2);
//       // v_samples = v_samples.elem();
//       for(int l = 0; l < l2; l++){
//         v_samples2[l] = v_samples[l];
//       }
//       
//       arma::vec y_all = arma::zeros(M_site[i]);
//       arma::vec v_all = arma::zeros(M_site[i]);
//       double x_all = 1;
//       for(int m = 0; m < M_site[i]; m++){
//         
//         y_all[m] = delta(sum_m + m, j);
//         v_all[m] = beta_theta(j, 0) +
//           Xw_beta_theta(m + sum_m, j);
//         
//       }
//       
//       // double prior_mean = 0;
//       double prior_mean = Xz_beta(i, j);
//       // double prior_var_x = arma::as_scalar(X_z.row(i) * Sigma_beta * arma::trans(X_z.row(i)));
//       // double prior_var = tau[j] * tau[j] + prior_var_x;
//       double prior_var = tau[j] * tau[j];
//       
//       double logz_star = findzero_cpp(logz_current - 50,
//                                       logz_current + 50,
//                                       .01,
//                                       x_all,
//                                       y_all, v_all, v_samples2, tau[j],
//                                                                    sigma[j], prior_mean);
//       
//       double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
//                                           x_all,
//                                           y_all, v_all, v_samples2, tau[j],
//                                                                        sigma[j], prior_mean)));
//       
//       
//       double logz_new = rt2(logz_star, sd_star, df_t);
//       
//       // if(i == 84 & j == 0) Rcout << logz_star << " - " << sd_star << " - " << logz_new << std::endl;
//       double loglikelihood = logposterior_logz_cpp(v_samples2, sigma[j],
//                                                    logz_current, logz_new,
//                                                    prior_mean, prior_var);
//       
//       double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
//                                                                           x_all,
//                                                                           v_all,
//                                                                           logz_current,
//                                                                           logz_new);
//       
//       double logproposal_ratio = log(dt2(logz_current, logz_star, sd_star, df_t)) - 
//         log(dt2(logz_new, logz_star, sd_star, df_t));
//       
//       double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
//       
//       if(R::runif(0,1) < exp(logposterior)){
//         
//         logz(i, j) = logz_new;
//         
//       }
//       
//     }
//     
//     sum_m += M_site[i];
//   }
//   
//   return logz;
//   
// }

// // [[Rcpp::export]]
// List update_logz_cpp(arma::mat logz, arma::vec beta0,
//                      arma::mat X_z, arma::mat beta_z,
//                      arma::vec mu, arma::mat v,
//                      arma::vec lambda, arma::mat beta_theta,
//                      arma::mat X_w,
//                      arma::mat beta_w,
//                      arma::vec tau, arma::mat delta,
//                      arma::mat gamma,
//                      arma::vec sigma, 
//                      arma::vec M_site,
//                      int S_star){
//   
//   double df_t = 3;
//   
//   List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
//                                        gamma, beta_theta, M_site, S_star);
//   arma::vec beta_bar = list_CP_cpp["beta_bar"];
//   arma::vec mu_bar = list_CP_cpp["mu_bar"];
//   arma::mat logz_bar = list_CP_cpp["logz_bar"];
//   arma::mat v_bar = list_CP_cpp["v_bar"];
//   
//   int n = M_site.size();
//   int S = beta_bar.size();
//   int ncov_w_theta = X_w.n_cols;
//   int ncov_z = beta_z.n_rows;
//   
//   arma::mat Xz_beta = X_z * beta_z;
//   arma::mat Xw_beta = X_w * beta_w;
//   arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
//   int sum_m = 0; 
//   for(int i = 0; i < n; i++){
//     for(int m = 0; m < M_site[i]; m++){
//       for(int j = 0; j < S; j++){
//         for(int cov_w = 0; cov_w < ncov_w_theta; cov_w++){
//           Xw_beta_theta(sum_m + m, j) += beta_theta(j, 2 + cov_w) * X_w(sum_m + m, cov_w);
//         }
//       }
//     }
//     sum_m += M_site[i];
//   }
//   
//   // update parameters
//   
//   sum_m = 0;
//   for(int i = 0; i < n; i++){
//     for(int j = 0; j < S; j++){
//       
//       double logzbar_current = logz_bar(i, j);
//       // double logzbar_star = R::rnorm(logz_bar(i, j), sigma_prop);
//       
//       arma::vec v_samples = arma::zeros(M_site[i]);
//       
//       int l2 = 0;
//       for(int m = 0; m < M_site[i]; m++){
//         if(delta(sum_m + m, j) == 1){
//           v_samples[l2] = v_bar(m + sum_m, j) - Xw_beta(m + sum_m, j);//r[m + sum_m] * alpha[j];
//           l2 += 1;
//         }
//       }
//       
//       arma::vec v_samples2 = arma::zeros(l2);
//       for(int l = 0; l < l2; l++){
//         v_samples2[l] = v_samples[l];
//       }
//       
//       arma::vec y_all = arma::zeros(M_site[i]);
//       arma::vec v_all = arma::zeros(M_site[i]);
//       double x_all = beta_theta(j, 1) / exp(lambda[j]);
//       for(int m = 0; m < M_site[i]; m++){
//         
//         y_all[m] = delta(sum_m + m, j);
//         v_all[m] = beta_theta(j, 0) +
//           Xw_beta_theta(m + sum_m, j);
//         // r[m + sum_m] * beta_theta(j, 2);
//         
//       }
//       
//       double prior_mean = beta_bar[j] + Xz_beta(i, j);
//       double prior_var = tau[j] * tau[j];
//       
//       double logz_star = findzero_cpp(logzbar_current - 50,
//                                       logzbar_current + 50,
//                                       .01,
//                                       x_all,
//                                       y_all, v_all, v_samples2, tau[j],
//                                                                    sigma[j], prior_mean);
//       
//       // Rcout << "i = " << i << " - j = " << j << " - zero = " << logz_star << std::endl;
//       double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
//                                           x_all,
//                                           y_all, v_all, v_samples2, tau[j],
//                                                                        sigma[j], prior_mean)));
//       
//       
//       // double logz_new = R::rnorm(logz_star, sd_star);
//       double logz_new = rt2(logz_star, sd_star, df_t);
//       
//       // findzero_cpp(logz_bar[i,j] - 4,
//       //              logz_bar[i,j] + 4,
//       //              tol = .01,
//       //              x_all[1], y_all, v_all, v_delta, tau[j],
//       //                                                  sigma[j], prior_mean)
//       
//       
//       // double posterior_var = 1 / (1 / prior_var + l2 / lik_var);
//       // double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
//       
//       double loglikelihood = logposterior_logz_cpp(v_samples2, sigma[j],
//                                                    logzbar_current, logz_new,
//                                                    prior_mean, prior_var);
//       
//       double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
//                                                                           x_all,
//                                                                           v_all,
//                                                                           logzbar_current,
//                                                                           logz_new);
//       
//       double logproposal_ratio = log(dt2(logzbar_current, logz_star, sd_star , df_t)) - 
//         log(dt2(logz_new, logz_star, sd_star , df_t));
//       // double logproposal_ratio = R::dnorm(logzbar_current, logz_star, sd_star, 1) - 
//       //   R::dnorm(logz_new, logz_star, sd_star, 1);
//       
//       double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
//       
//       if(R::runif(0,1) < exp(logposterior)){
//         
//         logz_bar(i, j) = logz_new;
//         
//       }
//       
//       if(i == 54 & j == 0){
//         Rcout << logzbar_current << " - " << logz_new << " - " << sd_star<< std::endl;  
//         Rcout << logposterior_ratio_logistic << std::endl;  
//         Rcout << logproposal_ratio << std::endl;  
//       }
//       
//     }
//     
//     sum_m += M_site[i];
//   }
//   
//   List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
//                                        lambda, mu_bar, 
//                                        logz_bar, v_bar, 
//                                        delta,
//                                        gamma, 
//                                        beta_theta,
//                                        M_site, S_star);
//   
//   arma::vec beta02 = list_SP_cpp["beta0"];
//   arma::vec mu2 = list_SP_cpp["mu"];
//   arma::mat logz2 = list_SP_cpp["logz"];
//   arma::mat v2 = list_SP_cpp["v"];
//   
//   
//   return List::create(_["logz"] = logz2,
//                       _["v"] = v2);
//   
// }
// 
// // [[Rcpp::export]]
// arma::mat update_logzbar_cpp(arma::mat logz_bar, arma::vec beta_bar,
//                      arma::mat X_z, arma::mat beta_z,
//                      arma::vec mu, arma::mat v_bar,
//                      arma::vec lambda, arma::mat beta_theta,
//                      arma::mat X_w,
//                      arma::mat beta_w,
//                      arma::vec tau, arma::mat delta,
//                      arma::mat gamma,
//                      arma::vec sigma, 
//                      arma::vec M_site,
//                      int S_star){
//   
//   double df_t = 3;
//   
//   int n = M_site.size();
//   int S = beta_bar.size();
//   int ncov_w_theta = X_w.n_cols;
//   int ncov_z = beta_z.n_rows;
//   
//   arma::mat Xz_beta = X_z * beta_z;
//   arma::mat Xw_beta = X_w * beta_w;
//   arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
//   int sum_m = 0; 
//   for(int i = 0; i < n; i++){
//     for(int m = 0; m < M_site[i]; m++){
//       for(int j = 0; j < S; j++){
//         for(int cov_w = 0; cov_w < ncov_w_theta; cov_w++){
//           Xw_beta_theta(sum_m + m, j) += beta_theta(j, 2 + cov_w) * X_w(sum_m + m, cov_w);
//         }
//       }
//     }
//     sum_m += M_site[i];
//   }
//   
//   // update parameters
//   
//   sum_m = 0;
//   for(int i = 0; i < n; i++){
//     for(int j = 0; j < S; j++){
//       
//       double logzbar_current = logz_bar(i, j);
//       
//       arma::vec v_samples = arma::zeros(M_site[i]);
//       
//       int l2 = 0;
//       for(int m = 0; m < M_site[i]; m++){
//         if(delta(sum_m + m, j) == 1){
//           v_samples[l2] = v_bar(m + sum_m, j) - Xw_beta(m + sum_m, j);//r[m + sum_m] * alpha[j];
//           l2 += 1;
//         }
//       }
//       
//       arma::vec v_samples2 = arma::zeros(l2);
//       for(int l = 0; l < l2; l++){
//         v_samples2[l] = v_samples[l];
//       }
//       
//       arma::vec y_all = arma::zeros(M_site[i]);
//       arma::vec v_all = arma::zeros(M_site[i]);
//       double x_all = beta_theta(j, 1) / exp(lambda[j]);
//       for(int m = 0; m < M_site[i]; m++){
//         
//         y_all[m] = delta(sum_m + m, j);
//         v_all[m] = beta_theta(j, 0) +
//           Xw_beta_theta(m + sum_m, j);
//         // r[m + sum_m] * beta_theta(j, 2);
//         
//       }
//       
//       double prior_mean = beta_bar[j] + Xz_beta(i, j);
//       double prior_var = tau[j] * tau[j];
//       
//       double logz_star = findzero_cpp(logzbar_current - 50,
//                                       logzbar_current + 50,
//                                       .01,
//                                       x_all,
//                                       y_all, v_all, v_samples2, tau[j],
//                                                                    sigma[j], prior_mean);
//       // Rcout << "i = " << i << " - j = " << j << " - zero = " << logz_star << std::endl;
//       
//       double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
//                                           x_all,
//                                           y_all, v_all, v_samples2, tau[j],
//                                                                        sigma[j], prior_mean)));
//       
//       
//       // double logz_new = R::rnorm(logz_star, sd_star);
//       double logz_new = rt2(logz_star, sd_star, df_t);
//       
//       double loglikelihood = logposterior_logz_cpp(v_samples2, sigma[j],
//                                                    logzbar_current, logz_new,
//                                                    prior_mean, prior_var);
//       
//       double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
//                                                                           x_all,
//                                                                           v_all,
//                                                                           logzbar_current,
//                                                                           logz_new);
//       
//       double logproposal_ratio = log(dt2(logzbar_current, logz_star, sd_star * sd_star, df_t)) - 
//         log(dt2(logz_new, logz_star, sd_star * sd_star, df_t));
//       
//       double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
//       
//       if(R::runif(0,1) < exp(logposterior)){
//         
//         logz_bar(i, j) = logz_new;
//         
//       }
//       
//     }
//     
//     sum_m += M_site[i];
//   }
//   
//  
//   
//   return logz_bar;
//   
// }


/////////////
/// LOGZ CORRELATED SAMPLER 
///////////////////////

double logposterior_logz_both_cpp(double logz, arma::vec v_samples, 
                                  double sigma, arma::vec y_all, 
                                  double x_all,
                                  arma::vec v_all){
  
  double loglikelihood1 = 0;
  for(int l = 0; (unsigned)l < v_samples.size(); l++){
    
    loglikelihood1 += R::dnorm(v_samples[l], logz, sigma, 1);
    
  }
  
  double logzstar_x = x_all * exp(logz);
  
  double loglikelihood2 = 0;
  for(int l = 0; (unsigned)l < y_all.size(); l++){
    
    double p_star = 1 / (1 + exp(- logzstar_x - v_all[l]));
    
    loglikelihood2 += R::dbinom(y_all[l], 1, p_star, 1);
    
  }
  
  return(loglikelihood1 + loglikelihood2);
  
  // sum(dpois(PCR_counts, lambda = lambda_j * u_im * wstar, log = T)) -
  //   sum(dpois(PCR_counts, lambda = lambda_j * u_im * wcurrent, log = T)) + 
  //   dnorm(w_star, logz, sigmasq, log = T) - 
  //   dnorm(wcurrent, logz, sigmasq, log = T)
}


double logf_likonly_corr_cpp(double l,
                             double x_all,
                             double sigmaj,
                             arma::vec v_samples,
                             arma::vec v,
                             arma::vec y){
  
  double x_all_l = exp(l) * x_all;
  
  double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
    sum(- (v_samples - l));
  
  // double loglikelihood_v = - sum(x_all * exp(v + x_all_l) / (1 + exp(v + x_all_l))) + 
  //   sum(y * x_all);
  double loglikelihood_v = 0;
  for(int i = 0; i < y.size(); i++){
    // loglikelihood_v += y[i] * x_all - 
    // x_all * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
    loglikelihood_v += y[i] * x_all - 
      x_all * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
  }
  
  return( loglikelihood + loglikelihood_v);
}

double findzero_likonly_cpp_corr(double a,
                                 double b,
                                 double tol,
                                 double x_all,
                                 arma::vec y_all,
                                 arma::vec v_all,
                                 arma::vec v_samples,
                                 double sigma_j){
  
  double c = (a + b) / 2;
  
  double fc = logf_likonly_corr_cpp(c,
                                    x_all,
                                    sigma_j,
                                    v_samples,
                                    v_all,
                                    y_all);
  
  int nsteps = 0;
  
  while( (b - a) / 2 > tol  & nsteps < 50){
    
    double fa = logf_likonly_corr_cpp(a,
                                      x_all,
                                      sigma_j,
                                      v_samples,
                                      v_all,
                                      y_all);
    
    if((fc < 0 & fa < 0) | (fc > 0 & fa > 0)){
      a = c;
    } else {
      b = c;
    }  
    
    c = (a + b) / 2;
    
    fc = logf_likonly_corr_cpp(c,
                               x_all,
                               sigma_j,
                               v_samples,
                               v_all,
                               y_all);
    
    nsteps++;
  }
  
  return (a + b) / 2;
}

double h_f_loglik_cpp(double l,
                      double x_all,
                      arma::vec y_all,
                      arma::vec v_all,
                      arma::vec v_samples,
                      double sigma_j){
  
  double x_all_l = exp(l) * x_all;
  
  double term1 = - 1 / pow(sigma_j,2) * v_samples.size();// + 
  
  // double term2_asscalar = 
  
  double term2_1 = sum(y_all * x_all_l);
  
  // double term2_2 = sum(x_all * x_all * exp(2 * (x_all_l + v_all + l)) /
  // ((1 + exp(v_all + x_all_l)) % (1 + exp(v_all + x_all_l))) );
  // double term2_2 = sum((exp(x_all_l + v_all + l + log(x_all)) / (1 + exp(v_all + x_all_l))) %
  //                      (exp(x_all_l + v_all + l + log(x_all)) / (1 + exp(v_all + x_all_l))) );
  // double term2_3 = sum(x_all * (x_all_l + 1) * exp(x_all_l + v_all + l) / (1 + exp(v_all + x_all_l)));
  
  double term2_2 = 0;
  for(int i = 0; i < v_all.size(); i++){
    term2_2 += pow(exp(l + log(x_all) ) / (1 + 1 / exp(x_all_l + v_all[i])), 2);
    // term2_2 += exp( 2 * ((x_all_l + v_all[i] + l + log(x_all)) - 
    //   log(1 + exp(v_all[i] + x_all_l))));
  }
  double term2_3 = 0;
  for(int i = 0; i < v_all.size(); i++){
    term2_3 += x_all * exp(l) / 
      (1 + 1 / exp(v_all[i] + x_all_l));
    // term2_3 += exp( log(x_all) + log(1 + x_all * exp(l)) + (x_all * exp(l) + v_all[i] + l)  - log(1 + exp(v_all[i] + x_all_l)));
  }
  
  double term2 = term2_1 + term2_2 - term2_3;
  
  // double term2 = -sum(x_all * x_all * exp(v_all + x_all_l) / 
  //                     ((1 + exp(v_all + x_all_l)) % (1 + exp(v_all + x_all_l)) ) );
  
  // double term3 = - 1 / pow(tauj,2);
  // -sum(x_all * exp(v_all + x_all_l)) + sum(y_all * x_all) -
  // sum(y_all * (v_all + x_all_l)) - sum(log(1 + exp(x_all_l + v_all))) - 
  
  Rcout << "y_all = " << y_all << 
    " - v_all = " << v_all << 
      " - x_all = " << x_all << " - logz_star_mean = " << 
        l << " - v_samples = " << v_samples <<
          " - term2_1 = " << term2_1 << 
            " - term2_2 = " << term2_2 << 
              " - term2_3 = " << term2_3 << 
                " - term2 = " << term2 << 
                  " - term1 = " << term1 << std::endl;
  
  return(term1 + term2);
  
}

// // [[Rcpp::export]]
// List update_logz_corr_cpp(arma::mat logz, arma::vec beta0,
//                           arma::mat X_z, arma::mat beta_z,
//                           arma::vec mu, arma::mat v,
//                           arma::vec lambda, arma::mat beta_theta,
//                           arma::mat X_w,
//                           arma::mat beta_w,
//                           arma::mat Tau, 
//                           arma::mat delta,
//                           arma::mat gamma,
//                           arma::vec sigma, 
//                           arma::vec M_site,
//                           int S_star){
//   
//   double df_t = 3;
//   
//   List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
//                                        gamma, beta_theta, M_site, S_star);
//   arma::vec beta_bar = list_CP_cpp["beta_bar"];
//   arma::vec mu_bar = list_CP_cpp["mu_bar"];
//   arma::mat logz_bar = list_CP_cpp["logz_bar"];
//   arma::mat v_bar = list_CP_cpp["v_bar"];
//   
//   int n = M_site.size();
//   int S = beta_bar.size();
//   int ncov_z = beta_z.n_rows;
//   int ncov_w_theta = X_w.n_cols;
//   arma::mat Xz_beta = X_z * beta_z;
//   arma::mat Xw_beta = X_w * beta_w;
//   arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
//   int sum_m = 0; 
//   for(int i = 0; i < n; i++){
//     for(int m = 0; m < M_site[i]; m++){
//       for(int j = 0; j < S; j++){
//         for(int cov_w = 0; cov_w < ncov_w_theta; cov_w++){
//           Xw_beta_theta(sum_m + m, j) += beta_theta(j, 2 + cov_w) * X_w(sum_m + m, cov_w);
//         }
//       }
//     }
//     sum_m += M_site[i];
//   }
//   // arma::mat invTau = arma::inv(Tau);
//   
//   // update parameters
//   
//   sum_m = 0;
//   for(int i = 0; i < n; i++){
//     
//     for(int j = 0; j < S; j++){
//       
//       double logzbar_current = logz_bar(i, j);
//       // double logzbar_star = R::rnorm(logz_bar(i, j), sigma_prop);
//       
//       arma::vec v_samples = arma::zeros(M_site[i]);
//       
//       int l2 = 0;
//       for(int m = 0; m < M_site[i]; m++){
//         if(delta(sum_m + m, j) == 1){
//           v_samples[l2] = v_bar(m + sum_m, j) - Xw_beta(m + sum_m, j);//r[m + sum_m] * alpha[j];
//           l2 += 1;
//         }
//       }
//       
//       arma::vec v_samples2 = arma::zeros(l2);
//       for(int l = 0; l < l2; l++){
//         v_samples2[l] = v_samples[l];
//       }
//       
//       arma::vec y_all = arma::zeros(M_site[i]);
//       arma::vec v_all = arma::zeros(M_site[i]);
//       double x_all = beta_theta(j, 1) / exp(lambda[j]);
//       for(int m = 0; m < M_site[i]; m++){
//         
//         y_all[m] = delta(sum_m + m, j);
//         v_all[m] = beta_theta(j, 0) +
//           Xw_beta_theta(m + sum_m, j);
//         // r[m + sum_m] * beta_theta(j, 2);
//         
//       }
//       
//       arma::vec Sigma_12 = arma::zeros(S - 1);
//       arma::mat Sigma_22 = arma::zeros(S - 1, S - 1);
//       arma::vec aminusmu = arma::zeros(S - 1);
//       
//       for(int l = 0; l < (S - 1); l++){
//         if(l < j){
//           Sigma_12[l] = Tau(j, l); 
//           aminusmu[l] = logz_bar(i, l) - (beta_bar[l] + Xz_beta(i, l));
//         } else {
//           Sigma_12[l] = Tau(j, l + 1);   
//           aminusmu[l] = logz_bar(i, l + 1) - (beta_bar[l + 1] + Xz_beta(i, l + 1));
//         }
//         for(int l2 = 0; l2 < (S - 1); l2++){
//           if(l < j){
//             if(l2 < j){
//               Sigma_22(l, l2) = Tau(l, l2); 
//             } else {
//               Sigma_22(l, l2) = Tau(l, l2 + 1);   
//             }
//           } else {
//             if(l2 < j){
//               Sigma_22(l, l2) = Tau(l + 1, l2); 
//             } else {
//               Sigma_22(l, l2) = Tau(l + 1, l2 + 1);   
//             }
//           }
//         }
//       }
//       
//       arma::mat invSigma_22 = arma::inv(Sigma_22);
//       
//       double prior_mean_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * aminusmu);
//       
//       double prior_var_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * Sigma_12);
//       
//       double prior_mean = beta_bar[j] + Xz_beta(i, j) + prior_mean_2;
//       
//       double prior_sd = sqrt(Tau(j, j) - prior_var_2);
//       
//       double logz_star = findzero_cpp(logzbar_current - 50,
//                                       logzbar_current + 50,
//                                       .01,
//                                       x_all,
//                                       y_all, v_all, v_samples2, prior_sd,
//                                       sigma[j], prior_mean);
//       // Rcout << "i = " << i << " - j = " << j << " - zero = " << logz_star << std::endl;
//       double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
//                                           x_all,
//                                           y_all, v_all, v_samples2, prior_sd,
//                                           sigma[j], prior_mean)));
//       
//       
//       // double logz_new = R::rnorm(logz_star, sd_star);
//       double logz_new = rt2(logz_star, sd_star * sd_star, df_t);
//       
//       // findzero_cpp(logz_bar[i,j] - 4,
//       //              logz_bar[i,j] + 4,
//       //              tol = .01,
//       //              x_all[1], y_all, v_all, v_delta, tau[j],
//       //                                                  sigma[j], prior_mean)
//       
//       
//       // double posterior_var = 1 / (1 / prior_var + l2 / lik_var);
//       // double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
//       
//       double loglikelihood = logposterior_logz_cpp(v_samples2, sigma[j],
//                                                    logzbar_current, logz_new,
//                                                    prior_mean, prior_sd * prior_sd);
//       
//       double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
//                                                                           x_all,
//                                                                           v_all,
//                                                                           logzbar_current,
//                                                                           logz_new);
//       
//       double logproposal_ratio = log(dt2(logzbar_current, logz_star, sd_star * sd_star, df_t)) - 
//         log(dt2(logz_new, logz_star, sd_star * sd_star, df_t));
//       // double logproposal_ratio = R::dnorm(logzbar_current, logz_star, sd_star, 1) - 
//       //   R::dnorm(logz_new, logz_star, sd_star, 1);
//       
//       double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
//       
//       if(R::runif(0,1) < exp(logposterior)){
//         
//         logz_bar(i, j) = logz_new;
//         
//       }
//       
//     }
//     
//     sum_m += M_site[i];
//   }
//   
//   List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
//                                        lambda, mu_bar, 
//                                        logz_bar, v_bar, 
//                                        delta,
//                                        gamma, 
//                                        beta_theta,
//                                        M_site, S_star);
//   
//   arma::vec beta02 = list_SP_cpp["beta0"];
//   arma::vec mu2 = list_SP_cpp["mu"];
//   arma::mat logz2 = list_SP_cpp["logz"];
//   arma::mat v2 = list_SP_cpp["v"];
//   
//   
//   return List::create(_["logz"] = logz2,
//                       _["v"] = v2);
//   
// }

arma::mat removeRowCol(arma::mat Sigma, int j){
  
  Sigma.shed_col(j);
  Sigma.shed_row(j);
  
  return(Sigma);
  
}

arma::vec removeRowColj(arma::mat Sigma, int j){
  
  arma::vec Sigmaj = Sigma.col(j);
  Sigmaj.shed_row(j);
  
  return(Sigmaj);
  
}

// // [[Rcpp::export]]
// arma::mat update_logz_NP_corr_cpp(arma::mat logz, arma::vec beta0,
//                                   arma::mat X_z, arma::mat beta_z,
//                                   arma::vec mu, arma::mat v,
//                                   arma::vec lambda, arma::mat beta_theta,
//                                   arma::mat X_w,
//                                   arma::mat beta_w,
//                                   arma::mat Tau, 
//                                   arma::mat delta,
//                                   arma::mat gamma,
//                                   arma::vec sigma, 
//                                   arma::vec M_site,
//                                   int S_star){
//   
//   double df_t = 3;
//   
//   int n = M_site.size();
//   int S = beta0.size();
//   int ncov_z = beta_z.n_rows;
//   int ncov_w_theta = X_w.n_cols;
//   
//   arma::mat Xz_beta = X_z * beta_z;
//   // arma::mat beta0Xz_beta = beta0 + X_z * beta_z;
//   arma::mat Xw_beta = X_w * beta_w;
//   arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
//   
//   if(ncov_w_theta > 0){
//     Xw_beta_theta = X_w * beta_theta.submat(0,2,S - 1,2 + ncov_w_theta);  
//   }
//   
//   arma::mat vtilde_mat = v - Xw_beta;
//   
//   // update parameters
//   
//   for(int j = 0; j < S; j++){
//     
//     arma::mat Sigma_22 = removeRowCol(Tau, j);
//     arma::mat invSigma_22 = arma::inv(Sigma_22);
//     arma::vec Sigma_12 = removeRowColj(Tau, j);
//     
//     arma::vec beta0Xz_beta = beta0[j] + Xz_beta.col(j);
//     
//     int sum_m = 0;
//     for(int i = 0; i < n; i++){
//       
//       double logz_current = logz(i, j);
//       
//       arma::uvec idxes = sum_m + linspace<uvec>(0, M_site[i] - 1, M_site[i]);
//       arma::vec y_all = delta.col(j);
//       y_all = y_all.elem(idxes);
//       arma::vec v_all = beta_theta(j, 0) + Xw_beta_theta.col(j);
//       v_all = v_all.elem(idxes);
//       double x_all = beta_theta(j, 1);
//       
//       arma::vec v_samples = vtilde_mat.col(j);
//       v_samples = v_samples.elem(idxes);
//       arma::uvec idxes_delta = find(y_all == 1);
//       v_samples = v_samples.elem(idxes_delta);
//       
//       // for(int m = 0; m < M_site[i]; m++){
//       //   
//       //   y_all[m] = delta(sum_m + m, j);
//       //   v_all[m] = beta_theta(j, 0) +
//       //     Xw_beta_theta(m + sum_m, j);
//       //   
//       // }
//       
//       arma::vec aminusmu = arma::zeros(S - 1);
//       
//       for(int l = 0; l < (S - 1); l++){
//         if(l < j){
//           aminusmu[l] = logz(i, l) - (beta0[l] + Xz_beta(i, l));
//         } else {   
//           aminusmu[l] = logz(i, l + 1) - (beta0[l + 1] + Xz_beta(i, l + 1));
//         }
//       }
//       
//       double prior_mean_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * aminusmu);
//       
//       double prior_var_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * Sigma_12);
//       
//       double prior_mean = beta0Xz_beta[i] + prior_mean_2;//beta0[j] + Xz_beta(i, j) + prior_mean_2;
//       
//       double prior_sd = sqrt(Tau(j, j) - prior_var_2);
//       
//       double logz_star = findzero_cpp(logz_current - 50,
//                                       logz_current + 50,
//                                       .01,
//                                       x_all,
//                                       y_all, v_all, v_samples, prior_sd,
//                                       sigma[j], prior_mean);
//       
//       double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
//                                           x_all,
//                                           y_all, v_all, v_samples, prior_sd,
//                                           sigma[j], prior_mean)));
//       
//       double logz_new = rt2(logz_star, sd_star, df_t);
//       
//       double loglikelihood = logposterior_logz_cpp(v_samples, sigma[j],
//                                                    logz_current, logz_new,
//                                                    prior_mean, prior_sd * prior_sd);
//       
//       double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
//                                                                           x_all,
//                                                                           v_all,
//                                                                           logz_current,
//                                                                           logz_new);
//       
//       double logproposal_ratio = log(dt2(logz_current, logz_star, sd_star, df_t)) - 
//         log(dt2(logz_new, logz_star, sd_star, df_t));
//       
//       double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;
//       
//       if(R::runif(0,1) < exp(logposterior)){
//         
//         logz(i, j) = logz_new;
//         
//       }
//       
//       sum_m += M_site[i]; 
//     }
//     
//   }
//   
//   return logz;
//   
// }

// [[Rcpp::export]]
arma::mat update_logz_corr_cpp(arma::mat logz, arma::vec beta0,
                               arma::mat X_z, arma::mat beta_z,
                               arma::vec mu, arma::mat v,
                               arma::vec lambda, arma::mat beta_theta,
                               arma::mat X_w,
                               arma::mat beta_w,
                               arma::mat Tau, 
                               arma::mat delta,
                               arma::mat gamma,
                               arma::vec sigma, 
                               arma::vec M_site,
                               int S_star,
                               int emptyTubes){
  
  double df_t = 3;
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  int n = M_site.size();
  int S = beta0.size();
  int ncov_z = beta_z.n_rows;
  int ncov_w_theta = X_w.n_cols;
  
  arma::mat Xz_beta = X_z * beta_z;
  arma::mat Xw_beta = X_w * beta_w;
  arma::mat Xw_beta_theta = arma::zeros(sum(M_site), S);
  
  if(ncov_w_theta > 0){
    Xw_beta_theta = X_w * arma::trans(beta_theta.submat(0,2,S - 1,2 + (ncov_w_theta - 1)));  
  }
  
  // v_bar = v_bar.submat(0,0,sum(M_site) - 1,S - 1);  
  // arma::mat vtilde_mat = v_bar - Xw_beta;
  
  arma::mat vtilde_mat = arma::zeros(sum(M_site), S);
  for(int l = 0; l < sum(M_site); l++){
    for(int j = 0; j < S; j++){
      vtilde_mat(l, j) = v_bar(l, j) - Xw_beta(l, j);
    }
  }
  
  // update parameters
  
  for(int j = 0; j < S; j++){
    
    arma::mat Sigma_22 = removeRowCol(Tau, j);
    arma::mat invSigma_22 = arma::inv(Sigma_22);
    arma::vec Sigma_12 = removeRowColj(Tau, j);
    
    arma::vec beta0Xz_beta = beta_bar[j] + Xz_beta.col(j);
    
    int sum_m = 0;
    for(int i = 0; i < n; i++){
      
      double logz_current = logz_bar(i, j);
      
      arma::uvec idxes = sum_m + linspace<uvec>(0, M_site[i] - 1, M_site[i]);
      arma::vec y_all = delta.col(j);
      y_all = y_all.elem(idxes);
      arma::vec v_all = beta_theta_bar(j) + Xw_beta_theta.col(j);
      v_all = v_all.elem(idxes);
      double x_all = beta_theta(j, 1);
      
      // arma::vec y_all = arma::zeros(M_site[i]);
      // arma::vec v_all = arma::zeros(M_site[i]);
      // // double x_all = beta_theta(j, 1) / exp(lambda[j]);
      // double x_all = beta_theta(j, 1);
      // for(int m = 0; m < M_site[i]; m++){
      //   
      //   y_all[m] = delta(sum_m + m, j);
      //   v_all[m] = beta_theta_bar(j) +
      //     // v_all[m] = beta_theta(j, 0) +
      //     Xw_beta_theta(m + sum_m, j);
      //   
      // }
      
      arma::vec v_samples = vtilde_mat.col(j);
      v_samples = v_samples.elem(idxes);
      arma::uvec idxes_delta = find(y_all == 1);
      v_samples = v_samples.elem(idxes_delta);
      
      // for(int m = 0; m < M_site[i]; m++){
      //   
      //   y_all[m] = delta(sum_m + m, j);
      //   v_all[m] = beta_theta(j, 0) +
      //     Xw_beta_theta(m + sum_m, j);
      //   
      // }
      
      arma::vec aminusmu = arma::zeros(S - 1);
      
      for(int l = 0; l < (S - 1); l++){
        if(l < j){
          aminusmu[l] = logz_bar(i, l) - (beta_bar[l] + Xz_beta(i, l));
        } else {   
          aminusmu[l] = logz_bar(i, l + 1) - (beta_bar[l + 1] + Xz_beta(i, l + 1));
        }
      }
      
      double prior_mean_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * aminusmu);
      
      double prior_mean = beta0Xz_beta[i] + prior_mean_2;
      
      double prior_var_2 = arma::as_scalar(arma::trans(Sigma_12) * invSigma_22 * Sigma_12);
      
      double prior_sd = sqrt(Tau(j, j) - prior_var_2);
      
      double logz_star = findzero_cpp(logz_current - 50,
                                      logz_current + 50,
                                      .01,
                                      x_all,
                                      y_all, v_all, v_samples, prior_sd,
                                      sigma[j], prior_mean);
      
      double sd_star = sqrt(1 / (-h_f_cpp(logz_star,
                                          x_all,
                                          y_all, v_all, v_samples, prior_sd,
                                          sigma[j], prior_mean)));
      
      double logz_new = rt2(logz_star, sd_star, df_t);
      
      double logprior = R::dnorm(logz_new, prior_mean, prior_sd, 1) -
        R::dnorm(logz_current, prior_mean, prior_sd, 1);
      
      double loglikelihood_v = logposterior_logz_cpp(v_samples, sigma[j],
                                                     logz_current, logz_new);
      
      double logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
                                                                          x_all,
                                                                          v_all,
                                                                          logz_current,
                                                                          logz_new);
      
      double logproposal_ratio = log(dt2(logz_current, logz_star, sd_star, df_t)) - 
        log(dt2(logz_new, logz_star, sd_star, df_t));
      
      double logposterior = logprior + loglikelihood_v + 
        logposterior_ratio_logistic + logproposal_ratio;
      
      if(R::runif(0,1) < exp(logposterior)){
        
        logz_bar(i, j) = logz_new;
        
      }
      
      sum_m += M_site[i]; 
    }
    
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, mu_bar,
                                       logz_bar, v_bar,
                                       delta,
                                       gamma,
                                       beta_theta_bar,
                                       M_site, S_star,
                                       emptyTubes);
  arma::mat logz2 = list_SP_cpp["logz"];
  
  return logz2;
  
}

//////////////////////////////////////////
///////// LAMBDA SAMPLER
//////////////////////////////////////////

// [[Rcpp::export]]
double logdpost_cpp(double lambda,
                    NumericVector X_l,
                    double betatheta1,
                    NumericVector Xwbetatheta,
                    arma::vec delta,
                    double beta_barj,
                    double sigma_beta,
                    double mu_barj,
                    double sigma_mu,
                    double lambda_priorj,
                    double sigma_lambda){
  
  NumericVector lmlambda = X_l - lambda;
  
  NumericVector Xbeta = lmlambda * betatheta1  + Xwbetatheta;
  
  double sum1 = //sum(dunif(nonPCRcounts, 0, exp(lambda), 1)) + // double sum1 = sum(dpois(nonPCRcounts, exp(lambda) * lambdatilde_j, 1)) +
    R::dnorm(beta_barj, lambda, sigma_beta, 1) +
    R::dnorm(mu_barj, lambda, sigma_mu, 1) +
    R::dnorm(lambda_priorj, lambda, sigma_lambda, 1);
  
  double sumdbern = 0;
  for(int i = 0; i < Xbeta.size(); i++){
    sumdbern += (- log(1 + exp(-Xbeta[i])));
    if(delta[i] == 0){
      sumdbern += (- Xbeta[i]);
    }
    // sumdbern += R::dbinom(delta[i], 1, logistic(Xbeta[i]), 1);
  }
  
  return(sum1 + sumdbern);
  
}

// // [[Rcpp::export]]
// double logdpost_cpp(double lambda,
//                     NumericVector X_l,
//                     double betatheta1,
//                     NumericVector Xwbetatheta,
//                     arma::vec delta,
//                     NumericVector nonPCRcounts,  // double lambdatilde_j,
//                     double beta_barj,
//                     double sigma_beta,
//                     double mu_barj,
//                     double sigma_mu,
//                     double lambda_priorj,
//                     double sigma_lambda){
//   
//   NumericVector lmlambda = X_l - lambda;
//   
//   NumericVector Xbeta = exp(lmlambda) * betatheta1  + Xwbetatheta ;
//   
//   Xbeta = pmax(pmin(Xbeta, 10),-10);
//   
//   double sum1 = //sum(dunif(nonPCRcounts, 0, exp(lambda), 1)) + // double sum1 = sum(dpois(nonPCRcounts, exp(lambda) * lambdatilde_j, 1)) +
//     R::dnorm(beta_barj, lambda, sigma_beta, 1) +
//     R::dnorm(mu_barj, lambda, sigma_mu, 1) +
//     R::dnorm(lambda_priorj, lambda, sigma_lambda, 1);
//   
//   // for(int i = 0; i < nonPCRcounts.size(); i++){
//   //   if(nonPCRcounts[i] > exp(lambda)){
//   //     sum1 = -exp(50);
//   //   }
//   // }
//   // if(any(nonPCRcounts > exp(lambda))){
//   //   sum1 = -exp(50);
//   // }
//   
//   double sumdbern = 0;
//   for(int i = 0; i < Xbeta.size(); i++){
//     sumdbern += (- log(1 + exp(-Xbeta[i])));
//     if(delta[i] == 0){
//       sumdbern += (- Xbeta[i]);
//     }
//     // sumdbern += R::dbinom(delta[i], 1, logistic(Xbeta[i]), 1);
//   }
//   
//   return(sum1 + sumdbern);
//   
// }

double der2_term1(double lmla,
                  double b,
                  double c){
  
  double term1 = b*(exp(lmla - 2*c - 2*b*exp(lmla)) - b*exp(2*lmla - c - b*exp(lmla)) + exp(lmla - c - b*exp(lmla)));
  double term2 = pow(1 + exp(-c - b*exp(lmla)), 2);
  
  return(term1 / term2);
}

double der2_term2(double lmla,
                  double b,
                  double c){
  
  double term1 = -b*(b*exp(2*lmla - c - b*exp(lmla)) + exp(lmla - c - b*exp(lmla)) + exp(lmla));
  double term2 = pow(1 + exp(-c - b*exp(lmla)), 2);
  
  return(term1 / term2);
}

// [[Rcpp::export]]
double der2_logdpost_cpp(double lambda,
                         NumericVector X_l,
                         double beta_theta1,
                         NumericVector Xwbetatheta,
                         arma::vec delta,
                         double beta_barj,
                         double sigma_beta,
                         double mu_barj,
                         double sigma_mu,
                         double lambda_priorj,
                         double sigma_lambda){
  
  // double sum1 = - exp(lambda) * lambdatilde_j * nonPCRcounts.size() + 
  double sum1 = 
    (- 1 / (sigma_beta * sigma_beta)) +
    (- 1 / (sigma_mu * sigma_mu)) +
    (- 1 / (sigma_lambda * sigma_lambda));
  
  NumericVector lminuslambda = X_l - lambda;
  
  // NumericVector Xbeta = pmax(pmin(lminuslambda, 10),-10);
  
  double sumdbern = 0;
  for(int i = 0; i < lminuslambda.size(); i++){
    sumdbern += ( - beta_theta1 * beta_theta1) * exp(- lminuslambda[i] * beta_theta1 - Xwbetatheta[i]) /
      pow(1 + exp(- lminuslambda[i] * beta_theta1 - Xwbetatheta[i]), 2);
    // if(delta[i] == 1){
    //   sumdbern += der2_term1(lminuslambda[i], beta_theta1, Xwbetatheta[i]);
    // } else {
    //   sumdbern += der2_term2(lminuslambda[i], beta_theta1, Xwbetatheta[i]);
    // }
  }
  
  return(sum1 + sumdbern);
  
}

//////////////////////////////////////////
///////// MU SAMPLER
//////////////////////////////////////////

// [[Rcpp::export]]
arma::vec update_mu_cpp(arma::vec mu, 
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
                        double sigma_mu, 
                        int S_star,
                        int emptyTubes){
  
  int n = M_site.size();
  int S = mu.size();
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
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
                                       beta_theta_bar,
                                       M_site, S_star,
                                       emptyTubes);
  arma::vec mu2 = list_SP_cpp["mu"];
  
  return mu2;
  
}

//////////////////////////////////////////
///////// TAU & SIGMA SAMPLER
//////////////////////////////////////////

// [[Rcpp::export]]
arma::vec update_sigma_cpp(arma::vec sigma, 
                           arma::vec lambda,
                           arma::mat beta_z,
                           arma::vec beta0,
                           arma::vec mu,
                           arma::mat logz,
                           arma::mat v, 
                           arma::mat X_w,
                           arma::mat beta_w,
                           arma::mat delta,
                           arma::mat gamma,
                           arma::mat beta_theta,
                           double a_sigma, 
                           arma::vec b_sigma, 
                           arma::vec M_site,
                           int S_star,
                           int emptyTubes){
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu,
                                       logz, v, delta, 
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  int S = sigma.size();
  int n = M_site.size();
  
  arma::mat Xbeta_w = X_w * beta_w;
  
  for(int j = 0; j < S; j++){
    
    int n_samples = 0;
    double sumsq = 0;
    
    int l = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(delta(l,j) == 1){
          sumsq += pow(v_bar(l,j) - (logz_bar(i,j) + Xbeta_w(l, j)), 2);
          
          n_samples++;
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
arma::vec update_tau_cpp(arma::vec tau, 
                         arma::mat logz, 
                         arma::mat X_z, 
                         arma::mat beta_z, 
                         arma::vec beta0, 
                         double a_tau, 
                         arma::vec b_tau){
  
  int S = beta0.size();
  int n = X_z.n_rows;
  
  for(int j = 0; j < S; j++){
    
    int n_samples = 0;
    double sumsq = 0;
    
    for (int i = 0; i < n; i++) {
      
      double Xz_betaz = arma::as_scalar(X_z.row(i) * beta_z.col(j));
      sumsq += pow((Xz_betaz + beta0[j]) - logz(i,j), 2);
      n_samples += 1;
      
    }
    
    tau[j] = sqrt(rinvgamma_cpp(a_tau + n_samples / 2.0, b_tau[j] + sumsq / 2));
  }
  
  return tau;
}


/// GRAPHICAL HORSESHOE


// [[Rcpp::export]]
arma::mat Matmimi(arma::mat M, int i){
  
  arma::mat M2 = arma::zeros(M.n_cols - 1,
                             M.n_cols - 1);
  
  for(int j1 = 0; j1 < M.n_cols; j1++){
    for(int j2 = 0; j2 < M.n_cols; j2++){
      if(j1 < i){
        if(j2 < i){
          M2(j1, j2) = M(j1, j2);
        } else if(j2 > i) {
          M2(j1, j2 - 1) = M(j1, j2);
        }  
      } else if(j1 > i) {
        if(j2 < i){
          M2(j1 - 1, j2) = M(j1, j2);
        } else if(j2 > i) {
          M2(j1 - 1, j2 - 1) = M(j1, j2);
        }
      }
    }  
  }
  
  return(M2);
}

// [[Rcpp::export]]
arma::mat Matimi(arma::mat M, int i){
  
  arma::mat M2 = arma::zeros(1, M.n_cols - 1);
  
  for(int j1 = 0; j1 < M.n_cols; j1++){
    
    if(j1 < i){
      
      M2(0, j1) = M(i, j1);
      
    } else if(j1 > i) {
      
      M2(0, j1 - 1) = M(i, j1);
      
    }
    
  }
  
  return(M2);
}

// [[Rcpp::export]]
List sample_GraphHorseshoe(int n, 
                           arma::mat S, 
                           List GH_params, 
                           double lambda_Y){
  
  arma::mat Omega = GH_params["Omega"];
  arma::mat Sigma = GH_params["Sigma"];
  arma::mat lambdasq = GH_params["lambdasq"];
  arma::mat nu = GH_params["nu"];
  double csi = GH_params["csi"];
  double tausq = GH_params["tausq"];
  
  int p = Sigma.n_cols;
  
  for(int i = 0; i < p; i++){
    
    double gamm = R::rgamma(n / 2.0 + 1, 2 / (S(i, i) + lambda_Y));
    
    double Sigma_ii = Sigma(i, i);
    arma::mat Sigma_mii = Matimi(Sigma, i); 
    arma::mat Sigma_mimi = Matmimi(Sigma, i);
    
    arma::mat invOmegaii = Sigma_mimi - arma::trans(Sigma_mii) * Sigma_mii / Sigma_ii;
    // Rcout << "here" << std::endl;
    // Sigma[-i,-i] - Sigma[-i,i] %*% t(Sigma[-i,i]) / Sigma[i,i]
    // 
    
    arma::mat lambdasq_mii = Matimi(lambdasq, i);
    arma::mat matrix1 = (S(i, i) + lambda_Y) * invOmegaii;
    for(int j = 0; j < (p - 1); j++){
      matrix1(j,j) += 1 / (lambdasq_mii(0, j) * tausq);
    }
    arma::mat C = arma::inv(matrix1);
    // arma::mat C = arma::inv( + arma::inv(diag(lambdasq[-i,i] * tausq)));
    
    arma::vec beta = mvrnormArma(- C * arma::trans(Matimi(S, i)), C);
    // // beta <- sampleNormFast(S[i,i] * invOmegaii + solve(S[i,i] * Lambda_star * tau^2), - S[-i,i])
    
    if(i == 0){
      // Rcout << beta << std::endl; 
    }
    
    for(int j = 0; j < p; j++){
      if(j < i){
        Omega(j, i) = beta[j];
        Omega(i, j) = beta[j];
      } else if (j > i){
        Omega(j, i) = beta[j - 1];
        Omega(i, j) = beta[j - 1];
      }
    }
    
    double tbetainvomegabeta = arma::as_scalar(arma::trans(beta) * invOmegaii * beta);
    Omega(i, i) = gamm + tbetainvomegabeta;
    // Omega[-i,i] <- beta
    //   Omega[i,-i] <- beta
    // 
    //   Omega[i,i] <- gamm + t(beta) %*% invOmegaii %*% beta
    
    arma::vec Lambda = arma::zeros(p - 1);
    for(int j = 0; j < p; j++){
      if(j < i){
        Lambda[j] = rinvgamma_cpp(1, (1 / nu(j,i)) + beta[j] * beta[j] / (2 * tausq));
        lambdasq(j, i) =  Lambda[j];
        lambdasq(i, j) =  Lambda[j];
      } else if(j > i){
        Lambda[j - 1] = rinvgamma_cpp(1, (1 / nu(j,i)) + beta[j - 1] * beta[j - 1] / (2 * tausq));
        lambdasq(j, i) =  Lambda[j - 1];
        lambdasq(i, j) =  Lambda[j - 1];
      }
    }
    
    //   Lambda_i <- sapply(1:(p-1), function(l){
    //     rinvgamma(1, (1 / nu[-i,i][l]) + beta[l]^2 / (2 * tausq))
    //   })
    //
    //   lambdasq[-i,i] <- Lambda_i
    //
    //   lambdasq[i,-i] <- lambdasq[-i,i]
    //
    for(int j = 0; j < p; j++){
      if(j < i){
        nu(j, i) = rinvgamma_cpp(1, 1 + 1 / Lambda[j]);
        nu(i, j) = nu(j, i);
      } else if(j > i){
        nu(j, i) = rinvgamma_cpp(1, 1 + 1 / Lambda[j - 1]);
        nu(i, j) = nu(j, i);
      }
    }
    // nu[-i,i] <- sapply(1:(p-1), function(l){
    //   rinvgamma(1, 1 + 1 / Lambda_i[l])
    // })
    //
    //   nu[i,-i] <- nu[-i,i]
    //
    
    // recompute Sigma
    //   {
    //     Sigma[-i,-i] <- invOmegaii + (invOmegaii %*% beta) %*% t(invOmegaii %*% beta) / gamm
    //
    //     Sigma[-i,i] <- - (invOmegaii %*% beta) / gamm
    //     Sigma[i,-i] <- Sigma[-i,i]
    //
    //     Sigma[i,i] <- 1 / gamm
    //   }
    
    Sigma(i, i) = 1 / gamm;
    
    arma::mat invOmegaiibeta = invOmegaii * beta;
    Sigma_mii = - invOmegaiibeta / gamm;
    Sigma_mimi = invOmegaii + invOmegaiibeta * arma::trans(invOmegaiibeta) / gamm;
    
    for(int j1 = 0; j1 < p; j1++){
      for(int j2 = 0; j2 < p; j2++){
        if(j1 < i){
          if(j2 < i){
            Sigma(j1, j2) = Sigma_mimi(j1, j2);
          } else if(j2 > i) {
            Sigma(j1, j2) = Sigma_mimi(j1, j2 - 1);
          }  
        } else if(j1 > i) {
          if(j2 < i){
            Sigma(j1, j2) = Sigma_mimi(j1 - 1, j2);
          } else if(j2 > i) {
            Sigma(j1, j2) = Sigma_mimi(j1 - 1, j2 - 1);
          }
        }
      }  
    }
    
    for(int j = 0; j < p; j++){
      if(j < i){
        Sigma(j, i) = Sigma_mii[j];
        Sigma(i, j) = Sigma(j, i);
      } else if(j > i){
        Sigma(j, i) = Sigma_mii[j - 1];
        Sigma(i, j) = Sigma(j, i);
      }
    }
    
  }
  
  // update tau
  
  double sumomegaijlambda = 0;
  for(int j1 = 0; j1 < p; j1++){
    for(int j2 = 0; j2 < j1; j2++){
      sumomegaijlambda += Omega(j1, j2) * Omega(j1, j2) / (lambdasq(j1, j2));
    }  
  }
  
  tausq = rinvgamma_cpp((R::choose(p, 2) + 1) / 2,
                        (1 / csi) + (0.5) * sumomegaijlambda);
  
  csi = rinvgamma_cpp(1, 1 + (1 / tausq));
  
  return List::create(_["Omega"] = Omega,
                      _["lambdasq"] = lambdasq,
                      _["tausq"] = tausq,
                      _["nu"] = nu,
                      _["csi"] = csi,
                      _["Sigma"] = Sigma);
  
}


//////////////////////////////////////////
///////// U SAMPLER
//////////////////////////////////////////

double logpost_u(double u, 
                 arma::vec lambdas, 
                 arma::vec v_current, 
                 arma::vec r_current,
                 double mean_u, 
                 double var_u){
  
  double logpost = 0;
  
  for(int i = 0; i < lambdas.size(); i++){
    logpost += R::dgamma(lambdas[i], r_current[i], 
                         exp(u + v_current[i]) / r_current[i], 
                                                          1);
  }
  
  logpost += R::dnorm(u, mean_u, sqrt(var_u), 1);
  
  return logpost;
}

// [[Rcpp::export]]
List update_u_poisgamma_cpp(arma::mat v,
                            arma::mat u,
                            arma::vec lambda,
                            arma::vec beta0,
                            arma::mat beta_z,
                            arma::mat logz,
                            arma::vec mu,
                            arma::cube lambda_ijk,
                            arma::vec r_nb,
                            arma::mat X_w,
                            arma::mat beta_w,
                            arma::cube c_imk, 
                            arma::mat delta, 
                            arma::mat gamma,
                            double sigma_u, 
                            arma::mat beta_theta,
                            arma::vec sigma,
                            double sigma_gamma, 
                            arma::vec M_site,
                            arma::vec K, 
                            int S_star,
                            int emptyTubes){
  
  double df_t = 3;
  
  int n = M_site.size();
  int S = mu.size();
  
  // UPDATE U CP ------
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta,
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  // UPDATE U CP --------------------------------------------
  
  int lk = 0;
  arma::mat u_vars = arma::zeros(sum(M_site), max(K));
  double sum_uvars = 0;
  double sum_u = 0;
  
  // update parameters
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        
        double u_current = u(l, k);
        
        double a = 0;
        double b = 0;
        arma::vec lambdas = arma::zeros(S + S_star);
        arma::vec v_present = arma::zeros(S + S_star);
        arma::vec r_present = arma::zeros(S + S_star);
        
        int l2 = 0;
        for(int j = 0; j < (S + S_star); j++){
          
          if(c_imk(l, k, j) == 1){
            
            // double a = sum(lambda_i * r[i] / expu_i);
            
            lambdas[l2] = lambda_ijk(l, k, j);
            v_present[l2] = v_bar(l, j);
            r_present[l2] = r_nb[j];
            
            a += lambda_ijk(l, k, j) * r_nb[j] / exp(v_bar(l, j));
            
            b += r_nb[j];
            
            l2 += 1;
            
          }
          
        }
        
        double prior_mean = 0;
        double prior_var = sigma_u * sigma_u;
        
        if(l2 > 0){
          
          arma::vec lambdas2 = arma::zeros(l2);
          arma::vec v_present2 = arma::zeros(l2);
          arma::vec r_present2 = arma::zeros(l2);
          for(int j = 0; j < l2; j++){
            lambdas2[j] = lambdas[j];
            v_present2[j] = v_present[j];
            r_present2[j] = r_present[j];
            
          }
          
          double mu = prior_mean;
          double s = 1 / prior_var;
          double c = mu * s - b;
          
          double log_lam_argument = log(a) - (c / s) - log(s);
          double u_star;
          if(log_lam_argument > 500){ // increment maybe
            u_star = c / s + log_lam_argument - log(log_lam_argument);
          } else {
            u_star = c / s + lambertW0_CS(exp(log_lam_argument));
          }
          
          double var_star = - 1 / (- a * exp(- u_star) - s);
          
          // double u_new = R::rnorm(u_star, sqrt(var_star));
          double u_new = rt2(u_star, sqrt(var_star), df_t);
          
          u_vars(l, k) = var_star;
          sum_uvars += u_vars(l, k);
          
          double logpost_new = logpost_u(u_new, lambdas2, v_present2, r_present2, 
                                         prior_mean, prior_var);
          double logpost_current = logpost_u(u_current,
                                             lambdas2, v_present2, r_present2, 
                                             prior_mean, prior_var);
          
          double logposterior = logpost_new - logpost_current;
          
          double logproposal = log(dt2(u_current, u_star, sqrt(var_star), df_t)) - 
            log(dt2(u_new, u_star, sqrt(var_star), df_t));
          
          if(R::runif(0,1) < exp(logposterior + logproposal)){
            
            u(l, k) = u_new;
          }
          
          sum_u += u(l, k);
          
        } else {
          u(l, k) = R::rnorm(prior_mean, sqrt(prior_var));
          
          u_vars(l, k) = prior_var;
          sum_u += u(l, k);
          sum_uvars += u_vars(l, k);
        }
        
      }
      l += 1;
    }
  }
  
  // // centering
  // // # u <- u - u_var * sum(u) / sum(u_var)
  // l = 0;
  // for(int i = 0; i < n; i++){
  //   for(int m = 0; m < M_site[i]; m++){
  //     for(int k = 0; k < K[l]; k++){
  //       // u(l, k) -= u_vars(l, k) * sum_u / sum_uvars;
  //     }
  //     l += 1;
  //   }
  // }
  
  // recompute lambda and ubar
  sum_u = 0;
  int num_u = 0;
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        sum_u += u(l, k);
        num_u += 1;
      }
      l += 1;
    }
  }
  
  double mean_u = sum_u / num_u;
  
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        // u(l, k) -= mean_u;
      }
      l += 1;
    }
  }
  // l = 0; 
  // for(int i = 0; i < n; i++){
  //   for(int m = 0; m < M_site[i]; m++){
  //     for(int j = 0; j < S; j++){
  //       v_bar(l, j) += mean_u;
  //     }
  //     l += 1;
  //   }
  // }
  for(int j = 0; j < (S + S_star); j++){
    // lambda[j] += mean_u;
    // lambdatilde[j] 
  }
  
  // define parameters in PX
  // arma::mat u_hat = u + zeta;
  // arma::mat v_hat = v_bar - zeta;
  
  // update zeta
  
  // arma::mat Xw_beta = X_w * beta_w;
  // 
  // l = 0;
  // for(int i = 0; i < n; i++){
  // 
  //   arma::mat u_hat = arma::zeros(M_site[i], max(K));
  //   arma::mat v_hat = arma::zeros(M_site[i], S);
  //   arma::vec logz_hat = arma::zeros(S);
  // 
  //   double sum_sigma = 0;
  //   double sum_mu = 0;
  // 
  //   int l2 = 0;
  //   for(int m = 0; m < M_site[i]; m++){
  // 
  //     for(int k = 0; k < K[l + l2]; k++)  {
  // 
  //       u_hat(l2, k) = u(l + l2, k) + zeta[i];
  // 
  //       sum_sigma += (1.0 / pow(sigma_u, 2));
  //       sum_mu += u_hat(l2, k) / pow(sigma_u, 2);
  // 
  //     }
  // 
  //     for(int j = 0; j < S; j++){
  // 
  // 
  // 
  //       if(delta(l, j) == 1){
  // 
  //         v_hat(l2, j) = v_bar(l + l2, j) - zeta[i];
  // 
  //         sum_mu += (logz_bar(i, j) + Xw_beta(l + l2, j) - (v_hat(l2, j))) / pow(sigma[j], 2);
  //         sum_sigma += 1 / pow(sigma[j], 2);
  // 
  //       } else if(gamma(l, j) == 1){
  // 
  //         v_hat(l2, j) = v_bar(l + l2, j) - zeta[i];
  // 
  //         sum_mu += (mu_bar[j] - v_hat(l2, j)) / pow(sigma_gamma, 2);
  //         sum_sigma += 1 / pow(sigma_gamma, 2);
  // 
  //       }
  //     }
  // 
  //     l2 += 1;
  //   }
  // 
  //   for(int j = 0; j < S; j++){
  // 
  //     logz_hat[j] = logz_bar(i, j) - zeta[i];
  // 
  //   }
  // 
  //   double prior_mean = 0;
  //   double prior_var = .5;
  // 
  //   double mean_likelihood = sum_mu;
  //   double sigma_likelihood = 1.0 / sum_sigma;
  // 
  //   double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
  //   double posterior_mean = ((prior_mean / prior_var) + sum_mu) * posterior_var;
  // 
  //   //
  // 
  //   zeta[i] = R::rnorm(posterior_mean, sqrt(posterior_var));
  // 
  //   l2 = 0;
  //   for(int m = 0; m < M_site[i]; m++){
  // 
  //     // reupdate parameters
  //     for(int k = 0; k < K[l + l2]; k++)  {
  // 
  //       u(l + l2, k) = u_hat(l2, k) - zeta[i];
  // 
  //     }
  // 
  //     for(int j = 0; j < S; j++){
  // 
  //       v_bar(l + l2, j) = v_hat(l2, j) + zeta[i];
  // 
  //     }
  // 
  //   }
  // 
  //   for(int j = 0; j < S; j++){
  // 
  //     logz_bar(i, j) = logz_hat[j] + zeta[i];
  // 
  //   }
  // 
  //   l += l2;
  // }
  
  
  // l = 0;
  // for(int i = 0; i < n; i++){
  //   for(int m = 0; m < M_site[i]; m++){
  // 
  //     arma::vec u_hat = arma::zeros(K[l]);
  //     arma::vec v_hat = arma::zeros(S);
  // 
  //     double sum_sigma = 0;
  //     double sum_mu = 0;
  // 
  //     for(int k = 0; k < K[l]; k++)  {
  // 
  //       u_hat[k] = u(l, k) + zeta[l];
  // 
  //       sum_sigma += (1.0 / pow(sigma_u, 2));
  //       sum_mu += u_hat[k] / pow(sigma_u, 2);
  // 
  //     }
  // 
  //     for(int j = 0; j < S; j++){
  // 
  //       if(delta(l, j) == 1){
  // 
  //         v_hat[j] = v_bar(l, j) - zeta[l];
  // 
  //         sum_mu += (logz_bar(i, j) + Xw_beta(l, j) - (v_hat[j])) / pow(sigma[j], 2);//r[l] * alpha[j]
  //         // sum_mu += (logz_bar(i, j) + r[l] * alpha[j] - v_hat(l, j)) / pow(sigma[j], 2);
  //         sum_sigma += 1 / pow(sigma[j], 2);
  //       } else if(gamma(l, j) == 1){
  // 
  //         v_hat[j] = v_bar(l, j) - zeta[l];
  // 
  //         sum_mu += (mu_bar[j] - (v_hat[j])) / pow(sigma_gamma, 2);
  //         // sum_mu += (mu_bar[j] - v_hat(l, j)) / pow(sigma_gamma, 2);
  //         sum_sigma += 1 / pow(sigma_gamma, 2);
  //       }
  //     }
  // 
  //     double prior_mean = 0;
  //     double prior_var = .5;
  // 
  //     double mean_likelihood = sum_mu;
  //     double sigma_likelihood = 1.0 / sum_sigma;
  // 
  //     double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
  //     double posterior_mean = ((prior_mean / prior_var) + sum_mu) * posterior_var;
  // 
  //     //
  // 
  //     zeta[l] = R::rnorm(posterior_mean, sqrt(posterior_var));
  // 
  //     // reupdate parameters
  //     for(int k = 0; k < K[l]; k++)  {
  //       u(l, k) = u_hat[k] - zeta[l];
  //     }
  // 
  //     for(int j = 0; j < S; j++){
  //       v_bar(l, j) = v_hat[j] + zeta[l];
  //     }
  // 
  // 
  //     l += 1;
  //   }
  // }
  
  
  // List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
  //                                      lambda, mu_bar,
  //                                      logz_bar, v_bar,
  //                                      delta,
  //                                      gamma,
  //                                      beta_theta,
  //                                      M_site, 
  //                                      S_star);
  // 
  // arma::vec beta02 = list_SP_cpp["beta0"];
  // arma::vec mu2 = list_SP_cpp["mu"];
  // arma::mat logz2 = list_SP_cpp["logz"];
  // arma::mat v2 = list_SP_cpp["v"];
  
  
  return List::create(_["u"] = u,
                      _["lambda"] = lambda);
  
}

///////////////////////////////
/// V SAMPLER
///////////////////////////////

double logpost_v(double v, 
                 arma::vec lambdas, 
                 arma::vec u_current, 
                 double r,
                 double mean_v, 
                 double var_v){
  
  double logpost = 0;
  
  for(int i = 0; i < lambdas.size(); i++){
    logpost += R::dgamma(lambdas[i], r, 
                         exp(v + u_current[i]) / r, 
                         1);
  }
  
  logpost += R::dnorm(v, mean_v, sqrt(var_v), 1);
  
  return logpost;
}

// [[Rcpp::export]]
arma::mat update_v_poisgamma_cpp(arma::mat v,
                                 arma::mat logz,
                                 arma::vec lambda,
                                 arma::mat X_z,
                                 arma::mat beta_theta,
                                 arma::mat u, arma::mat beta_z,
                                 arma::vec beta0,
                                 arma::vec r_nb,
                                 arma::vec mu, 
                                 arma::cube lambda_ijk,
                                 arma::cube c_imk, arma::mat delta, arma::mat gamma,
                                 arma::vec sigma,
                                 double sigma_gamma,
                                 arma::vec M_site,
                                 arma::mat X_w, 
                                 arma::mat beta_w,
                                 arma::vec K,
                                 int S_star,
                                 int emptyTubes){
  
  
  double df_t = 3;
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta,
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  int n = logz_bar.n_rows;
  int S = mu_bar.size();
  
  // update paramters
  
  arma::mat Xw_beta = X_w * beta_w;
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1 | 
           gamma(l, j) == 1){
          
          double v_current = v_bar(l, j);
          
          double a = 0;
          double b = 0;
          double rnb_current = r_nb[j];
          arma::vec lambdas = arma::zeros(K[l]);
          arma::vec u_present = arma::zeros(K[l]);
          
          int l2 = 0;
          for(int k = 0; k < K[l]; k++){
            
            if(c_imk(l, k, j) == 1){
              
              lambdas[l2] = lambda_ijk(l, k, j);
              u_present[l2] = u(l,k);
              
              a += lambda_ijk(l, k, j) * rnb_current / exp(u(l, k));
              
              b += rnb_current;
              
              l2 += 1;  
              
            }
            
          }
          
          double prior_mean;
          double prior_var;
          
          if(delta(l,j) == 1){
            
            prior_mean = logz_bar(i, j) + Xw_beta(l, j);//alpha[j] * r[l];
            prior_var = sigma[j] * sigma[j];
            
          } else {
            
            prior_mean = mu_bar[j];
            prior_var = sigma_gamma * sigma_gamma;
            
          }
          
          if(l2 > 0){
            
            arma::vec lambdas2 = arma::zeros(l2);
            arma::vec u_present2 = arma::zeros(l2);
            for(int j = 0; j < l2; j++){
              lambdas2[j] = lambdas[j];
              u_present2[j] = u_present[j];
              
            }
            
            double mu = prior_mean;
            double s = 1 / prior_var;
            double c = mu * s - b;
            
            double log_lam_argument = log(a)  - (c / s) - log(s);
            double v_star;
            if(log_lam_argument > 400){
              v_star = c / s + log_lam_argument - log(log_lam_argument);
            } else {
              v_star = c / s + lambertW0_CS(exp(log_lam_argument));
            }
            
            double var_star = - 1 / (- a * exp(- v_star) - s);
            
            double v_new = rt2(v_star, sqrt(var_star), df_t);
            
            double logpost_new = logpost_v(v_new, lambdas2, u_present2, r_nb[j], 
                                           prior_mean, prior_var);
            double logpost_current = logpost_v(v_current,
                                               lambdas2, u_present2, r_nb[j], 
                                                                         prior_mean, prior_var);
            
            double log_posterior = logpost_new - logpost_current;
            
            double log_proposal = log(dt2(v_current, v_star, sqrt(var_star), df_t)) - 
              log(dt2(v_new, v_star, sqrt(var_star), df_t));
            
            if(R::runif(0,1) < exp(log_posterior + log_proposal)){
              v_bar(l, j) = v_new;
            }
            
          } else {
            
            v_bar(l, j) = R::rnorm(prior_mean, sqrt(prior_var));
            
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
                                       beta_theta_bar,
                                       M_site, S_star,
                                       emptyTubes);
  
  arma::vec beta02 = list_SP_cpp["beta0"];
  arma::vec mu2 = list_SP_cpp["mu"];
  arma::mat logz2 = list_SP_cpp["logz"];
  arma::mat v2 = list_SP_cpp["v"];
  
  return v2;
  
}

///////////////////////
////// UV JOINT 
////////////////////////

// [[Rcpp::export]]
double logpost_gamma_prior(double v, 
                           arma::vec lambdas, 
                           arma::vec u_current, 
                           arma::vec r,
                           double mean_v, 
                           double var_v){
  
  double logpost = 0;
  
  for(int i = 0; i < lambdas.size(); i++){
    logpost += R::dgamma(lambdas[i], r[i], 
                         exp(v + u_current[i]) / r[i], 
                                                  1);
  }
  
  logpost += R::dnorm(v, mean_v, sqrt(var_v), 1);
  
  return logpost;
}

// [[Rcpp::export]]
double logpost_gamma(double v, 
                     arma::vec lambdas, 
                     arma::vec u_current, 
                     arma::vec r){
  
  double logpost = 0;
  
  for(int i = 0; i < lambdas.size(); i++){
    logpost += R::dgamma(lambdas[i], r[i], 
                         exp(v + u_current[i]) / r[i], 
                                                  1);
  }
  
  return logpost;
}

// [[Rcpp::export]]
double logpost_gamma_uv(double u, 
                        double v, 
                        arma::vec lambdas, 
                        arma::vec x_current, 
                        arma::vec r,
                        arma::vec lambdas2, 
                        arma::vec x_current2, 
                        arma::vec r2){
  
  double logpost = 0;
  
  for(int i = 0; i < lambdas.size(); i++){
    logpost += R::dgamma(lambdas[i], r[i], 
                         exp(u + v + x_current[i]) / r[i], 
                                                      1);
  }
  
  for(int i = 0; i < lambdas2.size(); i++){
    logpost += R::dgamma(lambdas2[i], r2[i], 
                         exp(u + x_current2[i]) / r2[i], 
                                                    1);
  }
  
  return logpost;
}

arma::mat H_fjoint(double u, 
                   double v, 
                   double a_1, 
                   double sigma_u, 
                   double sigma_v){
  
  arma::mat H = arma::zeros(2, 2);
  H(0, 0) = exp(- v - u) * a_1 - 1 / pow(sigma_u,2);
  H(0, 1) = exp(- v - u) * a_1;
  H(1, 0) = exp(- v - u) * a_1;
  H(1, 1) = exp(- v - u) * a_1 - 1 / pow(sigma_v,2);
  
  return(H);
}

arma::mat H_fjoint_noprior(double u, 
                           double v, 
                           double a_1,
                           double a_3){
  
  arma::mat H = arma::zeros(2, 2);
  H(0, 0) = exp(- v - u) * a_1 + a_3 * exp(-u);
  H(0, 1) = exp(- v - u) * a_1;
  H(1, 0) = exp(- v - u) * a_1;
  H(1, 1) = exp(- v - u) * a_1;
  
  return(H);
}

arma::vec update_coeff_joint(arma::vec xy_current, 
                             arma::vec lambdas, 
                             arma::vec x_present, 
                             arma::vec r,
                             double a_1, 
                             double a_2, 
                             double sigma_x, 
                             double sigma_y, 
                             double m_u, 
                             double m_v){
  
  double log_lamargu = log(-a_1) + log(sigma_x * sigma_x + sigma_y * sigma_y) + 
    (a_2*sigma_x * sigma_x + a_2*sigma_y * sigma_y - m_u - m_v);
  
  double lambert_logargu;
  
  if(log_lamargu > 400){
    lambert_logargu = log_lamargu - log(log_lamargu);
  } else {
    lambert_logargu = lambertW0_CS(exp(log_lamargu));
  }
  
  double x_star = (-a_2*pow(sigma_x, 4) - a_2*pow(sigma_x, 2)*pow(sigma_y, 2) + 
                   lambert_logargu*pow(sigma_x, 2) + m_u*pow(sigma_x, 2) + m_u*pow(sigma_y, 2))/
                     (pow(sigma_x, 2) + pow(sigma_y, 2));
  
  double y_star = (-a_2*pow(sigma_x, 2)*pow(sigma_y, 2) - a_2*pow(sigma_y, 4) + 
                   lambert_logargu*pow(sigma_y, 2) + m_v*pow(sigma_x, 2) + m_v*pow(sigma_y, 2))/(pow(sigma_x, 2) + pow(sigma_y, 2));
  
  arma::vec xy_star = arma::zeros(2);
  xy_star[0] = x_star;
  xy_star[1] = y_star;
  
  arma::mat Hf = H_fjoint(x_star, y_star, a_1, sigma_x, sigma_y);
  
  arma::vec xy_new = mvrnormArma(xy_star, arma::inv(-Hf));
  
  double loglik_new = logpost_gamma(sum(xy_new), lambdas, x_present, r);
  double loglik_old = logpost_gamma(sum(xy_current), lambdas, x_present, r);
  
  double loglik = loglik_new - loglik_old;
  
  double logpr_new = R::dnorm(xy_new[0], m_u, sigma_x, 1) +
    R::dnorm(xy_new[1], m_v, sigma_y, 1); 
  
  double logpr_old = R::dnorm(xy_current[0], m_u, sigma_x, 1) +
    R::dnorm(xy_current[1], m_v, sigma_y, 1);
  
  double logprior = logpr_new - logpr_old;
  
  double logproposal = dmvnorm_cpp(xy_new, xy_star, arma::inv(-Hf), 1) - 
    dmvnorm_cpp(xy_current, xy_star, arma::inv(-Hf), 1);
  
  if(R::runif(0, 1) < exp(loglik + logprior + logproposal)){
    xy_current = xy_new;
  } 
  
  return(xy_current);
}

arma::vec update_coeff_rnb(double x_current,
                           arma::vec lambdas, 
                           arma::vec v_present,
                           arma::vec r,
                           double mean_x, 
                           double sigma_x){
  
  double df_t = 3;
  
  arma::vec out = arma::zeros(2);
  
  if(lambdas.size() > 0){
    
    double a = 0;
    for(int l = 0; l < lambdas.size(); l++){
      a += lambdas[l] * r[l] / exp(v_present[l]);
    }
    
    double b = sum(r);
    
    double mu = mean_x;
    double s = 1 / sigma_x * sigma_x;
    double c_u = mu * s - b;
    
    double log_lam_argument = log(a) - (c_u / s) - log(s);
    
    double u_star;
    if(log_lam_argument > 500){
      u_star = c_u / s + log_lam_argument - log(log_lam_argument);
    } else {
      u_star = c_u / s + lambertW0_CS(exp(log_lam_argument));
    }
    
    double var_star = - 1 / (- a * exp(- u_star) - s);
    
    double x_new = rt2(u_star, sqrt(var_star), df_t);
    
    double logpost_new = logpost_gamma_prior(x_new, lambdas, v_present, r,
                                             mean_x, sigma_x * sigma_x);
    double logpost_current = logpost_gamma_prior(x_current, lambdas, v_present, r,
                                                 mean_x, sigma_x * sigma_x);
    
    double logposterior = logpost_new - logpost_current;
    
    double logproposal = log(dt2(x_current, u_star, sqrt(var_star), df_t)) -
      log(dt2(x_new, u_star, sqrt(var_star), df_t));
    
    if(R::runif(0, 1) < exp(logposterior + logproposal)){
      
      x_current = x_new;
      
    }
    
    out[0] = x_current;
    out[1] = var_star;
    
  } else {
    
    out[0] = R::rnorm(mean_x, sigma_x);
    out[1] = sigma_x * sigma_x;
    
  }
  
  return out;
}

// [[Rcpp::export]]
List update_uv_poisgamma_cpp(arma::mat u,
                             arma::mat v,
                             arma::mat logz,
                             arma::vec lambda,
                             arma::mat X_z,
                             arma::mat beta_theta,
                             arma::mat beta_z,
                             arma::vec beta0,
                             arma::vec r_nb,
                             arma::vec mu, 
                             arma::cube lambda_ijk,
                             arma::cube c_imk, 
                             arma::mat delta, 
                             arma::mat gamma,
                             arma::vec sigma,
                             double sigma_gamma,
                             double sigma_u,
                             arma::vec M_site,
                             arma::mat X_w, 
                             arma::mat beta_w,
                             arma::vec K,
                             int S_star,
                             int emptyTubes){
  
  
  double df_t = 3;
  
  int n = logz.n_rows;
  int S = mu.size();
  
  List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                       gamma, beta_theta, M_site, S_star, emptyTubes);
  arma::vec beta_bar = list_CP_cpp["beta_bar"];
  arma::vec mu_bar = list_CP_cpp["mu_bar"];
  arma::vec beta_theta_bar = list_CP_cpp["beta_theta_bar"];
  arma::mat logz_bar = list_CP_cpp["logz_bar"];
  arma::mat v_bar = list_CP_cpp["v_bar"];
  
  // update paramters
  
  arma::mat Xw_beta = X_w * beta_w;
  
  arma::vec mean_v = arma::zeros(sum(M_site));
  arma::vec var_v = arma::zeros(sum(M_site));
  arma::vec vbar_im = arma::zeros(sum(M_site));
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      int numSamples = 0;
      for(int j = 0; j < S; j++){
        if(delta(l, j) == 1){
          vbar_im[l] += v_bar(l, j);
          mean_v[l] += logz_bar(i, j) + Xw_beta(l, j);
          var_v[l] += sigma[j] * sigma[j];
          numSamples += 1;
        } else if(gamma(l, j) == 1){
          vbar_im(l) += v_bar(l, j);
          mean_v[l] += mu_bar[j];
          var_v[l] += sigma_gamma * sigma_gamma;
          numSamples += 1;
        }
      }
      vbar_im[l] = vbar_im[l] / numSamples;
      mean_v[l] = mean_v[l] / numSamples;
      var_v[l] = var_v[l] / (numSamples * numSamples);
      // var_v[l] = var_v[l] / (numSamples);
      l += 1;
    }
  }
  
  // compute vtilde
  l = 0;
  arma::mat vtilde = arma::zeros(sum(M_site), S + S_star);
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for (int j = 0; j < S; j++) {
        if(delta(l, j) == 1 |
           gamma(l, j) == 1){
          vtilde(l, j) = v_bar(l, j) - vbar_im[l];
        }
      }
      for (int j = S; j < (S + S_star); j++) {
        if(delta(l, j) == 1 |
           gamma(l, j) == 1){
          vtilde(l, j) = v_bar(l, j);
        }
      }
      l += 1;
    }
  }
  
  l = 0;
  arma::vec ubar_im = arma::zeros(sum(M_site));
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      ubar_im[l] = mean(u.row(l));
      l += 1;
    }
  }
  
  // compute utilde
  arma::mat utilde = arma::zeros(sum(M_site), max(K));
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        utilde(l, k) = u(l, k) - ubar_im[l];
      }
      l += 1;
    }
  }
  
  // update ubar vbar
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      // Rcout << "i = " << i << " - m = " << m << std::endl;
      arma::vec uv_current = arma::zeros(2);
      uv_current[0] = ubar_im[l];
      uv_current[1] = vbar_im[l];
      
      double m_v = mean_v[l];
      double m_u = 0;
      double v_v = var_v[l];
      double v_u = sigma_u * sigma_u / (K[l]);
      
      arma::vec mu_prior = arma::zeros(2);
      mu_prior[0] = m_u;
      mu_prior[1] = m_v;
      arma::mat D_var = arma::zeros(2, 2);
      D_var(0, 0) = v_u;
      D_var(1, 1) = v_v;
      
      double a_1 = 0;
      double a_2 = 0;
      double a_3 = 0;
      double a_4 = 0;
      
      arma::vec r_all = arma::zeros(S * K[l]);
      arma::vec lambdas = arma::zeros(K[l] * S);
      arma::vec x_present = arma::zeros(K[l] * S);
      arma::vec r_all_2 = arma::zeros(S_star * K[l]);
      arma::vec lambdas_2 = arma::zeros(K[l] * S_star);
      arma::vec x_present_2 = arma::zeros(K[l] * S_star);
      int l2 = 0;
      int l3 = 0;
      for (int k = 0; k < K[l]; k++) {
        for (int j = 0; j < S; j++) {
          if(c_imk(l, k, j) == 1){
            r_all[l2] = r_nb[j];
            x_present[l2] = vtilde(l, j) + utilde(l, k);
            lambdas[l2] = lambda_ijk(l, k, j);
            a_1 += ( - r_nb[j] * exp(- x_present[l2]) * lambdas[l2]);
            a_2 += r_nb[j];
            l2 += 1;
          }
        }
        for (int j = S; j < (S + S_star); j++) {
          if(c_imk(l, k, j) == 1){
            r_all_2[l3] = r_nb[j];
            x_present_2[l3] = vtilde(l, j) + utilde(l, k);
            lambdas_2[l3] = lambda_ijk(l, k, j);
            a_3 +=( - r_nb[j] * exp(- x_present_2[l3]) * lambdas_2[l3]);
            a_4 += r_nb[j];
            l3 += 1;
          }
        }
      }
      
      arma::vec r_all2 = arma::zeros(l2);
      arma::vec lambdas2 = arma::zeros(l2);
      arma::vec x_present2 = arma::zeros(l2);
      for(int l4 = 0; l4 < l2; l4++){
        x_present2[l4] = x_present[l4];
        lambdas2[l4] = lambdas[l4];
        r_all2[l4] = r_all[l4];
      }
      
      arma::vec r_all2_2 = arma::zeros(l3);
      arma::vec lambdas2_2 = arma::zeros(l3);
      arma::vec x_present2_2 = arma::zeros(l3);
      for(int l4 = 0; l4 < l3; l4++){
        x_present2_2[l4] = x_present_2[l4];
        lambdas2_2[l4] = lambdas_2[l4];
        r_all2_2[l4] = r_all_2[l4];
      }
      
      if(l2 == 0 & l3 == 0){ // no observations present
        
        ubar_im[l] = R::rnorm(m_u, sqrt(v_u));
        vbar_im[l] = R::rnorm(m_v, sqrt(v_v));
        
      } else {
        
        arma::vec mu_star = arma::zeros(2);
        arma::mat Sigma_star;
        
        if(l2 == 0){ // only spike-ins observations available
          
          double log_lamargu = log(-a_1) + log(v_u) + (a_2*v_u);
          
          double lambert_logargu;
          
          if(log_lamargu > 400){
            lambert_logargu = log_lamargu - log(log_lamargu);
          } else {
            lambert_logargu = lambertW0_CS(exp(log_lamargu));
          }
          
          double ubar_star = lambert_logargu - a_2 * v_u;
          double vbar_star = 0;
          
          mu_star[0] = ubar_star;
          mu_star[1] = vbar_star;
          
          Sigma_star = arma::zeros(2, 2);
          Sigma_star(0,0) = 1 / sqrt(-(exp(-ubar_star) * a_1 - 1 / v_u));
          Sigma_star(1,1) = exp(50);
          
        } else if(l3 == 0){ // no spike-ins observations available
          
          double log_lamargu = log(-a_1) + log(v_u + v_v) + (a_2*v_u + a_2*v_v - m_v);
          
          double lambert_logargu;
          
          if(log_lamargu > 400){
            lambert_logargu = log_lamargu - log(log_lamargu);
          } else {
            lambert_logargu = lambertW0_CS(exp(log_lamargu));
          }
          
          double ubar_star = v_u*(-a_2*v_u - a_2*v_v + lambert_logargu)/(v_u + v_v);
          double vbar_star = (-a_2*v_u*v_v - a_2*v_v*v_v +
                              lambert_logargu*v_v + m_v*v_u + m_v*v_v)/(v_u + v_v);
          
          // mu_star = arma::zeros(2);
          mu_star[0] = ubar_star;
          mu_star[1] = vbar_star;
          
          Sigma_star = arma::inv(- H_fjoint(ubar_star, vbar_star, a_1, sqrt(v_u), sqrt(v_v)));
          
        } else {
          
          double ubar_star = -log(-a_4/a_3);
          double vbar_star = log(a_1*a_4/(a_2*a_3));
          arma::vec mu_lik = arma::zeros(2);
          mu_lik[0] = ubar_star;
          mu_lik[1] = vbar_star;
          
          arma::mat Sigma_lik = arma::inv(- H_fjoint_noprior(ubar_star, vbar_star, a_1, a_3));
          
          Sigma_star = arma::inv(arma::inv(Sigma_lik) + arma::inv(D_var));
          
          arma::vec mu = arma::inv(Sigma_lik) * mu_lik + arma::inv(D_var) * mu_prior;
          mu_star = Sigma_star * mu;
          
        }
        
        // Sigma_star = 10 * Sigma_star;
        
        // arma::vec mu_new = mvrnormArma(mu_star, Sigma_star);
        arma::vec mu_new = mrt2(mu_star, Sigma_star, df_t);
        
        double loglik_new = logpost_gamma_uv(mu_new[0], mu_new[1],
                                             lambdas2, x_present2, r_all2,
                                             lambdas2_2, x_present2_2, r_all2_2);
        double loglik_current = logpost_gamma_uv(uv_current[0], uv_current[1],
                                                 lambdas2, x_present2, r_all2,
                                                 lambdas2_2, x_present2_2, r_all2_2);
        
        double loglik = loglik_new - loglik_current;
        
        double logprior_new = R::dnorm(mu_new[0], m_u, sqrt(v_u), 1) +
          R::dnorm(mu_new[1], m_v, sqrt(v_v), 1);
        double logprior_current = R::dnorm(uv_current[0], m_u, sqrt(v_u), 1) +
          R::dnorm(uv_current[1], m_v, sqrt(v_v), 1);
        
        double logprior = logprior_new - logprior_current;
        
        double logproposal_new = dmt_cpp(mu_new, df_t, mu_star, Sigma_star, 1);
        double logproposal_current = dmt_cpp(uv_current, df_t, mu_star, Sigma_star, 1);
        // double logproposal_new = dmvnorm_cpp(mu_new, mu_star, Sigma_star, 1);
        // double logproposal_current = dmvnorm_cpp(uv_current, mu_star, Sigma_star, 1);
        
        double logproposal = logproposal_current - logproposal_new;
        
        if(R::runif(0, 1) < exp(loglik + logprior + logproposal)){
          ubar_im[l] = mu_new[0];
          vbar_im[l] = mu_new[1];
        }
        
      }
      
      
      l += 1;
    }
  }
  
  // recompute v
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for (int j = 0; j < S; j++) {
        if(delta(l, j) == 1 |
           gamma(l, j) == 1){
          v_bar(l, j) = vtilde(l, j) + vbar_im[l];
        }
      }
      l += 1;
    }
  }
  
  // recompute u
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        u(l, k) = utilde(l, k) + ubar_im[l];
      }
      l += 1;
    }
  }
  
  // // update v given vbar (approx)
  // l = 0;
  // for(int i = 0; i < n; i++){
  //   for(int m = 0; m < M_site[i]; m++){
  //     
  //     double sumRow = 0;
  //     int l_j = 0;
  //     arma::vec vars_v = arma::zeros(S);
  //     for (int j = 0; j < S; j++) {
  //       
  //       if(delta(l, j) == 1 |
  //          gamma(l, j) == 1){
  //         
  //         double mean_v = 0;
  //         double sigma_v = 0;
  //         if(delta(l, j) == 1){
  //           mean_v = logz_bar(i, j) + Xw_beta(l, j);
  //           sigma_v = sigma[j];
  //         } else if (gamma(l, j) == 1){
  //           mean_v = mu_bar[j];
  //           sigma_v = sigma_gamma;
  //         }
  //         
  //         double v_current = v_bar(l, j);
  //         
  //         double a_1 = 0;
  //         double a_2 = 0;
  //         
  //         arma::vec r_all = arma::zeros(K[l]);
  //         arma::vec lambdas = arma::zeros(K[l]);
  //         arma::vec x_present = arma::zeros(K[l]);
  //         int l2 = 0;
  //         for (int k = 0; k < K[l]; k++) {
  //           if(c_imk(l, k, j) == 1){
  //             r_all[l2] = r_nb[j];
  //             x_present[l2] = u(l, k);
  //             lambdas[l2] = lambda_ijk(l, k, j);
  //             a_1 += r_nb[j] * exp(- x_present[l2]) * lambdas[l2];
  //             a_2 += r_nb[j];
  //             l2 += 1;
  //           }
  //         }
  //         
  //         arma::vec r_all2 = arma::zeros(l2);
  //         arma::vec lambdas2 = arma::zeros(l2);
  //         arma::vec x_present2 = arma::zeros(l2);
  //         for(int l3 = 0; l3 < l2; l3++){
  //           x_present2[l3] = x_present[l3];
  //           lambdas2[l3] = lambdas[l3];
  //           r_all2[l3] = r_all[l3];
  //         }
  //         
  //         arma::vec newvals = update_coeff_rnb(v_current, lambdas2, x_present2,
  //                                              r_all2, mean_v, sigma_v);
  //         v_bar(l, j) = newvals[0];
  //         vars_v[l_j] = newvals[1];
  //         sumRow += v_bar(l, j);
  //         l_j += 1;
  //       }
  //       
  //     }
  //     
  //     arma::vec vars_v2 = arma::zeros(l_j);
  //     for(int l3 = 0; l3 < l_j; l3++){
  //       vars_v2[l3] = vars_v[l3];
  //     }
  //     
  //     int l3 = 0;
  //     for(int j = 0; j < S; j++){
  //       if(delta(l, j) == 1 |
  //          gamma(l, j) == 1){
  //         v_bar(l, j) =  v_bar(l, j) -
  //           vars_v2[l3] * (sumRow - l_j * vbar_im[l]) / sum(vars_v2);
  //         l3 += 1;
  //       }
  //     }
  //     
  //     l += 1;
  //   }
  // }
  
  // update u given ubar (approx)
  // l = 0;
  // double sum_u = 0;
  // for(int i = 0; i < n; i++){
  //   for(int m = 0; m < M_site[i]; m++){
  //     
  //     arma::vec vars_u = arma::zeros(K[l]);
  //     
  //     for(int k = 0; k < K[l]; k++){
  //       
  //       double mean_u = 0;
  //       
  //       double u_current = u(l, k);
  //       
  //       double a_1 = 0;
  //       double a_2 = 0;
  //       
  //       arma::vec r_all = arma::zeros(S + S_star);
  //       arma::vec lambdas = arma::zeros(S + S_star);
  //       arma::vec x_present = arma::zeros(S + S_star);
  //       int l2 = 0;
  //       for (int j = 0; j < S; j++) {
  //         if(c_imk(l, k, j) == 1){
  //           r_all[l2] = r_nb[j];
  //           x_present[l2] = v_bar(l, j);
  //           lambdas[l2] = lambda_ijk(l, k, j);
  //           // a_1 += r_nb[j] * exp(- x_present[l2]) * lambdas[l2];
  //           // a_2 += r_nb[j];
  //           l2 += 1;
  //         }
  //       }
  //       for (int j = S; j < (S + S_star); j++) {
  //         if(c_imk(l, k, j) == 1){
  //           r_all[l2] = r_nb[j];
  //           x_present[l2] = v_bar(l, j);
  //           lambdas[l2] = lambda_ijk(l, k, j);
  //           // a_1 += r_nb[j] * exp(- x_present[l2]) * lambdas[l2];
  //           // a_2 += r_nb[j];
  //           l2 += 1;
  //         }
  //       }
  //       
  //       arma::vec r_all2 = arma::zeros(l2);
  //       arma::vec lambdas2 = arma::zeros(l2);
  //       arma::vec x_present2 = arma::zeros(l2);
  //       for(int l3 = 0; l3 < l2; l3++){
  //         x_present2[l3] = x_present[l3];
  //         lambdas2[l3] = lambdas[l3];
  //         r_all2[l3] = r_all[l3];
  //       }
  //       //
  //       // lambdas <- as.vector(lambda_ijk[(i - 1)*M + m,k,])
  //       // v_present <- lambda + v[(i - 1)*M + m,]
  //       
  //       arma::vec newvals = update_coeff_rnb(u_current, lambdas2, x_present2,
  //                                            r_all2, mean_u, sigma_u);
  //       
  //       u(l, k) = newvals[0];
  //       vars_u[k] = newvals[1];// 1 / (1 / list_val$u_vars + 1 / sigma_u^2)
  //       
  //     }
  //     
  //     for(int k = 0; k < K[l]; k++){
  //       u(l, k) =  u(l, k) - (vars_u[k] / sum(vars_u)) * (sum(u.row(l)) - K[l] * ubar_im[l]);
  //       sum_u += u(l, k);
  //     }
  //     
  //     l += 1;
  //   }
  // }
  
  // recompute lambda and ubar
  double sum_u = 0;
  int num_u = 0;
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        sum_u += u(l, k);
        num_u += 1;
      }
      l += 1;
    }
  }
  
  double mean_u = sum_u / num_u;
  
  l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        // u(l, k) -= mean_u;
      }
      l += 1;
    }
  }
  l = 0; 
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int j = 0; j < S; j++){
        // v_bar(l, j) += mean_u;
      }
      l += 1;
    }
  }
  // for(int j = S; j < (S + S_star); j++) {
  for(int j = 0; j < (S + S_star); j++) {
    // lambda[j] += mean_u;
    // lambdatilde[j] = lambdatilde[j] / exp(mean_u);
  }
  
  List list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                       lambda, 
                                       mu_bar,
                                       logz_bar, 
                                       v_bar,
                                       delta,
                                       gamma,
                                       beta_theta_bar,
                                       M_site, S_star,
                                       emptyTubes);
  
  arma::mat v2 = list_SP_cpp["v"];
  
  return List::create(
    _["lambda"] = lambda,
    // _["lambdatilde"] = lambdatilde,
    _["u"] = u,
    _["v"] = v2);
  
}

///////////////////////////////
////////////// R SAMPLER 
///////////////////////////////

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

double loglik_r(double r, 
                arma::vec y,
                arma::vec mu_nb){
  
  double sum = 0;
  for(int i = 0; i < y.size(); i++){
    double pi = mu_nb[i] / (mu_nb[i] + r);
    sum += R::dnbinom(y[i], r, 1 - pi, 1);
  }
  
  return sum;
} 


double obj_fun_rcpp(double& r, 
                    arma::vec& y, arma::vec& mean_nb){
  
  double loglik1 = sum(lgamma(r + y)) - y.size() * R::lgammafn(r);
  
  double loglik2 = mean_nb.size() * r * log(r);
  
  double loglik3 = - arma::as_scalar(arma::trans(y + r) * log(mean_nb + r));
  
  double loglik4 = r;
  
  // return (-(loglik1 + loglik2 + loglik3));
  return (-(loglik1 + loglik2 + loglik3 + loglik4));
}

double optim_r_rcpp(double& init_val,
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
                          arma::mat u, 
                          arma::mat v,
                          arma::cube &y,
                          arma::mat delta,
                          arma::mat gamma,
                          arma::cube &c_imk,
                          arma::vec M_site,
                          arma::vec K,
                          bool optimStep,
                          double sd_r_proposal){
  
  // List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
  //                                      gamma, beta_theta, M_site);
  // arma::vec beta_bar = list_CP_cpp["beta_bar"];
  // arma::vec mu_bar = list_CP_cpp["mu_bar"];
  // arma::mat logz_bar = list_CP_cpp["logz_bar"];
  // arma::mat v_bar = list_CP_cpp["v_bar"];
  // arma::mat beta_theta0_bar = list_CP_cpp["beta_theta0_bar"];
  
  int n = M_site.size();
  int S = r_nb.size();
  
  for (int j = 0; j < S; j++) {
    
    double sum_v = 0;
    
    arma::vec y_present0 = arma::zeros(sum(M_site) * max(K));
    arma::vec mean_uv0 = arma::zeros(sum(M_site) * max(K));
    
    int l = 0;
    int l2 = 0;
    for(int l = 0; l < v.n_rows; l++){
      for (int k = 0; k < K[l]; k++) {
        
        if(c_imk(l,k,j) == 1){
          y_present0[l2] = y(l,k,j);
          mean_uv0[l2] = exp(lambda[j] + v(l, j) + u(l, k));
          l2 += 1;
        } 
        
      }
    }
    
    if(l2 > 0){
      
      arma::vec y_present = arma::zeros(l2);
      arma::vec mean_uv = arma::zeros(l2);
      for(int l = 0; l < l2; l++){
        y_present[l] = y_present0[l];
        mean_uv[l] = mean_uv0[l];
      }
      
      double logproposal_diff = 0;
      
      double r_new;
      if(optimStep){
        
        double r_star = optim_r_rcpp(r_nb[j], y_present, mean_uv);
        double var_star = - 1 / d2loglik_r_cpp(r_star, y_present, mean_uv);
        
        // double prior_mean = mean_r;
        // double prior_var = sd_r * sd_r;
        // 
        // double mean_likelihood = r_star;
        // double sigma_likelihood = var_star;
        // 
        // double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
        // double posterior_mean = ((prior_mean / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var;
        // 
        // r_new = R::rnorm(posterior_mean, sqrt(posterior_var));
        
        r_new = R::rnorm(r_star, sqrt(var_star));
        
        logproposal_diff = R::dnorm(r_nb[j], r_star, sqrt(var_star), 1) - 
          R::dnorm(r_star, r_star, sqrt(var_star), 1);
        
      } else {
        
        r_new = R::rnorm(r_nb[j], sd_r_proposal);
        
      }
      
      double loglik_new = loglik_r(r_new, y_present, mean_uv);
      double loglik_current = loglik_r(r_nb[j], y_present, mean_uv);
      
      // double logprior_new = R::dnorm(r_new, mean_r, sd_r, 1);
      // double logprior_current = R::dnorm(r_nb[j], mean_r, sd_r, 1);
      double logprior_new = r_new;//R::dnorm(r_new, mean_r, sd_r, 1);
      // double logprior_new = R::dnorm(r_new, mean_r, sd_r, 1);
      double logprior_current = r_nb[j];//R::dnorm(r_nb[j], mean_r, sd_r, 1);
      // double logprior_current = R::dnorm(r_nb[j], mean_r, sd_r, 1);
      
      double logposterior_new = loglik_new + logprior_new;
      double logposterior_current = loglik_current + logprior_current;
      
      if(R::runif(0, 1) < exp(logproposal_diff + logposterior_new - logposterior_current)){
        r_nb[j] = r_new;
      }
      
    }
    
  }
  
  return(r_nb);
  
}

// struct SampleR : public Worker 
// {
// 
//   // inputs
//   const RVector<double> r_nb;
//   const RVector<double> lambda;
//   const RMatrix<double> u;
//   const RMatrix<double> v;
//   const RMatrix<double> y; // cube
//   const RMatrix<double> c_imk; // cube
//   const RMatrix<double> delta;
//   const RMatrix<double> gamma;
//   const RVector<double> M_site;
//   const RVector<double> K;
//   const RVector<double> sd_r_proposal;
// 
//   // output matrix to write to
//   RVector<double> r_nb_new;
// 
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   SampleR(NumericVector r_nb,
//           NumericVector lambda, 
//           NumericMatrix u, 
//           NumericMatrix v, 
//           NumericMatrix y,
//           NumericMatrix c_imk,
//           NumericMatrix delta,
//           NumericMatrix gamma,
//           NumericVector M_site,
//           NumericVector K,
//           NumericVector sd_r_proposal,
//           NumericVector r_nb_new)
//     : r_nb(r_nb),
//       lambda(lambda), 
//       u(u), 
//       v(v), 
//       y(y), 
//       c_imk(c_imk), 
//       delta(delta),
//       gamma(gamma),
//       M_site(M_site),
//       K(K), 
//       sd_r_proposal(sd_r_proposal), 
//       r_nb_new(r_nb_new) {}
// 
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t j = begin; j < end; j++) {
//       
//       
//       int n = M_site.length();
//       
//       double maxK = *std::max_element(K.begin(), K.end());
// 
//       double sum_v = 0;
// 
//       double sumMsite = std::accumulate(M_site.begin(), M_site.end(), 0.0);
// 
//       NumericVector y_present0(sumMsite * maxK);
//       NumericVector mean_uv0(sumMsite * maxK);
// 
//       int l = 0;
//       int l2 = 0;
//       for(int i = 0; i < n; i++){
//         for(int m = 0; m < M_site[i]; m++){
//           for (int k = 0; k < K[l]; k++) {
// 
//             if(c_imk(l * maxK + k,j) == 1){ // change
//               y_present0[l2] = y(l * maxK + k,j); // change
//               mean_uv0[l2] = exp(lambda[j] + v(l, j) + u(l, k));
//               l2 += 1;
//             }
// 
//           }
//           l += 1;
//         }
//       }
// 
//       
//       if(l2 > 0){
// 
//         arma::vec y_present = arma::zeros(l2);
//         arma::vec mean_uv = arma::zeros(l2);
//         for(int l = 0; l < l2; l++){
//           y_present[l] = y_present0[l];
//           mean_uv[l] = mean_uv0[l];
//         }
//         
//         Rcout << "j = " << j << " - l2 = " << l2 << " / " << std::endl;
// 
//         double logproposal_diff = 0;
// 
//         double r_new  = R::rnorm(r_nb[j], sd_r_proposal[0]);
// 
//         double loglik_new = loglik_r(r_new, y_present, mean_uv);
//         double loglik_current = loglik_r(r_nb[j], y_present, mean_uv);
// 
//         double logprior_new = r_new;
//         double logprior_current = r_nb[j];
//         // 
//         // double logposterior_new = loglik_new + logprior_new;
//         // double logposterior_current = loglik_current + logprior_current;
//         // 
//         // if(R::runif(0, 1) < exp(logproposal_diff + logposterior_new - logposterior_current)){
//         //   r_nb_new[j] = r_new;
//         // }
// 
//         r_nb_new[j] = 2;
//         
//       }
// 
//      
//     }
//   }
// };
// 
// // [[Rcpp::export]]
// NumericVector update_r_nb_cpp_parallel(NumericVector r_nb,
//                                        double mean_r,
//                                        double sd_r,
//                                        NumericVector lambda,
//                                        NumericMatrix u,
//                                        NumericMatrix v,
//                                        NumericMatrix &y,
//                                        NumericMatrix delta,
//                                        NumericMatrix gamma,
//                                        NumericMatrix &c_imk,
//                                        NumericVector M_site,
//                                        NumericVector K,
//                                        bool optimStep,
//                                        NumericVector sd_r_proposal) {
// 
//   // allocate the matrix we will return
//   NumericVector r_nb_new(r_nb.size());
// 
//   // create the worker
//   SampleR sampler(r_nb, lambda, u, v, y, c_imk, delta, gamma, M_site, 
//                   K, sd_r_proposal, r_nb_new);
// 
//   // call it with parallelFor
//   parallelFor(0, r_nb_new.size(), sampler);
// 
//   return r_nb_new;
// }

///////////////////////////////
//// LAMBDA IJK
///////////////////////////////

// [[Rcpp::export]]
arma::cube update_lambdaijk(arma::vec lambda, 
                            arma::cube lambda_ijk,
                            arma::mat v, 
                            arma::mat u,
                            arma::vec r_nb,
                            arma::cube c_imk,
                            arma::vec M_site,
                            arma::cube y,
                            arma::vec K,
                            int S_star,
                            int emptyTubes){
  
  int n = M_site.size();
  int S = y.n_slices - S_star;
  
  arma::cube lambda_ijk2 = arma::zeros(sum(M_site) + emptyTubes, max(K), S + S_star);
  
  int l = 0;
  for(int i = 0; i < n; i++){
    for(int m = 0; m < M_site[i]; m++){
      for(int k = 0; k < K[l]; k++){
        for(int j = 0; j < (S + S_star); j++){
          if(c_imk(l,k,j) == 1){
            
            double mean_lambdaijk = exp(lambda[j] + v(l, j) +
                                        u(l, k));
            
            lambda_ijk2(l, k, j) = 
              R::rgamma(r_nb[j] + y(l, k, j), 
                        1 / (1 + r_nb[j] / mean_lambdaijk) );
            
          } else {
            lambda_ijk2(l, k, j) = NA_REAL;
          }
          
        }
      }
      l += 1;
    }
  }
  
  for(int m = 0; m < emptyTubes; m++){
    for(int k = 0; k < K[l + m]; k++){
      for(int j = 0; j < (S + S_star); j++){
        if(c_imk(l + m,k,j) == 1){
          
          double mean_lambdaijk = exp(lambda[j] + v(l + m, j) +
                                      u(l + m, k));
          
          lambda_ijk2(l + m, k, j) = 
            R::rgamma(r_nb[j] + y(l + m, k, j), 
                      1 / (1 + r_nb[j] / mean_lambdaijk) );
          
        } else {
          lambda_ijk2(l + m, k, j) = NA_REAL;
        }
        
      }
    }
  }
  
  return lambda_ijk2;
}

///// OLD STUFF


// double loglikbetatheta11_cpp(arma::vec y, arma::mat X, arma::vec beta){
//   
//   arma::vec plogistic = logisticXb(X, beta);
//   
//   double loglikelihood = 0;
//   for(int l = 0; l < X.n_rows; l++){
//     loglikelihood += R::dbinom(y[l], 1, plogistic[l], 1);
//   }
//   
//   return(loglikelihood);
// }
// 
// List update_lambdajoint_cpp(arma::vec lambda,
//                             arma::mat v,
//                             arma::mat u, 
//                             arma::cube lambda_ijk,
//                             arma::cube c_imk, 
//                             arma::mat delta, 
//                             arma::mat gamma,
//                             arma::vec sigma,
//                             arma::mat logz,
//                             arma::vec mu,
//                             arma::vec r_nb,
//                             double sigma_gamma,
//                             arma::vec M_site,
//                             arma::mat X_w, 
//                             arma::mat beta_w,
//                             arma::vec lambda_prior,
//                             double sigma_la,
//                             arma::vec K,
//                             int emptyTubes,
//                             int S_star){
//   
//   
//   double df_t = 3;
//   
//   int n = logz.n_rows;
//   int S = mu.size();
//   
//   arma::mat Xw_beta = X_w * beta_w;
//   
//   // obtain centered variables
//   arma::vec vbar_j = arma::zeros(S);
//   arma::vec mean_v = arma::zeros(S);
//   arma::vec var_v = arma::zeros(S);
//   
//   for(int j = 0; j < S; j++){
//     int numSamples = 0;
//     int l = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         
//         if(delta(l, j) == 1){
//           mean_v[j] += logz(i, j) + Xw_beta(l, j);
//           var_v[j] += sigma[j] * sigma[j];
//           vbar_j[j] += v(l, j);
//           numSamples += 1;
//         } else if(gamma(l, j) == 1){
//           mean_v[j] += mu[j];
//           var_v[j] += sigma_gamma * sigma_gamma;
//           vbar_j[j] += v(l, j);
//           numSamples += 1;
//         }
//         l += 1;
//       }
//     }
//     vbar_j[j] = vbar_j[j] / numSamples;
//     mean_v[j] = mean_v[j] / numSamples;
//     var_v[j] = var_v[j] / (numSamples * numSamples);
//   }
//   
//   arma::mat vtilde = arma::zeros(sum(M_site) + emptyTubes, S);
//   for(int j = 0; j < S; j++){
//     int l = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         
//         if(delta(l, j) == 1 | 
//            gamma(l, j) == 1){
//           vtilde(l, j) = v(l, j) - vbar_j[j];
//         }
//         l += 1;
//       }
//     }
//     
//     for(int m = 0; m < emptyTubes; m++){
//       if(delta(l + m, j) == 1 | 
//          gamma(l + m, j) == 1){
//         vtilde(l + m, j) = v(l + m, j) - vbar_j[j];
//       }
//     }
//     
//   }
//   
//   // update parameters
//   for(int j = 0; j < S; j++){
//     
//     double m_v = mean_v[j];
//     double m_la = lambda_prior[j];
//     
//     arma::vec lambdavbar_current = arma::zeros(2);
//     lambdavbar_current[0] = lambda[j];
//     lambdavbar_current[1] = vbar_j[j];
//     
//     double a_1 = 0;
//     double a_2 = 0;
//     
//     arma::vec lambdas = arma::zeros(sum(M_site) * max(K));
//     arma::vec x_present = arma::zeros(sum(M_site) * max(K));
//     arma::vec r_all = arma::zeros(sum(M_site) * max(K));
//     
//     int l = 0;
//     int l2 = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         
//         if(delta(l, j) == 1 | 
//            gamma(l, j) == 1){
//           for(int k = 0; k < K[l]; k++){
//             
//             if(c_imk(l, k, j) == 1){
//               
//               lambdas[l2] = lambda_ijk(l, k, j);
//               x_present[l2] = u(l, k) + vtilde(l, j);
//               r_all[l2] = r_nb[j];
//               
//               a_1 += (- lambda_ijk(l, k, j) *  r_nb[j] / exp(x_present[l2]));
//               a_2 += r_nb[j];
//               
//               l2 += 1;  
//               
//             }
//             
//           }
//           
//         }
//         l += 1;
//       }
//     }
//     
//     arma::vec lambdas2 = arma::zeros(l2);
//     arma::vec x_present2 = arma::zeros(l2);
//     arma::vec r_all2 = arma::zeros(l2);
//     for(int l3 = 0; l3 < l2; l3++){
//       x_present2[l3] = x_present[l3];
//       lambdas2[l3] = lambdas[l3];
//       r_all2[l3] = r_all[l3];
//     }
//     
//     arma::vec xy_current = update_coeff_joint(lambdavbar_current, lambdas2, 
//                                               x_present2, r_all2, a_1, a_2, 
//                                               sigma_la, sqrt(var_v[j]), m_la, m_v);
//     lambda[j] = xy_current[0];
//     vbar_j[j] = xy_current[1];
//     
//     
//   }
//   
//   // recompute parameters
//   
//   for(int j = 0; j < S; j++){
//     int l = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         
//         if(delta(l, j) == 1 | 
//            gamma(l, j) == 1){
//           v(l, j)  = vtilde(l, j) + vbar_j[j];
//         }
//         l += 1;
//       }
//     }
//   }
//   
//   
//   return List::create(_["lambda"] = lambda,
//                       _["v"] = v);
//   
// }
// 
//////////////////////////////////////////
///////// BETA0BAR & SIGMA SAMPLER
//////////////////////////////////////////
// 
// arma::mat update_beta0_bar_cpp(arma::vec beta0_bar,
//                                arma::mat logz_bar, arma::vec tau,
//                                double sigma_beta, arma::vec emptySites){
//   
//   int n = logz_bar.n_rows;
//   int S = beta0_bar.size();
//   
//   // update parameters
//   for(int j = 0; j < S; j++){
//     
//     
//     double sum_logzbar = 0;
//     int numSamples = 0;
//     
//     for(int i = 0; i < n; i++){
//       
//       if(emptySites[i] != 1){
//         
//         sum_logzbar += logz_bar(i, j);
//         
//         numSamples++;
//         
//       }
//       
//     }
//     
//     double lik_mean;
//     double lik_var;
//     if(numSamples > 0){
//       
//       lik_mean = sum_logzbar;
//       lik_var = tau[j] * tau[j];
//       
//     } else {
//       
//       lik_mean = 0;
//       lik_var = exp(50);
//       
//     }
//     
//     double prior_mean = 0;
//     double prior_var = sigma_beta * sigma_beta;
//     
//     double posterior_var = 1 / (1 / prior_var + numSamples / lik_var);
//     double posterior_mean = ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var;
//     
//     beta0_bar[j] = R::rnorm(posterior_mean, sqrt(posterior_var));
//     
//     
//   }
//   
//   return beta0_bar;
//   
// }
// 
// double logposterior_betaz_cpp(arma::vec PCR_counts, arma::vec PCR_v,
//                               arma::colvec beta_current, arma::colvec beta_star,
//                               arma::mat X_z, double prior_var){
//   
//   double loglikelihood = 0;
//   for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
//     
//     double X_zbeta_star = arma::as_scalar(X_z.row(l) * beta_star);
//     double X_zbeta_current = arma::as_scalar(X_z.row(l) * beta_current);
//     
//     loglikelihood += (R::dpois(PCR_counts[l], exp(X_zbeta_star + PCR_v[l]), 1) - 
//       R::dpois(PCR_counts[l], exp(X_zbeta_current + PCR_v[l]), 1));
//     
//   }
//   
//   double logprior = 0;
//   for(int k = 0; k < beta_current.size(); k++){
//     logprior += (R::dnorm(beta_star[k], 0, sqrt(prior_var), 1) - 
//       R::dnorm(beta_current[k], 0, sqrt(prior_var), 1));
//   }
//   
//   return(loglikelihood + logprior);
//   
// }
// 
// double logposterior_betaz_logistic_cpp(arma::vec delta_l, arma::mat delta_x,
//                                        arma::vec delta_c,
//                                        arma::colvec beta_current, arma::colvec beta_star){
//   
//   double loglikelihood = 0;
//   for(int l = 0; (unsigned)l < delta_l.size(); l++){
//     
//     double X_zbeta_star = arma::as_scalar(delta_x.row(l) * beta_star);
//     double X_zbeta_current = arma::as_scalar(delta_x.row(l) * beta_current);
//     
//     double p_current = 1 / (1 + exp(- X_zbeta_current - delta_c[l]));
//     double p_star = 1 / (1 + exp(- X_zbeta_star - delta_c[l]));
//     
//     loglikelihood += (R::dbinom(delta_l[l], 1, p_star, 1) - 
//       R::dbinom(delta_l[l], 1, p_current, 1));
//     
//   }
//   
//   return(loglikelihood);
//   
// }
// 
// arma::mat update_beta_cpp(arma::mat v_bar, double zeta, arma::mat u, 
//                           arma::mat beta_z, arma::vec mu_bar, arma::cube y,
//                           arma::vec r, arma::mat X_z, arma::vec alpha,
//                           arma::mat logz_bar,
//                           arma::cube c_imk, arma::mat delta, arma::mat gamma,
//                           arma::mat Sigma_prop, arma::vec M_site,
//                           double sigma_beta,
//                           arma::mat beta_theta, 
//                           arma::vec K, arma::vec emptySites){
//   
//   int n = M_site.size();
//   int S = alpha.size();
//   int ncov = X_z.n_cols;
//   
//   // update parameters
//   for(int j = 0; j < S; j++){
//     
//     arma::colvec beta_current = beta_z.col(j);
//     
//     arma::vec PCR_counts = arma::zeros(sum(M_site) * max(K));
//     arma::vec PCR_v = arma::zeros(sum(M_site) * max(K));
//     arma::mat X_l = arma::zeros(sum(M_site) * max(K), ncov);
//     
//     int l = 0;
//     
//     int index_m = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         for(int k = 0; k < K[index_m]; k++){
//           
//           if(c_imk(index_m, k, j) == 1){
//             
//             X_l.row(l) = X_z.row(i);
//             PCR_counts[l] = y(index_m, k, j);
//             PCR_v[l] = v_bar(index_m, j) + 
//               alpha[j] * r[index_m] + u(index_m, k);
//             
//             l += 1;
//             
//           }
//           
//         }
//         index_m += 1;
//       }
//       
//     }
//     
//     arma::vec PCR_counts2 = arma::zeros(l);
//     arma::vec PCR_v2 = arma::zeros(l);
//     arma::mat X_l2 = arma::zeros(l, ncov);
//     for(int l2 = 0; l2 < l; l2++){
//       PCR_counts2[l2] = PCR_counts[l2];
//       PCR_v2[l2] = PCR_v[l2];
//       X_l2.row(l2) = X_l.row(l2);
//     }
//     
//     // double mean_prior = 0;
//     // double prior_var = sigma_u * sigma_u;
//     // 
//     // double mean_likelihood;
//     // double sigma_likelihood;
//     // if(l2 > 0){
//     //   
//     //   mean_likelihood = log(sum(PCR_counts2)) - log(sum(exp(PCR_v2)));
//     //   sigma_likelihood = 1.0 / sum(PCR_counts2);
//     //   
//     // } else {
//     //   
//     //   mean_likelihood = 0;
//     //   sigma_likelihood = exp(50);
//     //   
//     // }
//     // 
//     // double posterior_var = 1 / (1 / prior_var + 1 / sigma_likelihood);
//     // double posterior_mean = ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var;
//     
//     arma::vec beta_current_vec = arma::conv_to<arma::vec>::from(beta_current);
//     arma::vec beta_star_vec = mvrnormArma(beta_current_vec, Sigma_prop);
//     
//     arma::colvec beta_star = arma::conv_to<arma::colvec>::from(beta_star_vec);
//     
//     double logposterior_ratio_y = logposterior_betaz_cpp(PCR_counts2, PCR_v2,
//                                                          beta_current, beta_star,
//                                                          X_l2, sigma_beta);
//     
//     arma::vec delta_l = arma::zeros(sum(M_site));
//     arma::vec delta_c = arma::zeros(sum(M_site));
//     arma::mat delta_x = arma::zeros(sum(M_site), ncov);
//     index_m = 0;
//     l = 0;
//     for(int i = 0; i < n; i++){
//       for(int m = 0; m < M_site[i]; m++){
//         
//         delta_l[l] = delta(index_m, j);
//         delta_c[l] = beta_theta(j, 0) + 
//           beta_theta(j, 1) * logz_bar(i, j) +
//           beta_theta(j, 2) * r[index_m];
//         delta_x.row(l) = beta_theta(j, 1) * X_z.row(i);
//         
//         index_m++;
//         l++;
//       }
//     }
//     
//     double logposterior_ratio_logistic = logposterior_betaz_logistic_cpp(delta_l,
//                                                                          delta_x,
//                                                                          delta_c,
//                                                                          beta_current, 
//                                                                          beta_star);
//     
//     double logposterior_ratio = logposterior_ratio_y + logposterior_ratio_logistic;
//     // Rcout << exp(logposterior_ratio) << std::endl;
//     if(R::runif(0,1) < exp(logposterior_ratio)){
//       beta_z.col(j) = beta_star;
//     }
//     
//   }
//   
//   return beta_z;
//   
// }
// 
// // [[Rcpp::export]]
// List createMatricesNleq(int j, arma::mat X_z, arma::mat delta, 
//                         arma::mat beta_theta, arma::mat logz_bar, arma::vec r,
//                         arma::cube c_imk, arma::cube y, arma::vec alpha,
//                         arma::mat u, arma::mat v_bar,
//                         arma::vec M_site, arma::vec K){
//   
//   int n = M_site.size();
//   int ncov = X_z.n_cols;
//   
//   arma::vec y_all = arma::zeros(sum(M_site));
//   arma::vec v_all = arma::zeros(sum(M_site));
//   arma::mat X_all = arma::zeros(sum(M_site), ncov);
//   
//   // data from delta 
//   int index_m = 0;
//   for(int i = 0; i < n; i++){
//     for (int m = 0; m < M_site[i]; m++) {
//       X_all.row(index_m) = beta_theta(j, 1) * X_z.row(i);
//       y_all[index_m] = delta(index_m,j);
//       v_all[index_m] = beta_theta(j, 0) + 
//         beta_theta(j, 1) * logz_bar(i, j) +
//         beta_theta(j, 2) * r[index_m];
//       index_m++;
//     }
//   }
//   
//   // data from y 
//   arma::mat X_all2 = arma::zeros(sum(M_site) * max(K), ncov);
//   arma::vec y_all2 = arma::zeros(sum(M_site) * max(K));
//   arma::vec v_all2 = arma::zeros(sum(M_site) * max(K));
//   index_m = 0;
//   int l = 0;
//   for(int i = 0; i < n; i++){
//     for (int m = 0; m < M_site[i]; m++) {
//       for (int k = 0; k < K[index_m]; k++) {
//         if(c_imk(index_m,k,j) == 1){
//           X_all2.row(l) = X_z.row(i);
//           y_all2[l] = y(index_m,k,j);
//           v_all2[l] = v_bar(index_m,j) +
//             alpha[j] * r[index_m] +
//             u(index_m,k);
//           l++;
//         }
//       }
//       index_m++;
//     }
//   }
//   
//   arma::mat X_all3 = arma::zeros(l, ncov);
//   arma::vec y_all3 = arma::zeros(l);
//   arma::vec v_all3 = arma::zeros(l);
//   for(int l2 = 0; l2 < l; l2++){
//     X_all3.row(l2) = X_all2.row(l2);
//     y_all3[l2] = y_all2[l2];
//     v_all3[l2] = v_all2[l2];
//   }
//   
//   return List::create(_["X_all"] = X_all,
//                       _["y_all"] = y_all,
//                       _["v_all"] = v_all,
//                       _["X_all2"] = X_all3,
//                       _["y_all2"] = y_all3,
//                       _["v_all2"] = v_all3);
// }
// 
// arma::vec log_fp_cpp(arma::vec beta, 
//                      arma::mat X_all, 
//                      arma::vec v_all,
//                      arma::vec y_all,
//                      arma::mat X_all2, 
//                      arma::vec v_all2, 
//                      arma::vec y_all2){
//   
//   int ncov = beta.size();
//   
//   arma::mat Xbeta = X_all * beta;
//   arma::mat Xbeta2 = X_all2 * beta;
//   
//   arma::vec toReturn = arma::zeros(ncov);
//   
//   for(int k = 0; k < ncov; k++){
//     toReturn[k] = -sum(X_all.col(k) % exp(v_all + Xbeta)) + sum(y_all % X_all.col(k)) +
//       (-sum(X_all2.col(k) % exp(v_all2 + Xbeta2)) + sum(y_all2 % X_all2.col(k)));
//   }
//   
//   return(toReturn);
//   
// }
// 
// // [[Rcpp::export]]
// arma::mat H_f_cpp(arma::vec beta, 
//                   arma::mat X_all, 
//                   arma::vec v_all,
//                   arma::vec y_all,
//                   arma::mat X_all2, 
//                   arma::vec v_all2, 
//                   arma::vec y_all2){
//   
//   int ncov = beta.size();
//   
//   arma::mat Xbeta = X_all * beta;
//   arma::mat Xbeta2 = X_all2 * beta;
//   
//   arma::mat H = arma::zeros(ncov, ncov);
//   for(int l1 = 0; l1 < ncov; l1++){
//     for(int l2 = 0; l2 <= l1; l2++){
//       double H_fl1l2 = -sum(X_all.col(l1) % X_all.col(l2) % exp(v_all + Xbeta) / 
//                             ((1 + exp(v_all + Xbeta)) % (1 + exp(v_all + Xbeta))) ) + 
//                             (-sum(X_all2.col(l1) % X_all2.col(l2) % exp(v_all2 + Xbeta2)));
//       H(l1, l2) = H_fl1l2;
//       H(l2, l1) = H(l1, l2);
//     }
//   }
//   
//   return(H);
//   
// }
// 
// double logposterior_beta_cpp(arma::vec beta, 
//                              arma::mat X_all, 
//                              arma::vec v_all,
//                              arma::vec y_all,
//                              arma::mat X_all2, 
//                              arma::vec v_all2, 
//                              arma::vec y_all2){
//   
//   arma::mat Xbeta = X_all * beta;
//   arma::mat Xbeta2 = X_all2 * beta;
//   
//   double loglikelihood = 0;
//   for(int i = 0; i < y_all.size(); i++){
//     loglikelihood += R::dbinom(y_all[i], 1, 1 / (1 + exp(-Xbeta[i] - v_all[i])),
//                                1);
//   }
//   
//   for(int i = 0; i < y_all2.size(); i++){
//     loglikelihood += R::dpois(y_all2[i], exp(Xbeta2[i] + v_all2[i]), 1);
//   }
//   
//   return(loglikelihood);
// }
//
// 
// // [[Rcpp::export]]
// arma::vec update_lambda0_NB_cpp(arma::cube &y, arma::cube &c_imk, 
//                                 double mu0, double n0, 
//                                 double sd_mu0, double sd_n0,
//                                 arma::vec M_site, arma::vec K){
//   
//   int n = M_site.size();
//   int S = y.n_slices;
//   
//   arma::vec nonPCRcounts(sum(M_site) * max(K) * S);
//   int index_m = 0;
//   int l = 0;
//   for(int i = 0; i < n; i++){
//     for(int m = 0; m < M_site[i]; m++){
//       for(int k = 0; k < K[index_m]; k++){
//         for(int j = 0; j < S; j++){
//           if(c_imk(index_m,k,j) == 0){
//             nonPCRcounts[l] = y(index_m,k,j);
//             l++;
//           }
//         }
//       }
//       index_m++;
//     }
//   }
//   
//   NumericVector nonPCRcounts2 = NumericVector(l);
//   for(int i = 0; i < l; i++){
//     nonPCRcounts2[i] = nonPCRcounts[i];
//   }
//   
//   // nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
//   // nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
//   
//   // propose new sets of parameters
//   double mu0_star = R::rnorm(mu0, sd_mu0);
//   double n0_star = R::rnorm(n0, sd_n0);
//   
//   double p0_star = n0_star / (n0_star + mu0_star);
//   double p0 = n0 / (n0 + mu0);
//   
//   if(n0_star > 0){
//     
//     arma::vec lik_star_all = dnbinom(nonPCRcounts2, n0_star, p0_star, 1);
//     double lik_star = sum(lik_star_all);
//     arma::vec lik_current_all = dnbinom(nonPCRcounts2, n0, p0, 1);
//     double lik_current = sum(lik_current_all);
//     
//     if(R::runif(0,1) < exp(lik_star - lik_current)){
//       mu0 = mu0_star;
//       n0 = n0_star;
//     }
//     
//   }
//   
//   arma::vec toReturn(2);
//   toReturn[0] = mu0;
//   toReturn[1] = n0;
//   
//   return(toReturn);
// }
//
// 
// double logposterior_lambda_cpp(arma::vec PCR_counts, arma::vec PCR_v,
//                                double lambda_current, double lambda_star){
//   
//   double loglikelihood = 0;
//   for(int l = 0; (unsigned)l < PCR_counts.size(); l++){
//     
//     loglikelihood += (R::dpois(PCR_counts[l], exp(lambda_star + PCR_v[l]), 1) - 
//       R::dpois(PCR_counts[l], exp(lambda_current + PCR_v[l]), 1));
//     
//   }
//   
//   return(loglikelihood);
//   
// }
// 
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