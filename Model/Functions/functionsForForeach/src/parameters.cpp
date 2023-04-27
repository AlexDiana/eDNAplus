#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;

const double log2pi = std::log(2.0 * M_PI);

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

//

double dt2(double x, double mean, double scale, double df){
    double tstat = (x - mean) / scale;
    return(R::dt(tstat, df, 0) / scale);
}

double rt2(double mean, double scale, double df){
    double tstat = R::rt(df);
    return(tstat * scale + mean);
}

//

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

List convertCPtoSP_cpp(arma::mat beta0_bar,
                       arma::vec lambda, arma::vec mu_bar, 
                       arma::mat logz_bar, arma::mat v_bar, 
                       arma::mat delta, arma::mat gamma, 
                       arma::mat beta_theta_bar,
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
            logz(i,j) = logz_bar(i,j) - lambda[j];
        }  
    }
    
    arma::mat v = arma::zeros(sum(M_site), S);
    
    int l = 0;
    for (int i = 0; i < n; i++) {
        for (int m = 0; m < M_site[i]; m++) {
            for (int j = 0; j < S; j++) {
                if(delta(l, j) == 1){
                    v(l, j) = v_bar(l, j) - lambda[j];
                } else if (gamma(l,j) == 1){
                    v(l,j) = v_bar(l,j) - lambda[j];
                } 
                
            }
            l += 1;
        }
    }
    
    arma::mat beta_theta = arma::zeros(S, beta_theta_bar.n_cols);
    for(int j = 0; j < S; j++){
        beta_theta(j, 0) = beta_theta_bar(j, 0);
        beta_theta(j, 1) = beta_theta_bar(j, 1) * exp(lambda[j]);
        for(int k = 0; k < (beta_theta.n_cols - 2); k++){
            beta_theta(j, k + 2) = beta_theta_bar(j, k + 2);
        }
    }
    
    return List::create(_["beta0"] = beta0,
                        _["mu"] = mu,
                        _["logz"] = logz,
                        _["v"] = v,
                        _["beta_theta"] = beta_theta);
    
}

List convertSPtoCP_cpp(arma::vec lambda, 
                       arma::mat beta_z,
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
            logz_bar(i,j) = logz(i,j) + lambda[j];
        }  
    }
    
    arma::mat v_bar = arma::zeros(sum(M_site), S);
    
    int l = 0;
    for (int i = 0; i < n; i++) {
        for (int m = 0; m < M_site[i]; m++) {
            for (int j = 0; j < S; j++) {
                if(delta(l, j) == 1){
                    v_bar(l, j) = v(l, j) + lambda[j];
                } else if (gamma(l,j) == 1){
                    v_bar(l,j) = v(l,j) + lambda[j];
                } 
                
            }
            l += 1;
        }
    }
    
    arma::mat beta_theta_bar = arma::zeros(S, beta_theta.n_cols);
    for(int j = 0; j < S; j++){
        beta_theta_bar(j, 0) = beta_theta(j, 0);
        beta_theta_bar(j, 1) = beta_theta(j, 1) / exp(lambda[j]);
        for(int k = 0; k < (beta_theta.n_cols - 2); k++){
            beta_theta_bar(j, k + 2) = beta_theta(j, k + 2);
        }
    }
    
    return List::create(_["beta_bar"] = beta_bar,
                        _["beta_theta_bar"] = beta_theta_bar,
                        _["mu_bar"] = mu_bar,
                        _["logz_bar"] = logz_bar,
                        _["v_bar"] = v_bar);
    
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
                                 arma::vec K, arma::vec emptySites){
    
    
    double df_t = .05;
    
    List list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta,
                                         gamma, beta_theta, M_site);
    arma::vec beta_bar = list_CP_cpp["beta_bar"];
    arma::vec mu_bar = list_CP_cpp["mu_bar"];
    arma::mat logz_bar = list_CP_cpp["logz_bar"];
    arma::mat v_bar = list_CP_cpp["v_bar"];
    
    int n = logz_bar.n_rows;
    int S = mu_bar.size();
    
    // update paramters
    
    arma::mat Xw_beta = X_w * beta_w;
    
    int l = 0;
    for(int i = 0; i < n; i++){
        for(int m = 0; m < M_site[i]; m++){
            if(!(emptySites[i] == 1)){
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
                            if(log_lam_argument > 40){
                                v_star = c / s + log_lam_argument - log(log_lam_argument);
                            } else {
                                v_star = c / s + lambertW0_CS(exp(log_lam_argument));
                            }
                            
                            // double v_star = c / s + log_lam_argument - log(log_lam_argument);
                            
                            double var_star = - 1 / (- a * exp(- v_star) - s);
                            // var_star = 5 * var_star;
                            
                            // double v_new = R::rnorm(v_star, sqrt(var_star));
                            double v_new = rt2(v_star, var_star, df_t);
                            
                            double logpost_new = logpost_v(v_new, lambdas2, u_present2, r_nb[j], 
                                                           prior_mean, prior_var);
                            double logpost_current = logpost_v(v_current,
                                                               lambdas2, u_present2, r_nb[j], 
                                                                                         prior_mean, prior_var);
                            
                            double log_posterior = logpost_new - logpost_current;
                            
                            double log_proposal = log(dt2(v_current, v_star, var_star, df_t)) - 
                                log(dt2(v_new, v_star, var_star, df_t));
                            // double log_proposal = R::dnorm(v_current, v_star, sqrt(var_star), 1) - 
                            //   R::dnorm(v_new, v_star, sqrt(var_star), 1);
                            
                            if(R::runif(0,1) < exp(log_posterior + log_proposal)){
                                v_bar(l, j) = v_new;
                            }
                            
                        } else {
                            
                            v_bar(l, j) = R::rnorm(prior_mean, sqrt(prior_var));
                            
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
                                         beta_theta,
                                         M_site);
    
    arma::vec beta02 = list_SP_cpp["beta0"];
    arma::vec mu2 = list_SP_cpp["mu"];
    arma::mat logz2 = list_SP_cpp["logz"];
    arma::mat v2 = list_SP_cpp["v"];
    
    return v2;
    
}

//

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
