library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)

ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)

v_cbind <- foreach(m = 1:S, .combine=cbind, 
                   .packages = "functionsForForeach") %dopar% {
                     
                     v_current <- functionsForForeach::update_v_poisgamma_cpp(v[,m,drop=F], 
                                                                              logz[,m,drop=F], 
                                                                              lambda[m],  
                                                                              X_z, 
                                                                              beta_theta[m,,drop=F], 
                                                                              u, 
                                                                              beta_z[,m,drop=F],
                                                                              beta0[m], 
                                                                              r_nb[m], 
                                                                              mu[m], 
                                                                              lambda_ijk[,,m,drop=F],
                                                                              c_imk[,,m,drop=F], 
                                                                              delta[,m,drop=F], 
                                                                              gamma[,m,drop=F], 
                                                                              sigma[m], 
                                                                              sigma_gamma, 
                                                                              M_site, 
                                                                              X_w, 
                                                                              beta_w[,m,drop=F],
                                                                              K, emptySites)
                   }


v <- update_v_poisgamma_cpp(v, logz, lambda,  X_z, beta_theta, u, beta_z,
                            beta0, r_nb, mu, lambda_ijk,
                            c_imk, delta, gamma, sigma, sigma_gamma, M_site, 
                            X_w, beta_w,
                            K, emptySites)

library(microbenchmark)

microbenchmark({
  list_delta_all <- foreach(m = 1:S, .combine=c, 
                            .packages = "functionsForForeach") %dopar% {
                              
                              list_delta <- functionsForForeach::update_delta_c_d_rjmcmc(y[,,m,drop=F], 
                                                                                         v[,m,drop=F], 
                                                                                         lambda[m], 
                                                                                         r_nb[m],
                                                                                         M_site, K, 
                                                                                         lambdatilde[m],
                                                                                         mu0, n0, pi0, u, 
                                                                                         logz[,m,drop=F], X_w, 
                                                                                         beta_w[,m,drop=F],
                                                                                         sigma[m], 
                                                                                         mu[m], sigma_gamma, v_sd = .5,
                                                                                         p_11[m], 
                                                                                         p_10[m], 
                                                                                         theta11[,m,drop=F], 
                                                                                         theta10[m], emptySites)
                            }
  
  delta <- Reduce("cbind",list_delta_all[1 + 4 * 0:(S-1)])
  c_imk <- Reduce("abind",list_delta_all[2 + 4 * 0:(S-1)])
  gamma <- Reduce("cbind",list_delta_all[3 + 4 * 0:(S-1)])
  v <- Reduce("cbind",list_delta_all[4 + 4 * 0:(S-1)])
  
},{
  list_deltagammac <- update_delta_c_d_rjmcmc(y, v, lambda, r_nb,
                                              M_site, K, lambdatilde,
                                              mu0, n0, pi0, u, logz, X_w, beta_w,
                                              sigma, mu, sigma_gamma, v_sd = .5,
                                              p_11, p_10, theta11, theta10, emptySites)
}, times = 10)
