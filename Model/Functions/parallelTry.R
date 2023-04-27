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
  
},{
  v <- update_v_poisgamma_cpp(v, logz, lambda,  X_z, beta_theta, u, beta_z,
                              beta0, r_nb, mu, lambda_ijk,
                              c_imk, delta, gamma, sigma, sigma_gamma, M_site, 
                              X_w, beta_w,
                              K, emptySites)
}, times = 10)
