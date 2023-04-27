library(microbenchmark)

microbenchmark({
  logz <- update_logz_joint_cpp(logz, beta0, X_z, beta_z, mu,
                                v, lambda, beta_theta, X_w, beta_w,
                                Sigma_n, invSigma_n, Tau, delta, gamma, 
                                sigma, M_site, S_star, emptyTubes) 
},{
  list_beta_z <- update_betaz_CP_joint(beta0, beta_z, logz, 
                                       Tau, invSigma_n,
                                       X_z, sigma_beta, 
                                       !beta0equal0)
},{
  v <- update_v_poisgamma_cpp(v, logz,
                              lambda,  X_z,
                              beta_theta, u, beta_z,
                              beta0, r_nb, mu, lambda_ijk,
                              c_imk, delta, gamma, sigma,
                              sigma_gamma, M_site,
                              X_w, beta_w,
                              K, S_star, emptyTubes)
}, times= 5)
