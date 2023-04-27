list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta,
                                gamma, beta_theta, M_site)
beta_bar = list_CP_cpp$beta_bar
mu_bar = list_CP_cpp$mu_bar
logz_bar = list_CP_cpp$logz_bar
v_bar = list_CP_cpp$v_bar


i <- 10
j <- 1

#

Xz_beta <- X_z %*% beta_z
Xw_beta <- X_w %*% beta_w_true
Xw_beta_theta <- cbind(1, 
                       rep(exp(logz[,j]), each = M_site),
                       X_w) %*% t(beta_theta)

logzbar_current = logz_bar[i,j]

arma::vec v_samples = arma::zeros(M_site[i]);

v_samples <- v_bar[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] - 
  Xw_beta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]

v_samples <- v_samples[delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1]

y_all <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
v_all <- beta_theta[j,1] + Xw_beta_theta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
x_all = beta_theta[j, 2] / exp(lambda[j])

# for(int m = 0; m < M_site[i]; m++){
#   
#   y_all[m] = delta(sum_m + m, j);
#   v_all[m] = beta_theta(j, 0) +
#     Xw_beta_theta(m + sum_m, j);
#   // r[m + sum_m] * beta_theta(j, 2);
#   
# }

prior_mean = beta_bar[j] + Xz_beta[i,j]
prior_var = tau[j] * tau[j]

logz_star = findzero_cpp(logzbar_current - 50,
                         logzbar_current + 50,
                         .01,
                         x_all,
                         y_all, v_all, v_samples, tau[j],
                         sigma[j], prior_mean);

sd_star = sqrt(1 / (-h_f_cpp(logz_star,
                             x_all,
                             y_all, v_all, v_samples, tau[j],
                             sigma[j], prior_mean)))


logz_new = rt2(logz_star, sd_star * sd_star, .05)


loglikelihood = logposterior_logz_cpp(v_samples, sigma[j],
                                             logzbar_current, logz_new,
                                             prior_mean, prior_var);

logposterior_ratio_logistic = logposterior_logz_logistic_cpp(y_all,
                                                                  x_all,
                                                                    v_all,
                                                                    logzbar_current,
                                                                    logz_new);

double logproposal_ratio = log(dt2(logzbar_current, logz_star, sd_star * sd_star, df_t)) - 
  log(dt2(logz_new, logz_star, sd_star * sd_star, df_t));
// double logproposal_ratio = R::dnorm(logzbar_current, logz_star, sd_star, 1) - 
  //   R::dnorm(logz_new, logz_star, sd_star, 1);

double logposterior = loglikelihood + logposterior_ratio_logistic + logproposal_ratio;

v_grid <- seq(logz_star - 4, logz_star + 4, length.out = 100)
fun_values <- sapply(v_grid, function(x){
  logposterior_logz_cpp(v_samples, sigma[j],
                        logzbar_current, x,
                        prior_mean, prior_var)# +
    # logposterior_logz_logistic_cpp(y_all,
    #                                x_all,
    #                                v_all,
    #                                logzbar_current,
    #                                x)
})  

v_true[1:M_site[i] + sum(M_site[seq_len(i-1)]),j][delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1] + 
  lambda[j]

v_samples

qplot(v_grid, fun_values) + 
  geom_vline(aes(xintercept = logz_true[i,j] + lambda[j])) + 
  geom_vline(aes(xintercept = logz_star), color = "red")


Xw_beta_true <- X_w %*% beta_w_true
logz_true[i,j]
v_true[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] -
   Xw_beta_true[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]

