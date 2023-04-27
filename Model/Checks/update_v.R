Xw_beta = X_w %*% beta_w
i <- 4
m <- 5
l <- 35
m + sum(M_site[seq_len(i-1)])
j <- 1

list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta,
                                gamma, beta_theta, M_site)
beta_bar = list_CP_cpp$beta_bar
mu_bar = list_CP_cpp$mu_bar
logz_bar = list_CP_cpp$logz_bar
v_bar = list_CP_cpp$v_bar

v_current = v_bar[l,j]

a = 0
b = 0
rnb_current = r_nb[j]

lambdas <- lambda_ijk[l,c_imk[l,1:K[l],j] == 1,j]
u_present <- u[l,c_imk[l,1:K[l],j] == 1]
a <- sum(lambdas * rnb_current / exp(u_present))
b <- rnb_current * length(u_present)


if(delta[l,j] == 1){
  
  prior_mean = logz_bar[i,j] + Xw_beta[l,j]
  prior_var = sigma[j] * sigma[j]
  
} else {
  
  prior_mean = mu_bar[j];
  prior_var = sigma_gamma * sigma_gamma;
  
}


mu = prior_mean;
s = 1 / prior_var;
c = mu * s - b

log_lam_argument = log(a)  - (c / s) - log(s)

if(log_lam_argument > 40){
  v_star = c / s + log_lam_argument - log(log_lam_argument)
} else {
  v_star = c / s + lambertW0_CS(exp(log_lam_argument))
}

var_star = - 1 / (- a * exp(- v_star) - s)
var_star =  10 * var_star

# v_new = rnorm(1, v_star, sqrt(var_star))
v_new = rt2(v_star, var_star, .05)

logpost_new = logpost_v(v_new, lambdas, u_present, r_nb[j], 
                        prior_mean, prior_var);
logpost_current = logpost_v(v_current,
                            lambdas, u_present, r_nb[j], 
                            prior_mean, prior_var);

log_posterior = logpost_new - logpost_currentList list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                                                                  gamma, beta_theta, M_site);
arma::vec beta_bar = list_CP_cpp["beta_bar"];
arma::vec mu_bar = list_CP_cpp["mu_bar"];
arma::mat logz_bar = list_CP_cpp["logz_bar"];
arma::mat v_bar = list_CP_cpp["v_bar"];

# log_proposal = dnorm(v_current, v_star, sqrt(var_star), 1) - 
  # dnorm(v_new, v_star, sqrt(var_star), 1)
log_proposal = log(dt2(v_current, mean = v_star, scale = var_star, .05)) - 
  log(dt2(v_new, mean = v_star, scale = var_star, .05))


# log_proposal = dt(v_current, v_star, sqrt(var_star), 1) - 
  # dt(v_new, v_star, sqrt(var_star), 1)

exp(log_posterior + log_proposal)

v_grid <- seq(v_star - 4, v_star + 4, length.out = 100)
fun_values <- sapply(v_grid, function(x){
  logpost_v(x, lambdas, u_present, r_nb[j], 
            prior_mean, prior_var)
})  
lfun_values <- exp(fun_values - max(fun_values))
lfun_values <- lfun_values / sum(lfun_values)

qplot(v_grid, lfun_values) + 
  geom_vline(aes(xintercept = v_true[l,j] + lambda[j])) + 
  geom_vline(aes(xintercept = v_bar[l,j]), color = "red") + 
  geom_vline(aes(xintercept = v_new), color = "green")

ggplot(data = NULL) + 
  stat_function(fun = dt, args = list(df = 1000), color = "black") + 
  stat_function(fun = dnorm, args = list(mean = 0,
                                      sd = sqrt(var_star)), color = "red") +
  xlim(c( - 4,  + 4))

ggplot(data = NULL) + 
  stat_function(fun = Vectorize(dt2), args = list(df = .05,
                                       mean = v_star,
                                       scale = var_star), color = "black") + 
  # stat_function(fun = dnorm, args = list(mean = v_star,
                                      # sd = sqrt(var_star)), color = "red") +
  xlim(c(v_bar[l,j] - 4, v_bar[l,j] + 4))

dt2()
