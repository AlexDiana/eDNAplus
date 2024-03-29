# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dt2 <- function(x, mean, scale, df) {
    .Call(`_eDNAPlus_dt2`, x, mean, scale, df)
}

rt2 <- function(mean, scale, df) {
    .Call(`_eDNAPlus_rt2`, mean, scale, df)
}

mrt2 <- function(mean, Sigma, df) {
    .Call(`_eDNAPlus_mrt2`, mean, Sigma, df)
}

dmt_cpp <- function(x, nu, mu, Sigma, returnLog) {
    .Call(`_eDNAPlus_dmt_cpp`, x, nu, mu, Sigma, returnLog)
}

K2 <- function(x1, x2, a, l) {
    .Call(`_eDNAPlus_K2`, x1, x2, a, l)
}

rpg <- function(n, z) {
    .Call(`_eDNAPlus_rpg`, n, z)
}

lambertW0_CS <- function(x) {
    .Call(`_eDNAPlus_lambertW0_CS`, x)
}

lambertWm1_CS <- function(x) {
    .Call(`_eDNAPlus_lambertWm1_CS`, x)
}

sample_beta_cpp <- function(X, B, b, Omega, k) {
    .Call(`_eDNAPlus_sample_beta_cpp`, X, B, b, Omega, k)
}

sample_betaPG <- function(beta, X, b, B, n, k) {
    .Call(`_eDNAPlus_sample_betaPG`, beta, X, b, B, n, k)
}

dmvnorm_cpp <- function(data, m, Sigma, returnLog) {
    .Call(`_eDNAPlus_dmvnorm_cpp`, data, m, Sigma, returnLog)
}

rmtrnorm <- function(mu, U, V) {
    .Call(`_eDNAPlus_rmtrnorm`, mu, U, V)
}

rmtrnorm_chol <- function(mu, A, B) {
    .Call(`_eDNAPlus_rmtrnorm_chol`, mu, A, B)
}

update_betatheta11_cpp <- function(logz, beta_theta11, theta11, delta, X_w, M_site, b_theta11, B_theta11, updateBetaTheta0, updateBetaTheta1) {
    .Call(`_eDNAPlus_update_betatheta11_cpp`, logz, beta_theta11, theta11, delta, X_w, M_site, b_theta11, B_theta11, updateBetaTheta0, updateBetaTheta1)
}

rinvgamma_cpp <- function(a, b) {
    .Call(`_eDNAPlus_rinvgamma_cpp`, a, b)
}

update_betaw_cpp <- function(beta_w, v, delta, logz, X_w, sigma, sigma_beta, M_site) {
    .Call(`_eDNAPlus_update_betaw_cpp`, beta_w, v, delta, logz, X_w, sigma, sigma_beta, M_site)
}

convertSPtoCP_cpp <- function(lambda, beta_z, beta0, mu, logz, v, delta, gamma, beta_theta, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_convertSPtoCP_cpp`, lambda, beta_z, beta0, mu, logz, v, delta, gamma, beta_theta, M_site, S_star, emptyTubes)
}

convertCPtoSP_cpp <- function(beta0_bar, lambda, mu_bar, logz_bar, v_bar, delta, gamma, beta_theta_bar, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_convertCPtoSP_cpp`, beta0_bar, lambda, mu_bar, logz_bar, v_bar, delta, gamma, beta_theta_bar, M_site, S_star, emptyTubes)
}

logposterior_logz_cpp <- function(v_samples, sigma, logz_current, logz_star) {
    .Call(`_eDNAPlus_logposterior_logz_cpp`, v_samples, sigma, logz_current, logz_star)
}

logposterior_logz_logistic_cpp <- function(y_all, x_all, v_all, logz_current, logz_star) {
    .Call(`_eDNAPlus_logposterior_logz_logistic_cpp`, y_all, x_all, v_all, logz_current, logz_star)
}

logf_cpp <- function(l, x, sigmaj, v_samples, v, y, tauj, prior_mean) {
    .Call(`_eDNAPlus_logf_cpp`, l, x, sigmaj, v_samples, v, y, tauj, prior_mean)
}

findzero_cpp <- function(a, b, tol, x_all, y_all, v_all, v_samples, tauj, sigma_j, prior_mean) {
    .Call(`_eDNAPlus_findzero_cpp`, a, b, tol, x_all, y_all, v_all, v_samples, tauj, sigma_j, prior_mean)
}

h_f_cpp <- function(l, x, y_all, v_all, v_samples, tauj, sigma_j, prior_mean) {
    .Call(`_eDNAPlus_h_f_cpp`, l, x, y_all, v_all, v_samples, tauj, sigma_j, prior_mean)
}

update_logz_cpp <- function(logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, tau, delta, gamma, sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_logz_cpp`, logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, tau, delta, gamma, sigma, M_site, S_star, emptyTubes)
}

update_logz_corr_cpp <- function(logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Tau, delta, gamma, sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_logz_corr_cpp`, logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Tau, delta, gamma, sigma, M_site, S_star, emptyTubes)
}

kronProdVec <- function(A, B, v) {
    .Call(`_eDNAPlus_kronProdVec`, A, B, v)
}

update_logz_joint_cpp <- function(logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, Sigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_logz_joint_cpp`, logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, Sigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes)
}

update_logz_joint_fast_cpp <- function(logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, Sigma_S, chol_invSigma_n, chol_invSigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_logz_joint_fast_cpp`, logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, Sigma_S, chol_invSigma_n, chol_invSigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes)
}

update_logz_joint_speed_cpp <- function(logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, invSigma_n, Sigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_logz_joint_speed_cpp`, logz, beta0, X_z, beta_z, mu, v, lambda, beta_theta, X_w, beta_w, Sigma_n, invSigma_n, Sigma_S, delta, gamma, sigma, M_site, S_star, emptyTubes)
}

logdpost_cpp <- function(lambda, X_l, betatheta1, Xwbetatheta, delta, beta_barj, sigma_beta, mu_barj, sigma_mu, lambda_priorj, sigma_lambda) {
    .Call(`_eDNAPlus_logdpost_cpp`, lambda, X_l, betatheta1, Xwbetatheta, delta, beta_barj, sigma_beta, mu_barj, sigma_mu, lambda_priorj, sigma_lambda)
}

logdpost_cpp_beta0 <- function(lambda, X_l, betatheta1, Xwbetatheta, Xbetalogz, delta, logz_barj, tau, mu_barj, sigma_mu, lambda_priorj, sigma_lambda) {
    .Call(`_eDNAPlus_logdpost_cpp_beta0`, lambda, X_l, betatheta1, Xwbetatheta, Xbetalogz, delta, logz_barj, tau, mu_barj, sigma_mu, lambda_priorj, sigma_lambda)
}

der2_logdpost_cpp <- function(lambda, X_l, beta_theta1, Xwbetatheta, delta, beta_barj, sigma_beta, mu_barj, sigma_mu, lambda_priorj, sigma_lambda) {
    .Call(`_eDNAPlus_der2_logdpost_cpp`, lambda, X_l, beta_theta1, Xwbetatheta, delta, beta_barj, sigma_beta, mu_barj, sigma_mu, lambda_priorj, sigma_lambda)
}

der2_logdpost_cpp_beta0 <- function(lambda, X_l, beta_theta1, Xwbetatheta, Xbetalogz, delta, logz_barj, tau, mu_barj, sigma_mu, lambda_priorj, sigma_lambda) {
    .Call(`_eDNAPlus_der2_logdpost_cpp_beta0`, lambda, X_l, beta_theta1, Xwbetatheta, Xbetalogz, delta, logz_barj, tau, mu_barj, sigma_mu, lambda_priorj, sigma_lambda)
}

update_mu_cpp <- function(mu, lambda, delta, gamma, sigma, sigma_gamma, beta0, beta_z, logz, v, beta_theta, M_site, sigma_mu, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_mu_cpp`, mu, lambda, delta, gamma, sigma, sigma_gamma, beta0, beta_z, logz, v, beta_theta, M_site, sigma_mu, S_star, emptyTubes)
}

update_sigma_cpp <- function(sigma, lambda, beta_z, beta0, mu, logz, v, X_w, beta_w, delta, gamma, beta_theta, a_sigma, b_sigma, M_site, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_sigma_cpp`, sigma, lambda, beta_z, beta0, mu, logz, v, X_w, beta_w, delta, gamma, beta_theta, a_sigma, b_sigma, M_site, S_star, emptyTubes)
}

update_tau_cpp <- function(tau, logz_tilde, a_tau, b_tau) {
    .Call(`_eDNAPlus_update_tau_cpp`, tau, logz_tilde, a_tau, b_tau)
}

Matmimi <- function(M, i) {
    .Call(`_eDNAPlus_Matmimi`, M, i)
}

Matimi <- function(M, i) {
    .Call(`_eDNAPlus_Matimi`, M, i)
}

sample_GraphHorseshoe <- function(n, S, GH_params, lambda_Y) {
    .Call(`_eDNAPlus_sample_GraphHorseshoe`, n, S, GH_params, lambda_Y)
}

update_u_poisgamma_cpp <- function(v, u, offsets, lambda, beta0, beta_z, logz, mu, lambda_ijk, r_nb, X_w, beta_w, c_imk, delta, gamma, sigma_u, beta_theta, sigma, sigma_gamma, M_site, K, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_u_poisgamma_cpp`, v, u, offsets, lambda, beta0, beta_z, logz, mu, lambda_ijk, r_nb, X_w, beta_w, c_imk, delta, gamma, sigma_u, beta_theta, sigma, sigma_gamma, M_site, K, S_star, emptyTubes)
}

update_v_poisgamma_cpp <- function(v, logz, lambda, X_z, beta_theta, u, offsets, beta_z, beta0, r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, M_site, X_w, beta_w, K, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_v_poisgamma_cpp`, v, logz, lambda, X_z, beta_theta, u, offsets, beta_z, beta0, r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, M_site, X_w, beta_w, K, S_star, emptyTubes)
}

logpost_gamma_prior <- function(v, lambdas, u_current, r, mean_v, var_v) {
    .Call(`_eDNAPlus_logpost_gamma_prior`, v, lambdas, u_current, r, mean_v, var_v)
}

logpost_gamma <- function(v, lambdas, u_current, r) {
    .Call(`_eDNAPlus_logpost_gamma`, v, lambdas, u_current, r)
}

logpost_gamma_uv <- function(u, v, lambdas, x_current, r, lambdas2, x_current2, r2) {
    .Call(`_eDNAPlus_logpost_gamma_uv`, u, v, lambdas, x_current, r, lambdas2, x_current2, r2)
}

update_uv_poisgamma_cpp <- function(u, v, offsets, logz, lambda, X_z, beta_theta, beta_z, beta0, r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, sigma_u, M_site, X_w, beta_w, K, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_uv_poisgamma_cpp`, u, v, offsets, logz, lambda, X_z, beta_theta, beta_z, beta0, r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, sigma_u, M_site, X_w, beta_w, K, S_star, emptyTubes)
}

update_r_nb_cpp <- function(r_nb, lambda, u, offsets, v, y, delta, gamma, c_imk, M_site, K, mean_r, sd_r, optimStep, sd_r_proposal, S, S_star) {
    .Call(`_eDNAPlus_update_r_nb_cpp`, r_nb, lambda, u, offsets, v, y, delta, gamma, c_imk, M_site, K, mean_r, sd_r, optimStep, sd_r_proposal, S, S_star)
}

update_lambdaijk <- function(lambda, lambda_ijk, v, u, offsets, r_nb, c_imk, M_site, y, K, S_star, emptyTubes) {
    .Call(`_eDNAPlus_update_lambdaijk`, lambda, lambda_ijk, v, u, offsets, r_nb, c_imk, M_site, y, K, S_star, emptyTubes)
}

dnbinom_mean <- function(x, n, mu) {
    .Call(`_eDNAPlus_dnbinom_mean`, x, n, mu)
}

sample_cpp <- function(x, probs) {
    .Call(`_eDNAPlus_sample_cpp`, x, probs)
}

compute_logprob_y_delta0_cpp <- function(y_counts, c_imk_current, currentK, n0, mu0, pi0, n_tilde, mu_tilde, lambda) {
    .Call(`_eDNAPlus_compute_logprob_y_delta0_cpp`, y_counts, c_imk_current, currentK, n0, mu0, pi0, n_tilde, mu_tilde, lambda)
}

DecToBin_cpp <- function(l, currentK) {
    .Call(`_eDNAPlus_DecToBin_cpp`, l, currentK)
}

compute_logprob_y_delta1_rnb_cpp <- function(y_counts, c_imk_current, currentK, n0, mu0, pi0, r_nb, v_im, n_tilde, mu_tilde, lambda, u_im) {
    .Call(`_eDNAPlus_compute_logprob_y_delta1_rnb_cpp`, y_counts, c_imk_current, currentK, n0, mu0, pi0, r_nb, v_im, n_tilde, mu_tilde, lambda, u_im)
}

update_delta_c_d_rjmcmc <- function(v_pres, y, v, lambda, r_nb, M_site, K, mu0, n0, pi0, mu_tilde, n_tilde, u, offsets, logz, X_w, beta_w, sigma, mu, sigma_gamma, v_sd, p11, p10, theta11, theta10, spikedSample, emptyTubes, S_star) {
    .Call(`_eDNAPlus_update_delta_c_d_rjmcmc`, v_pres, y, v, lambda, r_nb, M_site, K, mu0, n0, pi0, mu_tilde, n_tilde, u, offsets, logz, X_w, beta_w, sigma, mu, sigma_gamma, v_sd, p11, p10, theta11, theta10, spikedSample, emptyTubes, S_star)
}

convertDeltaIndexes <- function(delta, gamma, c, K) {
    .Call(`_eDNAPlus_convertDeltaIndexes`, delta, gamma, c, K)
}

convertIndexToDeltaGammaC <- function(index, K) {
    .Call(`_eDNAPlus_convertIndexToDeltaGammaC`, index, K)
}

update_delta_c_d_proposals <- function(v_pres, c_imk, delta, gamma, A, y, v, lambda, r_nb, M_site, K, mu0, n0, pi0, mu_tilde, n_tilde, u, logz, X_w, beta_w, sigma, mu, sigma_gamma, v_sd, p11, p10, theta11, theta10, spikedSample, emptyTubes, S_star) {
    .Call(`_eDNAPlus_update_delta_c_d_proposals`, v_pres, c_imk, delta, gamma, A, y, v, lambda, r_nb, M_site, K, mu0, n0, pi0, mu_tilde, n_tilde, u, logz, X_w, beta_w, sigma, mu, sigma_gamma, v_sd, p11, p10, theta11, theta10, spikedSample, emptyTubes, S_star)
}

update_theta10_cpp <- function(theta_10, delta, gamma, M_site, a0, b0) {
    .Call(`_eDNAPlus_update_theta10_cpp`, theta_10, delta, gamma, M_site, a0, b0)
}

update_p_11_cpp <- function(p_11, delta, gamma, c_imk, M_site, K, a0, b0) {
    .Call(`_eDNAPlus_update_p_11_cpp`, p_11, delta, gamma, c_imk, M_site, K, a0, b0)
}

update_p_10_cpp <- function(p_10, delta, gamma, c_imk, M_site, K, a_p1, b_p1, emptyTubes) {
    .Call(`_eDNAPlus_update_p_10_cpp`, p_10, delta, gamma, c_imk, M_site, K, a_p1, b_p1, emptyTubes)
}

