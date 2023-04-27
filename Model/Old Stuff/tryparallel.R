y_mat <- apply(y, 3, function(x){
  c(t(x))
})
c_imk_mat <- apply(c_imk, 3, function(x){
  c(t(x))
})

update_r_nb_cpp_parallel(r_nb, prior_r, sqrt(r_var), lambda, u, 
                                 v, y_mat, delta, gamma, c_imk_mat, 
                                 M_site, K, optimStep = F,
                                 sd_r_proposal = .05)


r_nb <- update_r_nb_cpp(r_nb, prior_r, sqrt(r_var), lambda, u, v, 
                        y, delta, gamma, c_imk, M_site, K, 
                        optimStep = 
                          # T,
                          F,
                        # ((iter - 1) %% 5 == 0),
                        sd_r_proposal = .05)
