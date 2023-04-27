a_t <- y_all # 1
b_t <- x_all #.05
c_t <- v_all # 0
x <- logz_star_mean[j] #9

a_t <- 1
b_t <- .05
c_t <- .01
x <- 3

term1 <- sum(a_t * b_t * exp(x))

term2_2 = sum(b_t * b_t * exp(2 * (b_t * exp(x) + c_t + x)) /
                             ((1 + exp(c_t + b_t * exp(x)))^2  ) )

b_t * b_t * (exp( (b_t * exp(x) + c_t + x)) /
                             ((1 + exp(c_t + b_t * exp(x)))  ) )^2

b_t * b_t * (   exp( (x )) /
                             ((1 + 1 / exp(c_t + b_t * exp(x)))  )  )^2
 (   exp( (x + log(b_t))) /
                             ((1 + 1 / exp(c_t + b_t * exp(x)))  )  )^2


b_t * (exp(x) * b_t + 1) * exp( x) / (1 + 1 / exp(c_t + exp(x) * b_t))
 

(b_t * exp(x)) / (1 + (1 / exp(b_t * exp(x) + c_t)))

  // double term2_3 = sum(x_all * (x_all_l + 1) * exp(x_all_l + v_all + l) / (1 + exp(v_all + x_all_l)));