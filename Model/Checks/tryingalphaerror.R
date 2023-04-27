j <- 1

true_v_diff <- v_true[,j] - rep(logz_true[,j], each = M_site)
v_diff <- v[,j] - rep(logz[,j], each = M_site)

ggplot(data = NULL) +  geom_point(aes(x = true_v_diff[delta[,j] == 1],
                    y = v_diff[delta[,j] == 1])) + 
  ylim(c(-6,6)) + 
  xlim(c(-6,6)) + stat_function(fun = function(x) x)

ggplot(data = NULL) +  geom_point(aes(x = v_true[,j],
                    y =  v[,j])) +
  ylim(c(-6,6)) + 
  xlim(c(-6,6)) + stat_function(fun = function(x) x)

ggplot(data = NULL
       ) +  geom_point(aes(x = logz_true[,j],
                           y = logz[,j])) +
  ylim(c(-6,6)) +
  xlim(c(-6,6)) + stat_function(fun = function(x) x)



