library(R.matlab)
data = readMat('../data/DataSW.mat')
name = sapply(1:7, \(i) data$ShortDescr[[i]][[1]])
data = data$y
colnames(data) = name

# data = optimism

spec = specify_bsvarSIGN$new(data, p = 4)

start_time = Sys.time()

spec$prior$estimate_hyper(TRUE, TRUE, TRUE, TRUE, S = 20000, burn = 10000, start = 2000)
# post = estimate(spec, S = 1000)
# irf  = compute_impulse_responses(post, horizon = 40)
# plot(irf, probability = 0.68)

end_time = Sys.time()
end_time - start_time

plot.ts(t(spec$prior$hyper))




data(optimism)

# optimism shock
# no contemporaneous effect on productivity
zero_irf          = matrix(0, nrow = 5, ncol = 5)
zero_irf[1, 1]    = 1
# positive contemporaneous effect on stock prices
sign_irf          = array(0, dim = c(5, 5, 1))
sign_irf[2, 1, 1] = 1

specification     = specify_bsvarSIGN$new(optimism,
                                          p        = 4,
                                          sign_irf = sign_irf,
                                          zero_irf = zero_irf)
specification$prior$estimate_hyper(FALSE, TRUE, TRUE, TRUE)

posterior         = estimate(specification, S = 10000)
irf               = compute_impulse_responses(posterior, horizon = 40)
plot(irf, probability = 0.68)