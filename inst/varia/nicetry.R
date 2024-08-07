library(R.matlab)
data = readMat('../data/DataSW.mat')
name = sapply(1:7, \(i) data$ShortDescr[[i]][[1]])
data = data$y
colnames(data) = name

# data = optimism

spec = specify_bsvarSIGN$new(data, p = 4)

start_time = Sys.time()

spec$prior$estimate_hyper(S = 10000, burn_in = 5000,
                          mu = TRUE, delta = TRUE, lambda = TRUE, psi = TRUE)
post = estimate(spec, S = 100)
class(post) = "PosteriorBSVAR"
irf  = compute_impulse_responses(post, horizon = 40)
plot(irf, probability = 0.68)

end_time = Sys.time()
end_time - start_time

plot.ts(t(spec$prior$hyper))

# hist(spec$prior$hyper[1, ], breaks = 200)
