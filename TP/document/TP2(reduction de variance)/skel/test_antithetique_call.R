## Fixons la graine pour reproduire les memes resultats qu'avant
 set.seed(123)
S0 <- 100
K <- 0.95
sigma <- 0.25
r <- 0.02
t <- 0

## Calcul du prix theorique
dplus <- (log(S0 / K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t))
dminus <- (log(S0 / K) + (r - 0.5 * sigma^2) * t) / (sigma * sqrt(t))
call.theo <- S0 * pnorm(dplus) - exp(-r * t) * K * pnorm(dminus)

n.sim <- 10^4
Z <- rnorm(n.sim, 0, 1)
St <- S0 * exp((r - 0.5 * sigma^2) * t + sigma * sqrt(t) * Z)
St.anti <- S0 * exp((r - 0.5 * sigma^2) * t - sigma * sqrt(t) * Z)

call <- exp(-r * t) * pmax(St - K, 0)
call.anti <- exp(-r * t) * pmax(St.anti - K, 0)

 cat(cor(call, call.anti))

precision.anti <- 1.96 * sqrt(0.5 * var(call) / n.sim *   (1 + cor(call, call.anti)))

c(call.theo, mean(0.5 * (call + call.anti)), precision.anti)
precision.anti 

