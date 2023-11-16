# Whittle -----------------------------------------------------------------

N <- 169
S <- 100
M <- trunc((length(Xt) - N) / S + 1)
table <- c()
for (j in 1:M) {
  x <- Xt[(1 + S * (j - 1)):(N + S * (j - 1))]
  table <- rbind(table, nlminb(
    start = c(0, 0), N = N,
    objective = LS.whittle.loglik,
    series = x, order = c(p = 1, q = 1)
  )$par)
}
u <- (N / 2 + S * (1:M - 1)) / length(malleco)
table <- as.data.frame(cbind(u, table))
colnames(table) <- c("u", "phi", "sigma")
# Start parameters
phi <- smooth.spline(table$phi,
                     spar = 1,
                     tol = 0.01)$y
fit.1 <- nls(phi ~ a0 + a1 * u, start = list(a0 = 0.65, a1 = 0.00))
sigma <- smooth.spline(table$sigma, spar = 1)$y
fit.2 <- nls(sigma ~ b0 + b1 * u, start = list(b0 = 0.65, b1 = 0.00))
fit_whittle <- LS.whittle(series = malleco, start = c(coef(fit.1), coef(fit.2)), order = c(p = 1, q = 0),  ar.order = 1, sd.order = 1, N = 180, n.ahead = 10)


library(LSTS)
LSTS::LS.whittle(series = Xt,
                 start = c(ar = 1, ma = 1), 
                 ar.order = 1, ma.order = 1)

summary(fit)
summary(forecast::Arima(Xt, order = c(1, 0, 1), method = "ML")) # Es mejor

# El de whittle gana 4  -  3 al ARMA solito
?forecast::Arima




nlminb(
  start = c(fit$coef[1], abs(fit$coef[2])), N = N,
  objective = LS.whittle.loglik,
  series = x, order = c(p = 1, q = 1)
  )

whittle.loglik = function(x, serie, p = 1, q = 1){
  Y = serie
  n = length(Y)
  ar  = x[1:p]
  ma  = x[(p+1):(p+q)]
  if(p == 0){ar = numeric()}
  if(q == 0){ma = numeric()}
  sigma = x[p+q+1]
  aux = LSTS::periodogram(Y, plot = F)
  I = aux$periodogram
  lambda = aux$lambda
  f = LSTS::spectral.density(ar = ar, ma = ma, sd = sigma, lambda = lambda)
  aux = 0.5*(sum(log(f)) + sum(I/f))/n
  aux
}

whittle <- nlminb(start = c(0.5, 0.5,0.5), objective = whittle.loglik, serie = Xt, p = 1, q = 1)$par

fit_whittle <-  forecast::Arima(Xt-mean(Xt), order = c(1,0,1), fixed = whittle[1:2], include.mean = F)

fitted_whittle.auto.arima <- fit_whittle$fitted

sum_w <- summary(fit_whittle)
