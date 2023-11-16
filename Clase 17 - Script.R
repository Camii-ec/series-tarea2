#####################
## Script Clase 17 ##
#####################

## Ejemplo: 
X <- LSTS::malleco

par(mfrow = c(1,1), bty = "n", las = 1)
plot(X, xlim = c(1200, 2000), ylim = c(0,2), xlab = "", ylab = "", lwd = 2)
acf(X, lag.max = 10, ylim = c(-1,+1), main = "", lwd = 3)

## 4. Algotirmo de Durbin - Levinson (EMV)

source("Durbin_Levinson.R")

## Ej: LSTS::malleco

X <- LSTS::malleco
n <- length(X)
plot(X)
par(mfrow = c(1,2))
ACF <- acf(X, ylim = c(-1,+1), lag.max = 10)
pacf(X, ylim = c(-1,+1), lag.max = 10, xlim = c(0,10))

## Estimación de Máxima Verosimilitud

## Estimador de Máxima Verosimilitud (Durbin Levinson)
dl.loglik = function(x, serie, p = 1, q = 1){
  ar = NULL
  ma = NULL 
  if(p > 0){
    ar  = x[1:p]
  }
  if(q > 0){
    ma  = x[(p+1):(p+q)]
  }
  fit   = DurbinLevinson(serie = serie, ar = ar, ma = ma)
  e     = serie-fit$fitted
  nu    = fit$nu 
  aux = 0.5*(sum(log(fit$nu)) +  sum(e^2/fit$nu))
  if(is.nan(aux)){aux = Inf}
  aux
}

fit <- nlminb(start = c(0.8), objective = dl.loglik, serie = X-mean(X), p = 1, q = 0, lower = c(-0.99999), upper = c(+0.99999))
fit

############################################
## Método de Máxima Verosimilitid Whittle ##
############################################

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
coef(fit)
nlminb(start = c(fit$coef[1], fit$coef[2]),
       objective = whittle.loglik, 
       serie = Xt, p = 1, q = 1)$par


############################
## Inferencia Estadística ##
############################

## Gamma teorica para AR(1):
f = function(x,phi){
  ((-2*cos(x)+2*phi)/(1-2*phi*cos(x)+phi^2))^2/(4*pi)
}

integrate(f, lower = -pi, upper = +pi, phi = 0.5)$value
phi = 0.5
1/(1-phi^2)

## Gamma teórica para ARMA(1,1)

var.coef.arma11 <- function(phi, theta){
  f11 = function(x, theta){
    aux = (2*cos(x)+2*theta)/(1+2*theta*cos(x)+theta^2)
    aux = aux*aux/(4*pi)
    aux
  }
  
  f22 = function(x, phi){
    aux = (-2*cos(x)+2*phi)/(1-2*phi*cos(x)+phi^2)
    aux = aux*aux/(4*pi)
    aux
  }
  
  f12 = function(x, phi, theta){
    aux1 = (-2*cos(x)+2*phi)/(1-2*phi*cos(x)+phi^2)
    aux2 = (2*cos(x)+2*theta)/(1+2*theta*cos(x)+theta^2)
    aux = aux1*aux2/(4*pi)
    aux
  }
  
  G11 = integrate(f11, lower = -pi, upper = pi, theta = theta)$value
  1/(1-theta^2)
  
  G22 = integrate(f22, lower = -pi, upper = pi, phi = phi)$value
  1/(1-phi^2)
  
  G12 = G21 = integrate(f12, lower = -pi, upper = pi, phi = phi, theta = theta)$value
  -1/(1+phi*theta)
  
  G = matrix(c(G22,G12,G21,G11), ncol = 2, byrow = T)
  solve(G)
}

var.coef.arma11(phi = 0.5, theta = 0.8)

## Aplicación Malleco:
n = length(X)

## AR(1)
fit <- forecast::Arima(X-mean(X), order = c(1,0,0), include.mean = F)
fit

## coef
phi <- fit$coef[1]
## s.e.
se <- sqrt((1-phi^2)/n)
## z value
Z0 = fit$coef[1]/se
## p value
2*(1-pnorm(abs(Z0)))

## ARMA(1,1)
fit <- forecast::Arima(X-mean(X), order = c(1,0,1), include.mean = F)
fit

## z value
Z0 <- fit$coef/sqrt(diag(fit$var))
## p value
2*(1-pnorm(abs(Z0)))

source("summary.arima.R")

summary_arima(fit, fixed = c(NA,NA))
