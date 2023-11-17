library(lmtest)
library(splines)
library(forecast)
library(tidyverse)
library(dplyr)
library(tseries)
library(patchwork)

## Funciones 
source("TS.diag.R")
source("summary.arima.R")
source("salida.arima.R")

# Carga de datos ----
datos <- read.table("norw001x-rwl-noaa.txt", header = TRUE)

datos <- datos %>%
  filter(!is.na(X540011_raw)) 

datos <- datos %>%
  filter(age_CE < 1890) %>% 
  select(age_CE,X540021_raw) %>% 
  rename("tiempo" = age_CE, "ancho" = X540021_raw)

str(datos)

# Volviendo a variables temporales ----
Xt <- ts(datos$ancho, start = 1721, frequency = 1)
t <- as.numeric(time(Xt, start = 1721))


# Estacionalidad ----------------------------------------------------------

per <- LSTS::periodogram(Xt)
plot(per$periodogram ~ per$lambda, type = "l", lwd = 2)

which(per$periodogram  == max(per$periodogram))

2*pi/per$lambda[1] # Periodo d=169, que es el total de datos
# No necesitamos diferenciar (?)
abline(v = 2*pi/169, col = "red")


# Box-Cox -----------------------------------------------------------------

lambda <- forecast::BoxCox.lambda(Xt, method = "guerrero")
plot(forecast::BoxCox(Xt, lambda = lambda), col = "steelblue")

MASS::boxcox()

# Ajuste sin estacionalidad -----------------------------------------------

auto.arima(Xt, d = 0)
# ARMA(1,1)

fit <- forecast::Arima(Xt, c(1,0,1)) #incluye la media, porque ajusta mejor

fit2 <- forecast::Arima(Xt, c(1,0,1), lambda = "auto")

AIC(fit); BIC(fit)
AIC(fit2); BIC(fit2)

# Tests -------------------------------------------------------------------

salida.arima(Xt, fit) # significancia de los coeficientes, normalidad y homocedasticidad. Las weás son normales, son homocedásticas y los coeficientes son significativos

TS.diag(fit$res, lag = 10) 

adf.test(fit$res) # SON ESTACIONARIOS

autoplot(fit) # Valor dentro del circulito
all(abs(coef(fit)) < 1) # Es invertible :D

Box.test(fit$res)

# Queda mejor el modelo sin box-cox
salida.arima(Xt, fit2) # significancia de los coeficientes, normalidad y homocedasticidad. Las weás son normales, son homocedásticas y los coeficientes son significativos

TS.diag(fit2$res, lag = 30) 

adf.test(fit2$res) # SON ESTACIONARIOS

autoplot(fit2) # Valor dentro del circulito
all(abs(coef(fit2)) < 1)


# Predicción --------------------------------------------------------------

source("Durbin_Levinson.R")

## Durbin Levinson
fitted.durbinlevinson <- DurbinLevinson(Xt-mean(Xt), ma = fit$coef[2],
                                        ar = fit$coef[1])$fitted

# Estimación puntual
forecast::forecast(fitted.durbinlevinson, h = 1)$mean + mean(Xt)

# Bandas de confianza
forecast::forecast(fitted.durbinlevinson, h = 1)$lower + mean(Xt)
forecast::forecast(fitted.durbinlevinson, h = 1)$upper + mean(Xt)
# Con un 95% de confianza, está entre 0.61 y 0.72


# ACF ---------------------------------------------------------------------


#### Bandas de confianza:
# rho \in (rho(h) +- z*sqrt(w/n))
Bartlett = function(ma = NULL, ar = NULL, m = 3000, lag.max = 10){
  rho = ARMAacf(ma = ma, ar = ar, lag.max = lag.max+m)
  j = 1:m
  w = c()
  for(h in 1:lag.max){
    w[h] = sum((rho[abs(h+j)+1]+rho[abs(h-j)+1]-2*rho[h+1]*rho[j+1])^2)
  }
  w
}
n <- length(Xt)
w <- Bartlett(ar = coef(fit)[1], ma = coef(fit)[2], lag.max = 22)

# Banda de confianza para el acf empirico
rhoh <- acf(Xt, lag.max = 22)$acf[2:23]

IC_inf <- rhoh - qnorm(1 - 0.05/2)*sqrt(w/n)
IC_sup <- rhoh + qnorm(1 - 0.05/2)*sqrt(w/n)

acf_teo <- ARMAacf(ar = coef(fit)[1], ma = coef(fit)[2], lag.max = 22)
IC <- data.frame(IC_inf = IC_inf, IC_sup = IC_sup)

# ACF teórico
data.frame(ACF = acf_teo, Lag = 0:(length(acf_teo) - 1)) %>%
  ggplot(aes(x = Lag,
             y = ACF)) +
  geom_point() +
  geom_segment(aes(x = Lag, xend = Lag,
                   y = ACF, yend = 0)) +
  geom_hline(yintercept = 1.96/sqrt(169),
             col = "red", linetype = "dashed") +
  geom_hline(yintercept = -1.96/sqrt(169),
             col = "red", linetype = "dashed") +
  geom_line(data = IC, aes(x = 1:22,
                           y = IC_inf),
            col = "blue") +
  geom_line(data = IC, aes(x = 1:22,
                           y = IC_sup),
            col = "blue") +
  labs(y = "ACF teórico", 
       x = "Lag") +
  theme_bw()


# ACF EMPÍRICO
aux <- acf(Xt, plot = FALSE)
data.frame(ACF = aux$acf, Lag = aux$lag) %>%
  ggplot(aes(x = Lag,
             y = ACF)) +
  geom_point() +
  geom_hline(yintercept = 1.96/sqrt(169),
             col = "red", linetype = "dashed") +
  geom_hline(yintercept = -1.96/sqrt(169),
             col = "red", linetype = "dashed") +
  geom_segment(aes(x = Lag, xend = Lag,
                   y = ACF, yend = 0)) +
  labs(y = "ACF empírico", 
       x = "Lag") +
  theme_bw()

## FALTAN LOS INTERVALOS DE CONFIANZA PARA ACF TEÓRICO   -> No sé si está bueno


# Densidad espectral ------------------------------------------------------


spec <- LSTS::spectral.density(ar = fit$coef[1], ma = fit$coef[2],
                               sd = sqrt(fit$sigma2))

#  IC_sup <- spec + 1.96 / sqrt(length(spec))
#  IC_inf <- spec - 1.96 / sqrt(length(spec))
par(mfrow=c(1,2))
plot(spec, type = "l")
# lines(IC_sup, type = "l", col = "red", lty = 2)
# lines(IC_inf, type = "l", col = "red", lty = 2)
plot(per$periodogram, type = "l")

# FALTAN LAS BANDAS DE CONFIANZA x2    <-    en el enunciado no lo dice :V
 

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
fit_whittle <- LS.whittle(
  series = malleco, start = c(coef(fit.1), coef(fit.2)), order = c(p = 1, q = 0),
  ar.order = 1, sd.order = 1, N = 180, n.ahead = 10
)


library(LSTS)
LSTS::LS.whittle(series = Xt,
                 start = c(ar = 1, ma = 1, sigma=1), 
                 order = c(1,1))

summary(fit)
summary(forecast::Arima(Xt, order = c(1, 0, 1), method = "ML")) # Es mejor

# El de whittle gana 4  -  3 al ARMA solito
?forecast::Arima




nlminb(
  start = c(fit$coef[1], abs(fit$coef[2])), N = N,
  objective = LS.whittle.loglik,
  series = x, order = c(p = 1, q = 1)
)
