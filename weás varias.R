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

acf_teo <- ARMAacf(ar = fit$coef[1], ma = fit$coef[2], lag.max = 22)

# Gráfico piñufla del acf estimado vs empírico
par(mfrow=c(1,2))
plot(acf_teo, type = "h")
abline(a = 1.96/sqrt(169), b = 0)
acf(Xt)
# Creo que hay weás raras aquí, cuidao

## FALTAN LOS INTERVALOS DE CONFIANZA PARA ACF TEÓRICO


# Densidad espectral ------------------------------------------------------

par(mfrow=c(1,2))
spec <- LSTS::spectral.density(ar = fit$coef[1], ma = fit$coef[2], 
                               sd = sqrt(fit$sigma2))
plot(spec, type = "l")
plot(per$periodogram, type = "l")

# FALTAN LAS BANDAS DE CONFIANZA x2


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