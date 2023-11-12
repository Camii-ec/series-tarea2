# Librerías ----
library(lmtest)
library(splines)
library(forecast)
library(tidyverse)
library(dplyr)
library(tseries)

## Funciones 
source("TS.diag.R")
source("summary.arima.R")
source("salida.arima.R")

# Carga de datos ----
datos <- read.table("norw001x-rwl-noaa.txt", header = TRUE)

datos <- datos %>%
  filter(!is.na(X540011_raw)) 

datos <- datos %>%
  filter(age_CE < 1890)

# Gráfico de serie de tiempo ----
datos %>%
  ggplot(aes(x = age_CE, y = X540021_raw)) +
  geom_line(color = "brown4")+
  geom_hline(yintercept = mean(datos$X540021_raw), col = "blue3") +
  labs(x = "Año", y = "Ancho anillo") +
  theme_bw()

# ACF y PACF ----

acf(datos$X540021_raw, main = "ACF")
pacf(datos$X540021_raw, main = "PACF")

# Volviendo a variables temporales ----
Xt <- ts(datos$X540021_raw, start = 1721, frequency = 1)
t <- as.numeric(time(Xt, start = 1721)) # ayuda no empieza en 1721 aa

# Tests ----
ks.test(Xt,"pnorm") # Rechaza normalidad y existen "ties"
bptest(Xt ~ t) # no son homocedasticos

lambda <- forecast::BoxCox.lambda(Xt, method = "guerrero")

fit <- forecast::auto.arima(Xt, lambda = lambda)
TS.diag(fit$res)
salida.arima(Xt, fit)

adf.test(resid(fit)) # No estacionaridad

acf(Xt) 
pacf(Xt)

Xt_new <- fitted(fit)



