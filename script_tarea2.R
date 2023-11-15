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
  filter(age_CE < 1890) %>% 
  select(age_CE,X540021_raw) %>% 
  rename("tiempo" = age_CE, "ancho" = X540021_raw)

str(datos)

# Volviendo a variables temporales ----
Xt <- ts(datos$ancho, start = 1721, frequency = 1)
t <- as.numeric(time(Xt, start = 1721))

# Gráfico de serie de tiempo ----
datos %>%
  ggplot(aes(x = tiempo, y = ancho)) +
  geom_line(color = "brown4", lwd = 0.8)+
  geom_hline(yintercept = mean(datos$ancho),
             col = "darkgreen",lwd = 0.8) +
  labs(x = "Año", y = "Ancho anillo [mm]", title = "Variación del ancho del anillo de pinos silvestres de Noruega en mm") +
  theme_bw()

# Box-Plot ----

box_plot = datos %>%
  ggplot(aes(x = ancho, y = 0)) +
  # geom_violin(color = "brown3") +
  geom_boxplot(fill = "brown4",
               width = 0.3,color = "darkgreen") +
  ylim(c(-0.5, 0.5)) +
  theme_bw() +
  xlab("CO ppm") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10,
                                  #face = "bold",
                                  color = "black",
                                  hjust = 0.5),
        plot.subtitle = element_text(size = 9,
                                     #face = "bold",
                                     color = "black",
                                     hjust = 0.5),
        axis.title.x = element_text(size = 8)) +
  labs(title = "CAMBIAR ESTA WEA",
       subtitle = "ESTE TMB")

box_plot

# ACF y PACF ----

acf(datos$X540021_raw, main = "ACF")
pacf(datos$X540021_raw, main = "PACF")

# ACF base
acf <- acf(Xt, plot = FALSE)
acf_data <- data.frame(Lag = acf$lag,
                       ACF = acf$acf)

# PACF base
pacf <- pacf(Xt, plot = FALSE)
pacf_data <- data.frame(Lag = pacf$lag,
                        PACF = pacf$acf)

n = nrow(datos)

ACF = ggplot(acf_data, aes(x = Lag, y = ACF)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") + 
  geom_hline(yintercept = 1.96/sqrt(n),
             linetype = "dashed",
             col = "darkgreen") + 
  geom_hline(yintercept = -1.96/sqrt(n),
             linetype = "dashed",
             col = "darkgreen") + 
  geom_segment(aes(xend = Lag,
                   yend = 0),
               color = "brown",
               size = 1) + # Líneas verticales
  geom_point(size = 2.5,
             shape = 17,
             col = "darkgreen") + # Puntos de color en el extremo de las líneas
  labs(x = "Lag",
       y = "ACF",
       title = "Gráfico de Autocorrelación de los datos") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        # axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 10,
                                  #face = "bold",
                                  color = "black",
                                  hjust = 0.5),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8))

PACF = pacf_data %>%
  ggplot(aes(x = Lag, y = PACF)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") + 
  geom_hline(yintercept = 1.96/sqrt(n),
             linetype = "dashed",
             col = "darkgreen") + 
  geom_hline(yintercept = -1.96/sqrt(n),
             linetype = "dashed",
             col = "darkgreen") + 
  geom_segment(aes(xend = Lag,
                   yend = 0),
               color = "brown",
               size = 1) + # Líneas verticales
  geom_point(size = 2.5,
             shape = 17,
             col = "darkgreen") + # Puntos de color en el extremo de las líneas
  labs(x = "Lag",
       y = "ACF",
       title = "Gráfico de Autocorrelación Parcial de los datos") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        # axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 10,
                                  #face = "bold",
                                  color = "black",
                                  hjust = 0.5),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8))


PACF
ACF

# ARIMA ----

auto.arima(Xt)
## AR(1)
fit <- forecast::Arima(Xt, order = c(0,1,1), include.mean = F)


# Tests ----

shapiro.test(fit$res) # son normales 
bptest(fit$res ~ t) # son homocedasticos

TS.diag(fit$res) # pasa test de blacura
salida.arima(Xt, fit) # significancia de los coeficientes

adf.test(fit$res) # SON ESTACIONARIOS

acf(Xt) 
pacf(Xt)

Xt_new <- fitted(fit)



