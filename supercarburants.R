# ==============================================================
# Modelisation ARMA des ventes de supercarburants automobiles
# Source : annexe A du rapport (M2 Actuariat, ISFA)
# Extrait du PDF puis nettoye - a relire avant utilisation
# ==============================================================

# #########################################################
# SECTION A  -  Script R : Modélisation des ventes de supercarburants
# #########################################################

# =========================================================
# A.1  Chargement des bibliothèques
# =========================================================

library (readr)
library (forecast)
library (tseries)
library (ggplot2)
library (AST)
library (lmtest)
library (FinTS)

# =========================================================
# A.2  Chargement et préparation des données
# =========================================================

data <- read.csv ("C:/Users/dylan/Documents/Disque local/ISFA/Data.csv",
sep = ";",
header = TRUE,
check.names = FALSE)

data <- data [-1,]
data <- data [nrow (data) :1,]
data <- data [-(1:12),]
data [, -1] <- lapply (data [, -1], as.numeric)

# =========================================================
# A.3  Création des séries temporelles
# =========================================================

ts_super <- ts (data [["5.2 Ventes de supercarburants auto (en kt)"]],
start = c (1981, 1),
frequency = 12)

log_super <- log (ts_super)

# =========================================================
# A.4  Décomposition STL et extraction de la composante résiduelle
# =========================================================

stl_super <- stl (log_super, s.window = "periodic")
super_deseason_detrend <- stl_super$time.series [, "remainder"]

# =========================================================
# A.5  Tests de stationnarité
# =========================================================

adf.test (na.omit (super_deseason_detrend))
kpss.test (na.omit (super_deseason_detrend))

# =========================================================
# A.6  Estimation initiale du modèle ARMA
# =========================================================

modele_arma <- arima (super_deseason_detrend, order = c (5,0,3))
modele_arma

# =========================================================
# A.7  Comparaison des modèles ARMA candidats
# =========================================================

results <- data.frame ()

for (p in 0:5) {
for (q in 0:3) {

model <- try (arima (super_deseason_detrend,
order = c (p,0, q)),
silent = TRUE)

if (! inherits (model, "try - error")) {

lb <- Box.test (residuals (model),
lag = 108,
type = "Ljung - Box",
fitdf = p + q)

bp <- Box.test (residuals (model),
lag = 108,
type = "Box - Pierce",
fitdf = p + q)

results <- rbind (results,
data.frame (
p = p,
q = q,
AIC = AIC (model),
LjungBox_pvalue = lb$p.value,
BoxPierce_pvalue = bp$p.value
))
}
}
}

results

# =========================================================
# A.8  Sélection des modèles admissibles
# =========================================================

subset (results,
LjungBox_pvalue > 0.05 &
BoxPierce_pvalue > 0.05)

# =========================================================
# A.9  Validation hors échantillon
# =========================================================

train <- window (super_deseason_detrend, end = c (2019,12))
test <- window (super_deseason_detrend, start = c (2020,1))

modele_arma <- arima (train, order = c (5,0,3))

pred <- forecast (modele_arma, h = length (test))
errors <- test - pred$mean

MAE <- mean (abs (errors), na.rm = TRUE)
RMSE <- sqrt (mean (errors ^2, na.rm = TRUE))

MAE
RMSE

# =========================================================
# A.10  Diagnostic global des résidus
# =========================================================
modele_arma <- arima (super_deseason_detrend, order = c (5,0,3))
checkresiduals (modele_arma)

# =========================================================
# A.11  Racines du modèle
# =========================================================
roots_ar <- polyroot (c (1, - modele_arma$model$phi))
roots_ma <- polyroot (c (1, modele_arma$model$theta))

Mod (roots_ar)
Mod (roots_ma)

# =========================================================
# A.12  Reconstruction de la série initiale
# =========================================================
trend_super <- stl_super$time.series [, "trend"]
season_super <- stl_super$time.series [, "seasonal"]

modele_arma_super <- arima (super_deseason_detrend, order = c (5,0,3))

arma_fitted_super <- super_deseason_detrend - residuals (modele_arma_super)

log_super_hat <- trend_super + season_super + arma_fitted_super
super_hat <- exp (log_super_hat)

# =========================================================
# A.13  Indicateurs d’ajustement
# =========================================================
obs <- ts_super
pred <- super_hat

valid <- complete.cases (obs, pred)
obs <- obs [valid]
pred <- pred [valid]

errors <- obs - pred

MAE    <-   mean (abs (errors))
RMSE   <-   sqrt (mean (errors ^2))
MAPE   <-   mean (abs (errors / obs)) * 100
R2     <-   1 - sum (errors ^2) / sum ((obs - mean (obs)) ^2)
BIAS   <-   mean (errors)

results <- data.frame (
MAE = MAE,
RMSE = RMSE,
MAPE_percent = MAPE,
R2 = R2,
Bias = BIAS
)

results

# =========================================================
# A.14  Prévisions à 12 mois
# =========================================================

trend_super <- stl_super$time.series [, "trend"]
season_super <- stl_super$time.series [, "seasonal"]
remainder_super <- stl_super$time.series [, "remainder"]

modele_arma_super <- arima (remainder_super, order = c (5,0,3))

h <- 12
pred_arma <- predict (modele_arma_super, n.ahead = h)

arma_mean <- as.numeric (pred_arma$pred)
arma_se   <- as.numeric (pred_arma$se)

n_trend_window <- 60
n <- length (trend_super)

trend_train <- trend_super [(n - n_trend_window + 1) : n]
t_train <- 1: length (trend_train)

trend_lm <- lm (trend_train~t_train)

t_future <- (length (trend_train) + 1) :(length (trend_train) + h)
pred_trend <- predict (trend_lm,
newdata = data.frame (t_train = t_future),
se.fit = TRUE)

trend_mean <- as.numeric (pred_trend$fit)
trend_se   <- as.numeric (pred_trend$se.fit)

season_pattern <- tapply (season_super, cycle (season_super), mean)

end_year <- end (ts_super) [1]
end_month <- end (ts_super) [2]

future_months <- ((end_month + (1: h) - 1) %% 12) + 1
season_mean <- as.numeric (season_pattern [future_months])

log_forecast_mean <- trend_mean + season_mean + arma_mean
log_forecast_se <- sqrt (trend_se ^2 + arma_se ^2)

log_lower95 <- log_forecast_mean - 1.96 * log_forecast_se
log_upper95 <- log_forecast_mean + 1.96 * log_forecast_se

forecast_mean <- exp (log_forecast_mean + 0.5 * log_forecast_se ^2)
lower95 <- exp (log_lower95)
upper95 <- exp (log_upper95)

future_dates <- seq (as.Date ("2025 -12 -01"), by = "month", length.out = h)

table_previsions <- data.frame (
Date = future_dates,
Forecast_kt = round (forecast_mean, 2),
Lower_95_kt = round (lower95, 2),
Upper_95_kt = round (upper95, 2)
)

table_previsions

