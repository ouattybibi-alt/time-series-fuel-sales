# ==============================================================
# Modelisation ARMA/SARIMA des ventes de gazole automobile
# Source : annexe B du rapport (M2 Actuariat, ISFA)
# Extrait du PDF puis nettoye - a relire avant utilisation
# ==============================================================

# #########################################################
# SECTION B  -  Script R : Modélisation des ventes de Gazole
# #########################################################

# =========================================================
# B.1  Chargement des données
# =========================================================

library (readr)
library (dplyr)

df <- read.csv ("C:\COURS ET ANNALES ISFA\Series temporelles\Data.csv",
skip = 1,        # saute la 1 re ligne
header = TRUE,
sep = ";")

str (df)
head (df)
summary (df)

# =========================================================
# B.2  Récupération de la variable à prédire
# =========================================================

library (lubridate)

df <- df |>
mutate (PERIODE = ym (PERIODE)) |>
arrange (PERIODE)

y_gaz    <- ts (df$VGAZOL_RAFF,
start = c (1980, 1), frequency = 12)

# =========================================================
# B.3  Visualisation graphique et quelques statistiques
# =========================================================

plot (y_gaz, main = "Ventes de gazole auto", ylab = "kt")
boxplot (y_gaz~cycle (y_gaz), main = "Saisonnalite", xlab = "Mois", ylab = "kt")

summary (y_gaz)
tapply (as.numeric (y_gaz), cycle (y_gaz), mean, na.rm = TRUE)

# =========================================================
# B.4  Décomposition de la série
# =========================================================

y_gaz <- log (y_gaz) # modele multiplicatif donc on prend log

decomp_gaz <- decompose (y_gaz, type = "additive") # Donc decomposition additive
         comme vu plus haut
plot (decomp_gaz)

decomp_gaz$trend
sum (decomp_gaz$seasonal [1:12]) # Verifier de la condition d ' identifiabilite de la
          saisonnalite

mean (decomp_gaz$random, na.rm = TRUE) # Le bruit est cense etre centre en 0
sd (decomp_gaz$random, na.rm = TRUE)

# =========================================================
# B.5  Vérification de la stationnarité de la série
# =========================================================

y_gaz_clean <- y_gaz [! is.na (y_gaz)]
acf (y_gaz_clean, lag.max = 120,
main = "ACF des residus (sans NA) de y_gaz")

# =========================================================
# B.6  Méthode 1 : ajustement d’un ARMA(p,q) sur le bruit
# =========================================================

# =========================================================
# B.6.1  Tansformation de la série : retrait du trend et de la saisonnalité
# =========================================================

dy_gaz <-y_gaz - decomp_gaz$trend
plot (dy_gaz, main = "Difference premiere de y_gaz", ylab = "dy_gaz", xlab = "Temps"
)

# 3. ACF et PACF de la serie differenciee
dy_gaz_clean <- dy_gaz [! is.na (dy_gaz)]
acf (dy_gaz_clean, lag.max = 120, main = "ACF de dy_gaz")
pacf (dy_gaz_clean, lag.max = 120, main = "PACF de dy_gaz")

# =========================================================
# B.6.2  Vérification de la non stationnarité de la série
# =========================================================

library (tseries)
adf_dy_gaz <- adf.test (dy_gaz_clean, k = 120)
adf_dy_gaz

# =========================================================
# B.6.3  Traitement de la saisonnalité
# =========================================================

dsy_gaz <- dy_gaz - decomp_gaz$seasonal
# 2. Representation de la serie differenciee
plot (dsy_gaz, main = "dy_gaz desaisonnalisee", ylab = "dsy_gaz", xlab = "Temps")

# =========================================================
# B.6.4  ACF et PACF de la série détendancialisée et désaisonnalisée
# =========================================================

# 3. ACF de la serie desaisonnalisee
dsy_gaz_clean <- dsy_gaz [! is.na (dsy_gaz)]
acf (dsy_gaz_clean, lag.max = 120, main = "ACF de dsy_y_gaz")
pacf (dsy_gaz_clean, lag.max = 120, main = "PACF de dsy_gaz")

# =========================================================
# B.6.5  Test de stationnarité
# =========================================================

library (tseries)
adf_dsy_gaz <- adf.test (dsy_gaz_clean, k =60)
adf_dsy_gaz

library (tseries)
library (zoo)
# 1. Nettoyage des NA (interpolation)
y_gaz_pp_vec <- na.approx (as.numeric (dsy_gaz))

y_gaz_pp_ts <- ts (y_gaz_pp_vec, start      = start (dsy_gaz), frequency = frequency (
           dsy_gaz))

# 2. Test de Phillips - Perron
pp_y_gaz <- pp.test (y_gaz_pp_ts,
alternative = "stationary", # H1 : stationnaire
type = "Z (t_alpha)") # statistique disponible
pp_y_gaz

library (tseries)
library (zoo)

# 1. Nettoyer les NA
y_gaz_kpss_vec <- na.approx (as.numeric (dsy_gaz))
y_gaz_kpss_ts <- ts (y_gaz_kpss_vec, start = start (dsy_gaz), frequency = frequency (
           dsy_gaz))

# 2. Test KPSS avec stationnarite autour d ' une tendance
kpss_y_gaz_trend <- kpss.test (y_gaz_kpss_ts, null = "Trend")

kpss_y_gaz_trend

# =========================================================
# B.6.6  Vérification des autocorrélations des résidus
# =========================================================

Box.test (dsy_gaz_clean, lag = 120, type = "Box - Pierce")

# =========================================================
# B.6.7  Analyse automatique du meilleur ARMA avec le moins de paramètres à estimer
# =========================================================

library (forecast)
library (dplyr)
library (knitr)

x <- dsy_gaz_clean
x <- x [! is.na (x)]

p_max <- 5
q_max <- 10

results <- expand.grid (p = 0: p_max,
q = 0: q_max) |>
rowwise () |>
mutate (
fit = list (
tryCatch (
Arima (x, order = c (p, 0, q), seasonal = c (0, 0, 0)),
error = function (e) NULL
)
),

AIC = if (! is.null (fit)) AIC (fit) else NA_real _,
BIC = if (! is.null (fit)) BIC (fit) else NA_real _
) |>
ungroup () |>
filter (! is.na (AIC), ! is.na (BIC)) |>
select (- fit) |>
arrange (BIC)

# Les 5 meilleurs modeles selon le BIC (avec leur AIC)
top_models <- head (results, 5)

library (forecast)
library (dplyr)

# Serie stationnaire
x <- dsy_gaz_clean
x <- x [! is.na (x)]

# 1. Estimation des          3 modeles
arma_44 <- Arima (x,      order = c (4, 0, 4), seasonal = c (0, 0, 0))
arma_45 <- Arima (x,      order = c (4, 0, 5), seasonal = c (0, 0, 0))
arma_54 <- Arima (x,      order = c (5, 0, 4), seasonal = c (0, 0, 0))

# 2. Fonction d ' extraction RMSE / MAE / AIC / BIC
get_stats <- function (mod) {
acc <- accuracy (mod)
# La premiere ligne correspond au training set
rmse <- acc [1, "RMSE"]
mae <- acc [1, "MAE"]
aic <- AIC (mod)
bic <- BIC (mod)
c (RMSE = rmse, MAE = mae, AIC = aic, BIC = bic)
}

# 3. Construction du tibble
errors_tbl <- tibble (
modele = c ("ARMA (4,4)", "ARMA (4,5)", "ARMA (5,4)")
) |>
bind_cols (
as_tibble (
t(
sapply (list (arma_44, arma_45, arma_54), get_stats)
)
)
)

errors_tbl

# =========================================================
# B.6.8  Analyse du modèle final
# =========================================================

coeftest (arma_44)
arma_44 <- Arima (x, order = c (4, 0, 4), seasonal = c (0, 0, 0))
# Residus
res_m2 <- residuals (arma_44)

par (mfrow = c (2,2))

plot (res_m2, main = "Residus du modele ARIMA")
# acf (res_m2, lag.max = 48, main = "ACF des residus")
# pacf (res_m2, lag.max = 48, main = "PACF des residus")
qqnorm (res_m2) ; qqline (res_m2)
par (mfrow = c (1,1))

# Test de Ljung - Box : autocorrelation residuelle globale
Box.test (res_m2, lag = 120, type = "Ljung - Box", fitdf = length (coef (arma_44)))
Box.test (res_m2, lag = 120, type = "Box - Pierce", fitdf = length (coef (arma_44)))

# ou avec forecast :
checkresiduals (arma_44, lag =120)
mean (sarima$residuals)

# =========================================================
# B.6.9  Prévisions du modèle final
# =========================================================

#   ============================================================
#   PREVISIONS A 12 MOIS DES VENTES DE GAZOLE
#   en reprenant les variables du script actuel
#   ============================================================

library (forecast)

#   ----------------------------
#   1) Objets deja construits
#   ----------------------------
#   y_gaz           : serie log - transformee
#   decomp_gaz      : decomposition additive de y_gaz
#   dsy_gaz_clean : serie log, detrendee et desaisonnalisee
#   arma_44         : modele ARMA (4,4) retenu

# ----------------------------
# 2) Extraire tendance et saisonnalite
# ----------------------------
trend_gaz <- decomp_gaz$trend
season_gaz <- decomp_gaz$seasonal

# ----------------------------
# 3) Prevision de la composante ARMA a 12 mois
# ----------------------------
h <- 12
pred_arma <- predict (arma_44, n.ahead = h)

arma_mean <- as.numeric (pred_arma$pred)
arma_se   <- as.numeric (pred_arma$se)

# ----------------------------
# 4) Prevision de la tendance
#     Regression lineaire sur les 60 derniers mois non manquants
# ----------------------------
trend_gaz_clean <- trend_gaz [! is.na (trend_gaz)]

n_trend_window <- 60
n <- length (trend_gaz_clean)

trend_train <- trend_gaz_clean [(n - n_trend_window + 1) : n]
t_train <- 1: length (trend_train)

trend_lm <- lm (trend_train~t_train)

t_future <- (length (trend_train) + 1) :(length (trend_train) + h)
pred_trend <- predict (trend_lm,
newdata = data.frame (t_train = t_future),
se.fit = TRUE)

trend_mean <- as.numeric (pred_trend$fit)
trend_se   <- as.numeric (pred_trend$se.fit)

# ----------------------------
# 5) Prevision de la saisonnalite
#    On repete le profil mensuel moyen
# ----------------------------
season_pattern <- tapply (season_gaz, cycle (y_gaz), mean, na.rm = TRUE)

end_year <- end (y_gaz) [1]
end_month <- end (y_gaz) [2]

future_months <- ((end_month + (1: h) - 1) %% 12) + 1
season_mean <- as.numeric (season_pattern [future_months])

# ----------------------------
# 6) Reconstruction sur l ' echelle log
# ----------------------------
log_forecast_mean <- trend_mean + season_mean + arma_mean

# Approximation de l ' ecart - type total sur l ' echelle log
log_forecast_se <- sqrt (trend_se ^2 + arma_se ^2)

# Bornes 95%
log_lower95 <- log_forecast_mean - 1.96 * log_forecast_se
log_upper95 <- log_forecast_mean + 1.96 * log_forecast_se

# ----------------------------
# 7) Retour a l ' echelle initiale
#    avec correction du biais de retransformation
# ----------------------------
forecast_mean <- exp (log_forecast_mean + 0.5 * log_forecast_se ^2)
lower95 <- exp (log_lower95)
upper95 <- exp (log_upper95)

# ----------------------------
# 8) Tableau des previsions
# ----------------------------
future_dates <- seq (as.Date ("2025 -12 -01"), by = "month", length.out = h)

table_previsions_gaz <- data.frame (
Date = future_dates,
Forecast_kt = round (forecast_mean, 2),
Lower_95_kt = round (lower95, 2),
Upper_95_kt = round (upper95, 2)
)

print (table_previsions_gaz)

# ----------------------------
# 9) Objets ts pour le graphique

# ----------------------------
forecast_ts <- ts (forecast_mean, start = c (2025, 12), frequency = 12)
lower_ts    <- ts (lower95,        start = c (2025, 12), frequency = 12)
upper_ts    <- ts (upper95,        start = c (2025, 12), frequency = 12)

# ----------------------------
# 10) Graphique : serie observee + previsions
# ----------------------------
ylim_all <- range (c (exp (y_gaz), forecast_ts, lower_ts, upper_ts), na.rm = TRUE)

plot (exp (y_gaz),
xlim = c (start (exp (y_gaz)) [1], end (forecast_ts) [1] + (end (forecast_ts) [2] -1) /
                  12),
ylim = ylim_all,
main = "Ventes observees et previsions de gazole sur 12 mois",
ylab = "kt",
xlab = "Temps",
col = "black",
lwd = 1.5)

# relier la derniere observation a la premiere prevision
last_obs_time <- time (exp (y_gaz)) [length (exp (y_gaz))]
last_obs_value <- tail (exp (y_gaz), 1)

first_fc_time <- time (forecast_ts) [1]
first_fc_value <- forecast_ts [1]

segments (x0 = last_obs_time, y0 = last_obs_value,
x1 = first_fc_time, y1 = first_fc_value,
col = "red", lwd = 2)

lines (forecast_ts, col = "red", lwd = 2)
lines (lower_ts, col = "blue", lty = 2, lwd = 1.5)
lines (upper_ts, col = "blue", lty = 2, lwd = 1.5)

legend ("topleft",
legend = c ("Observee", "Prevision", "IC 95%"),
col = c ("black", "red", "blue"),
lty = c (1, 1, 2),
lwd = c (1.5, 2, 1.5),
bty = "n")

#   ============================================================
#   GRAPHIQUE ZOOME A PARTIR DE 2017
#   serie observee + previsions + IC 95%
#   ============================================================

# fenetre observee a partir de 2017
obs_2017 <- window (exp (y_gaz), start = c (2017, 1))

# bornes de l ' axe y
ylim_all <- range (c (obs_2017, forecast_ts, lower_ts, upper_ts), na.rm = TRUE)

plot (obs_2017,
xlim = c (2017, end (forecast_ts) [1] + (end (forecast_ts) [2] -1) / 12),
ylim = ylim_all,
main = "Ventes observees et previsions de gazole a partir de 2017",
ylab = "kt",
xlab = "Temps",

col = "black",
lwd = 1.5)

# relier la derniere observation a la premiere prevision
last_obs_time <- time (exp (y_gaz)) [length (exp (y_gaz))]
last_obs_value <- tail (exp (y_gaz), 1)

first_fc_time <- time (forecast_ts) [1]
first_fc_value <- forecast_ts [1]

segments (x0 = last_obs_time, y0 = last_obs_value,
x1 = first_fc_time, y1 = first_fc_value,
col = "red", lwd = 2)

# tracer les previsions et les IC
lines (forecast_ts, col = "red", lwd = 2)
lines (lower_ts, col = "blue", lty = 2, lwd = 1.5)
lines (upper_ts, col = "blue", lty = 2, lwd = 1.5)

# ligne verticale au debut des previsions
abline (v = first_fc_time, col = "darkgray", lty = 3)

legend ("bottomleft",
legend = c ("Observee", "Prevision", "IC 95%"),
col = c ("black", "red", "blue"),
lty = c (1, 1, 2),
lwd = c (1.5, 2, 1.5),
bty = "n")

# =========================================================
# B.7  Méthode 2 : Ajustement d’un SARIMA(p,d,q)(P,D,Q)[s] sur la série
# =========================================================

# =========================================================
# B.7.1  Transformation de la série : différenciation ordinaire
# =========================================================

# 1. Difference d ' ordre 1 (detendancialisation)
dy_gaz <- diff (y_gaz, differences = 1)
dy_gaz

# 2. Representation de la serie differenciee
plot (dy_gaz,
main = "Difference premiere de y_gaz",
ylab = "dy_gaz", xlab = "Temps")

# 3. ACF de la serie differenciee
dy_gaz_clean <- dy_gaz [! is.na (dy_gaz)]

acf (dy_gaz_clean,
lag.max = 120,
main = "ACF de dy_gaz")

# =========================================================
# B.7.2  Vérification de la non stationnarité de la série différenciée(test ADF)
# =========================================================

library (tseries)
adf_dy_gaz <- adf.test (dy_gaz_clean, k = 120)
adf_dy_gaz

# =========================================================
# B.7.3  Transformation de la série : différenciation saisonnière
# =========================================================

dsy_gaz <- diff (dy_gaz, lag = 12, differences = 1)
dsy_gaz

# 2. Representation de la serie differenciee
plot (dsy_gaz,
main = "dy_gaz desaisonnalisee",
ylab = "dsy_gaz", xlab = "Temps")

# 3. ACF de la serie desaisonnalisee
dsy_gaz_clean <- dsy_gaz [! is.na (dsy_gaz)]

acf (dsy_gaz_clean,
lag.max = 120,
main = "ACF de dsy_gaz")

pacf (dsy_gaz_clean,
lag.max = 120,
main = "PACF de dsy_gaz")

# =========================================================
# B.7.4  Tests de stationnarité
# =========================================================

adf_dsy_gaz <- adf.test (dsy_gaz_clean, k =60)
adf_dsy_gaz

library (zoo)

# 1. Nettoyage des NA (interpolation)
y_gaz_pp_vec <- na.approx (as.numeric (dsy_gaz))

y_gaz_pp_ts <- ts (y_gaz_pp_vec,
start            = start (dsy_gaz),
frequency = frequency (dsy_gaz))

# 2. Test de Phillips - Perron
pp_y_gaz <- pp.test (y_gaz_pp_ts,
alternative = "stationary",             # H1 : stationnaire
type = "Z (t_alpha)")               # statistique disponible
pp_y_gaz

library (zoo)

# # 1. Nettoyer les NA (comme pour ADF / PP)
y_gaz_kpss_vec <- na.approx (as.numeric (dsy_gaz))

y_gaz_kpss_ts <- ts (y_gaz_kpss_vec,
start          = start (dsy_gaz),
frequency = frequency (dsy_gaz))

# # 2. Test KPSS avec stationnarite autour d ' une tendance
kpss_y_gaz_trend <- kpss.test (y_gaz_kpss_ts,
null = "Trend")      # H0 : stationnaire autour d ' une
                                                  tendance
kpss_y_gaz_trend

# # 3. (Optionnel) KPSS avec stationnarite autour d ' une constante

kpss_y_gaz_level <- kpss.test (y_gaz_kpss_ts,
null = "Level")      # H0 : stationnaire autour d ' une
                                                  moyenne constante
kpss_y_gaz_level

# =========================================================
# B.7.5  Recherche du modèle initial
# =========================================================

sarima <- Arima (y_gaz_clean, order = c (2, 1, 4), seasonal = list (order = c (2, 1, 2),
         period =12)) # modele d t e r m i n e par analyse de l ' ACF et du PACF
coeftest (sarima)

sarima <- Arima (y_gaz_clean, order = c (2, 1, 4), seasonal = list (order = c (1, 1, 2),
         period =12)) # modele ajuste conformement           l ' ancien coeftest
coeftest (sarima)

# =========================================================
# B.7.6  Recherche d’un modèle éfficace comportant moins de paramètres à estimer
# =========================================================

library (forecast)
library (dplyr)
library (knitr)

x <- y_gaz_clean
x <- x [! is.na (x)]

p_max   <-   2
q_max   <-   4
P_max   <-   1
Q_max   <-   2

results <- expand.grid (p = 0: p_max,
q = 0: q_max,
P = 0: P_max,
Q = 0: Q_max) |>
rowwise () |>
mutate (
fit = list (
tryCatch (
Arima (x, order = c (p, 1, q), seasonal = list (order = c (P, 1, Q), period =
                  12)),
error = function (e) NULL
)
),
AIC = if (! is.null (fit)) AIC (fit) else NA_real _,
BIC = if (! is.null (fit)) BIC (fit) else NA_real _,

LjungBox_pvalue = if (! is.null (fit)) tryCatch (
Box.test (residuals (fit), lag = 120, type = "Ljung - Box")$p.value,
error = function (e) NA_real _
) else NA_real _,

BoxPierce_pvalue = if (! is.null (fit)) tryCatch (
Box.test (residuals (fit), lag = 120, type = "Box - Pierce")$p.value,
error = function (e) NA_real _
) else NA_real _
) |>
ungroup () |>

filter (! is.na (AIC), ! is.na (BIC)) |>
select (- fit) |>
arrange (BIC)

# Les 5 meilleurs modeles selon le BIC (avec leur AIC)
top_models <- head (results, 5)
results

results_BB <- subset (results,
LjungBox_pvalue > 0.05 &
BoxPierce_pvalue > 0.05)
results_BB

library (forecast)
library (dplyr)

# Serie stationnaire
x <- y_gaz_clean
x <- x [! is.na (x)]

# 1. Estimation des 4 modeles
arima_2401 <- Arima (x, order =           c (2, 1, 4), seasonal = list (order = c (0, 1, 1),
         period =12))
arima_2411 <- Arima (x, order =           c (2, 1, 4), seasonal = list (order = c (1, 1, 1),
         period =12))
arima_2402 <- Arima (x, order =           c (2, 1, 4), seasonal = list (order = c (0, 1, 2), period
         =12))
arima_2412 <- Arima (x, order =           c (2, 1, 4), seasonal = list (order = c (1, 1, 2), period
         =12))

# 2. Fonction d ' extraction RMSE / MAE / AIC / BIC
get_stats <- function (mod) {
acc <- accuracy (mod)
# La premiere ligne correspond au training set
rmse <- acc [1, "RMSE"]
mae <- acc [1, "MAE"]
aic <- AIC (mod)
bic <- BIC (mod)
c (RMSE = rmse, MAE = mae, AIC = aic, BIC = bic)
}

# 3. Construction du tibble
errors_tbl <- tibble (
modele = c ("ARIMA (2,1,4) (0,1,1)", "ARIMA (2,1,4) (1,1,1)", "ARIMA (2,1,4) (0,1,2)",
            "ARIMA (2,1,4) (1,1,2)")
) |>
bind_cols (
as_tibble (
t(
sapply (list (arima_2401, arima_2411, arima_2402, arima_2412), get_stats)
)
)
)

errors_tbl

# =========================================================
# B.7.7  Analyse du modèle retenu
# =========================================================

sarima <- Arima (y_gaz_clean, order = c (2, 1, 4), seasonal = list (order = c (0, 1, 1),
         period =12))

coeftest (sarima)

# Residus
res_m2 <- residuals (sarima)

par (mfrow = c (2,2))
plot (res_m2, main = "Residus du modele ARIMA")
qqnorm (res_m2) ; qqline (res_m2)
par (mfrow = c (1,1))

# Test de Ljung - Box : autocorrelation residuelle globale
Box.test (res_m2, lag = 120, type = "Ljung - Box", fitdf = length (coef (sarima)))
Box.test (res_m2, lag = 120, type = "Box - Pierce", fitdf = length (coef (sarima)))

checkresiduals (sarima, lag =120)

mean (sarima$residuals)

# =========================================================
# B.7.8  Prévisions du modèle retenu
# =========================================================

library (forecast)
prevision = forecast (sarima, h =12)
plot (prevision)

previsions_kt <- exp (prevision$mean)

# DF final
df_final <- data.frame (
mois_ahead = 1:12,
prev_kt = round (previsions_kt, 0)
)

print (df_final)

