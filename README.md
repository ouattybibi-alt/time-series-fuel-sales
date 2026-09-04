# Time Series Analysis — French Automotive Fuel Sales

ARMA and SARIMA modelling of monthly gasoline and diesel sales in France, with in-sample diagnostics, out-of-sample validation, and 12-month forecasts.

Academic project (M2 Actuariat, ISFA), supervised by Prof. Christian Robert. Co-authored with C. A. D. Kouamé and S. Ouattara.

## What this does

Two independent scripts, each following the same modelling workflow on a different series:

- **STL decomposition** of the log-transformed series to separate trend, seasonality and residuals
- **Stationarity tests** on the residual component (Augmented Dickey–Fuller, KPSS)
- **ARMA(p, q) grid search** with joint selection on AIC, Ljung–Box and Box–Pierce p-values
- **Out-of-sample validation** on the 2020–onward period (MAE, RMSE)
- **SARIMA(p, d, q)(P, D, Q)[s]** as an alternative to detrending + deseasonalisation
- **12-month forecasts** with confidence intervals

## Files

| File | Content |
|---|---|
| `supercarburants.R` | Gasoline sales (~1981–present) — ARMA approach on STL residuals |
| `gazole.R` | Diesel sales — comparison of two approaches: ARMA on residuals vs SARIMA on the raw series |

## Requirements

R with the following packages:

```r
install.packages(c("readr", "forecast", "tseries", "ggplot2",
                   "lmtest", "FinTS", "dplyr", "lubridate"))
```

The scripts read `Data.csv` from a local path (hardcoded to the original author's machine). Replace the `read.csv(...)` path at the top of each script with your own location before running.

## Status

Reconstructed from the PDF appendix of the report — original R scripts not preserved as standalone files.

- Structural validation OK: all brackets, quotes and calls balance correctly.
- Four calls in the original code had commented-out arguments (`# lag.max = 120, ...`) that made the enclosing parenthesis invisible to R. Restored to functional form.
- Section titles were converted to R comment blocks so the file structure mirrors the report's appendix.

## Full write-up

Full methodology, plots and statistical results are in the project report (available on request).
