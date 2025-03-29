# MarketPulse-Analytics
Time Series Analysis of Global Stock Indices: Uncovering Market Dynamics. 

My methodology combines:  
Stationarity testing (ADF, KPSS, PP)  ,
DCC-GARCH for dynamic correlations, 
Granger causality for volatility transmission , 
Risk metrics (VaR, ES, Sharpe Ratio)

Visualizations including correlation plots and density distributions make complex relationships accessible, while a comprehensive CSV report summarizes key metrics. This end-to-end analysis demonstrates how advanced time series techniques can extract actionable intelligence from market data, supporting better investment decisions and risk management.

The project showcases professional-grade financial analytics using R's rugarch, vars, and ggplot2 packages, providing a template for analyzing any set of financial time series.


######################  CODE ########################################################################################################################


=======================================================================================================================================================
=======================================================================================================================================================


# Install and load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(xts, zoo, tseries, urca, rugarch, rmgarch, vars, 
               PerformanceAnalytics, quantmod, ggplot2, lubridate, 
               tidyr, dygraphs, riskmetrics)

# 1. Data Loading and Cleaning --------------------------------------------
# Read the data with proper NA handling
data <- read.table("shareprice.txt", header = TRUE, stringsAsFactors = FALSE, 
                   na.strings = c("NA", "N/A", ""))

# Convert Date column with multiple format handling
data$Date <- parse_date_time(data$Date, orders = c("mdy", "ymd", "dmy"))

# Remove rows with NA dates or prices
data <- data[complete.cases(data), ]

# Verify data structure
if (ncol(data) < 4) stop("Data should have at least 4 columns (Date + 3 indices)")
if (nrow(data) < 100) warning("Small sample size may affect analysis results")

# Convert to xts object
prices <- xts(data[, -1], order.by = data$Date)

# Handle zeros/negative prices that would cause problems with log returns
prices[prices <= 0] <- NA
prices <- na.omit(prices)
prices

# 2. Returns Calculation --------------------------------------------------
calculate_returns <- function(prices) {
  returns <- diff(log(prices))[-1, ]
  colnames(returns) <- c("N225", "BVSP", "SSCI")
  returns <- returns[complete.cases(returns), ]
  return(returns)
}

returns <- calculate_returns(prices)
returns

# 3. Stationarity Testing -------------------------------------------------
run_stationarity_tests <- function(returns) {
  tests <- list()
  for (i in 1:ncol(returns)) {
    series <- returns[, i]
    tests[[colnames(returns)[i]]] <- list(
      ADF = adf.test(na.omit(series)),
      KPSS = kpss.test(na.omit(series)),
      PP = pp.test(na.omit(series))
    )
  }
  return(tests)
}

stationarity_results <- run_stationarity_tests(returns)
stationarity_results

# 4. DCC-GARCH Modeling --------------------------------------------------
fit_dcc_garch <- function(returns) {
  # Univariate GARCH specifications
  uspec <- ugarchspec(
    mean.model = list(armaOrder = c(0, 0)),
    variance.model = list(garchOrder = c(1, 1)),
    distribution.model = "norm"
  )
  
  # Multivariate specification
  mspec <- multispec(replicate(ncol(returns), uspec))
  
  # DCC specification
  dccspec <- dccspec(
    uspec = mspec,
    dccOrder = c(1, 1),
    distribution = "mvnorm"
  )
  
  # Fit the model with error handling
  tryCatch({
    dccfit <- dccfit(dccspec, data = returns)
    return(dccfit)
  }, error = function(e) {
    message("DCC-GARCH fitting failed: ", e$message)
    return(NULL)
  })
}

dcc_fit <- fit_dcc_garch(returns)

if (!is.null(dcc_fit)) {
  # Plot conditional correlations
  plot(dcc_fit, which = 4)
  
  # Extract and plot dynamic correlations
  corr <- rcor(dcc_fit)
  time_correlations <- xts(corr[1, 2, ], order.by = index(returns))
  
  print(
    ggplot(data.frame(Date = index(time_correlations), 
                      Correlation = coredata(time_correlations)),
           aes(x = Date, y = Correlation)) +
      geom_line() +
      labs(title = "Dynamic Conditional Correlations", y = "Correlation") +
      theme_minimal()
  )
}

# 5. Volatility Spillover Analysis ----------------------------------------
analyze_spillovers <- function(returns) {
  # Determine optimal lag length
  lag_selection <- VARselect(returns, lag.max = 10)
  optimal_lag <- lag_selection$selection["AIC(n)"]
  
  # Estimate VAR model
  var_model <- VAR(returns, p = optimal_lag)
  
  # Calculate generalized FEVD (alternative to spilloverDY12)
  tryCatch({
    # Using fevd from vars package as alternative
    fevd_results <- fevd(var_model, n.ahead = 10)
    
    # Convert to spillover table format
    spillover_table <- matrix(0, ncol(returns), ncol(returns))
    colnames(spillover_table) <- rownames(spillover_table) <- colnames(returns)
    
    for (i in 1:length(fevd_results)) {
    spillover_table[i, ] <- colMeans(fevd_results[[i]]) * 100
  }
  
  # Calculate directional spillovers
  to_spillover <- colSums(spillover_table) - diag(spillover_table)
  from_spillover <- rowSums(spillover_table) - diag(spillover_table)
  net_spillover <- to_spillover - from_spillover
  
  # Create spillover index
  total_spillover <- mean(to_spillover)
  
  # Print results
  cat("\nSpillover Table (% contributions to variance):\n")
  print(round(spillover_table, 2))
  
  cat("\nDirectional Spillovers:\n")
  print(data.frame(
    To = to_spillover,
    From = from_spillover,
    Net = net_spillover
  ))
  
  cat("\nTotal Spillover Index:", round(total_spillover, 2), "%\n")
  
  return(list(
    spillover_table = spillover_table,
    to_spillover = to_spillover,
    from_spillover = from_spillover,
    net_spillover = net_spillover,
    total_spillover = total_spillover
  ))
  }, error = function(e) {
    message("Spillover analysis failed: ", e$message)
    return(NULL)
  })
}

spillover_results <- analyze_spillovers(returns)

# 6. Risk Metrics --------------------------------------------------------
calculate_risk_stats <- function(returns) {
  metrics <- list()
  for (i in 1:ncol(returns)) {
    series <- returns[, i]
    metrics[[colnames(returns)[i]]] <- list(
      VaR_95 = VaR(series, p = 0.95, method = "historical"),
      ES_95 = ES(series, p = 0.95, method = "historical"),
      Sharpe = SharpeRatio(series, Rf = 0, FUN = "StdDev")
    )
  }

# Portfolio metrics (equal-weighted)
portfolio_returns <- xts(rowMeans(returns), order.by = index(returns))
metrics$Portfolio <- list(
  VaR_95 = VaR(portfolio_returns, p = 0.95, method = "historical"),
  ES_95 = ES(portfolio_returns, p = 0.95, method = "historical"),
  Sharpe = SharpeRatio(portfolio_returns, Rf = 0, FUN = "StdDev")
)

return(metrics)
}

risk_metrics <- calculate_risk_stats(returns)
risk_metrics
# 7. Print All Results ---------------------------------------------------
cat("\n\n=== STATIONARITY TEST RESULTS ===\n")
print(stationarity_results)

cat("\n\n=== RISK METRICS ===\n")
print(risk_metrics)

if (!is.null(spillover_results)) {
  cat("\n\n=== SPILLOVER ANALYSIS RESULTS ===\n")
  print(spillover_results)
}

if (!is.null(dcc_fit)) {
  cat("\n\n=== DCC-GARCH MODEL SUMMARY ===\n")
  print(dcc_fit)
}
# 8. Comparative Analysis and Investment Recommendation ------------------

generate_investment_recommendation <- function(returns, risk_metrics, spillover_results) {
  
  # Create comparison dataframe
  comparison_df <- data.frame(
    Index = colnames(returns),
    Avg_Return = colMeans(returns, na.rm = TRUE),
    Volatility = apply(returns, 2, sd, na.rm = TRUE),
    VaR_95 = sapply(risk_metrics[1:3], function(x) x$VaR_95),
    ES_95 = sapply(risk_metrics[1:3], function(x) x$ES_95),
    Sharpe = sapply(risk_metrics[1:3], function(x) x$Sharpe),
    Spillover_To = if(!is.null(spillover_results)) spillover_results$to_spillover else NA,
    Spillover_From = if(!is.null(spillover_results)) spillover_results$from_spillover else NA
  )
  
  # Calculate composite score (higher is better)
  comparison_df$Composite_Score <- with(comparison_df, {
    # Normalize each metric (all converted to "higher is better" logic)
    norm_return <- (Avg_Return - min(Avg_Return)) / (max(Avg_Return) - min(Avg_Return))
    norm_sharpe <- (Sharpe - min(Sharpe)) / (max(Sharpe) - min(Sharpe))
    norm_volatility <- 1 - ((Volatility - min(Volatility)) / (max(Volatility) - min(Volatility)))
    norm_var <- 1 - ((VaR_95 - min(VaR_95)) / (max(VaR_95) - min(VaR_95)))
    norm_es <- 1 - ((ES_95 - min(ES_95)) / (max(ES_95) - min(ES_95)))
    norm_spillover <- if(!is.null(spillover_results)) {
      (Spillover_To - min(Spillover_To)) / (max(Spillover_To) - min(Spillover_To))
    } else 0
    
    # Weighted sum (adjust weights as needed)
    0.25*norm_return + 0.30*norm_sharpe + 0.15*norm_volatility + 
      0.10*norm_var + 0.10*norm_es + 0.10*norm_spillover
  })
  
  # Rank the indices
  comparison_df$Rank <- rank(-comparison_df$Composite_Score)
  comparison_df$Rank
  comparison_df
  # Generate recommendations
  best_index <- comparison_df[which.max(comparison_df$Composite_Score), "Index"]
  worst_index <- comparison_df[which.min(comparison_df$Composite_Score), "Index"]
  
  cat("\n\n=== INVESTMENT RECOMMENDATION ===\n")
  cat("Based on comprehensive analysis of returns, risk, and spillover effects:\n")
  cat("\nTOP PERFORMING INDEX:", best_index, "\n")
  cat("REASONS:\n")
  cat("- Highest risk-adjusted returns (Sharpe Ratio:", 
      round(comparison_df[comparison_df$Index == best_index, "Sharpe"], 3), ")\n")
  cat("- Favorable return-to-risk profile\n")
  if(!is.null(spillover_results)) {
    cat("- Contributes", round(comparison_df[comparison_df$Index == best_index, "Spillover_To"], 1),
        "% to market volatility while receiving only", 
        round(comparison_df[comparison_df$Index == best_index, "Spillover_From"], 1), "%\n")
  }
  
  cat("\nLEAST FAVORABLE INDEX:", worst_index, "\n")
  cat("REASONS:\n")
  cat("- Lower risk-adjusted returns (Sharpe Ratio:", 
      round(comparison_df[comparison_df$Index == worst_index, "Sharpe"], 3), ")\n")
  cat("- Higher downside risk (ES 95%:", 
      round(comparison_df[comparison_df$Index == worst_index, "ES_95"], 3), ")\n")
  
  cat("\nDETAILED COMPARISON:\n")
  print(comparison_df[order(comparison_df$Rank), ])
  
  cat("\nINVESTMENT STRATEGY SUGGESTIONS:\n")
  cat("1. Consider overweight allocation to", best_index, "in your portfolio\n")
  cat("2. Use", best_index, "as a core holding for stable returns\n")
  cat("3. If holding", worst_index, ", consider hedging strategies\n")
  cat("4. Monitor spillover effects from other markets regularly\n")
  
  # Visualization
  print(
    ggplot(comparison_df, aes(x = reorder(Index, Composite_Score), y = Composite_Score, fill = Index)) +
      geom_bar(stat = "identity") +
      labs(title = "Composite Performance Score by Index", 
           x = "Index", y = "Composite Score") +
      theme_minimal() +
      theme(legend.position = "none")
  )
  
  return(comparison_df)
}

# Generate the final investment recommendation
investment_decision <- generate_investment_recommendation(returns, risk_metrics, spillover_results)
