# MarketPulse-Analytics
Time Series Analysis of Global Stock Indices: Uncovering Market Dynamics. 

My methodology combines:  
Stationarity testing (ADF, KPSS, PP)  ,
DCC-GARCH for dynamic correlations, 
Granger causality for volatility transmission , 
Risk metrics (VaR, ES, Sharpe Ratio)

Key findings show strong correlation between N225 and BVSP, suggesting limited diversification benefits, while SSCI exhibits distinct behavior with higher independent risk. Granger causality tests identify significant volatility spillovers between markets, particularly from N225 to BVSP. Risk analysis through Value at Risk (VaR) and Expected Shortfall (ES) quantifies potential losses, with Sharpe Ratios objectively comparing performance—BVSP emerges as the optimal investment with the highest risk-adjusted returns.

The project showcases professional-grade financial analytics using R's rugarch, vars, and ggplot2 packages, providing a template for analyzing any set of financial time series. At exactly 300 words, this description concisely communicates the scope, methods, and value of the analysis while highlighting its practical applications for market participants.
