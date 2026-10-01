# FACAViaR-MIDAS-JBF
These are sample codes to estimate FACAViaR-MIDAS for Huang and Liu (2026) published in Journal of Banking and Finance: 
"Real-Time Macroeconomic Quantile Factors for Daily Value-at-Risk and Expected Shortfall Forecasting: Statistical Significance and Economic Value" (with Xiao Huang), Journal of Banking and Finance, 2026, 192: 107831

1. Code for estimating two quantile factors using CDG method
   a. FRED_MD_Data_Clean.R: this is a R file used to clean FRED-MD data "2026-02.csv". Cleaned data are saved as "FRED_MD_Clean"
   b. CDG_QFactors_Monthly.m: This is Matlab file to estimate 2 quantile factors using CDG method and "FRED_MD_Clean" data.
   c. ConvertCSVdataToRdata_Monthly.R: This is a R file to convert the estimated quantile factor csv files into a merged R data file, named as "CDG_QFactors2_Monthly.Rdata"
2. Code for estimating FACAViaR-MIDAS with 2 CDG quantile factors
    a. CAViaRMIDAS_CDG_2Factor.R: this is the R file to estimate the FACAViaR-MIDAS with 2 CDF quantile factors
    b. DailyData.Rdata: this is a Rdata file which provides daily financial returns.
    c. CAViaR_MIDAS_Joint_VaR_ES_Rcpp.R: This is a R file which provides R help functions.
    d. CAViaR_MIDAS_cpp.cpp: this is a C++ file which provides the core estimation functions.
