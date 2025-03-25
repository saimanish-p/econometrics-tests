# Loading required libraries
library(vars)
library(tseries)
library(haven)

# Reading the STATA data file
Model_3 <- read_dta("/Users/saimanishprabhakar/Desktop/Dissertation/Dissertation Data and Estimation/STATA Dissertation files/STATA data files/Model_3 (USEPU Endog & VIX Exog).dta")

# Examining the structure and first few rows of the data
str(Model_3)
head(Model_3)

# Converting the data to a time series object
Model_3_ts <- ts(Model_3, start = c(2004, 1), frequency = 12)  

# Checking for missing values 
colSums(is.na(Model_3))

# Removing rows with missing values
Model_3 <- na.omit(Model_3)

# Defining endogenous variables
endog_vars <- cbind(Model_3$dlogWTI, 
                    Model_3$dlogUSFPCO, 
                    Model_3$dlogUSINDPROD, 
                    Model_3$ddlogUSINV, 
                    Model_3$dlogDXY, 
                    Model_3$dlogUSEPU)

# Defining exogenous variables
exog_vars <- cbind(Model_3$dlogVIX, Model_3$L1_dlogVIX, Model_3$L2_dlogVIX)

# Assigning column names to the endogenous and exogenous variables (optional, but helpful)
colnames(endog_vars) <- c("dlogWTI", "dlogUSFPCO", "dlogUSINDPROD", "ddlogUSINV", "dlogDXY", "dlogUSEPU")
colnames(exog_vars) <- c("dlogVIX", "L1_dlogVIX", "L2_dlogVIX")

# Estimating the VAR model with 4 lags and exogenous variables
VAR_model <- VAR(endog_vars, p = 4, exogen = exog_vars)

# Extracting residuals from the VAR model
var_residuals <- residuals(VAR_model)

# Looping through all residuals and performing the BDS test for each variable
for (i in 1:ncol(var_residuals)) {
  bds_result <- bds.test(var_residuals[,i])
  cat("\nBDS Test for residuals of variable", colnames(var_residuals)[i], ":\n")
  print(bds_result)
}