# Loading required libraries
library(haven)
library(vars)
library(strucchange)

# Reading the STATA data file
Model_3 <- read_dta("/Users/saimanishprabhakar/Desktop/Dissertation/Dissertation Data and Estimation/STATA Dissertation files/STATA data files/Model_3 (USEPU Endog & VIX Exog).dta")

# Checking the structure and first few rows of the data
str(Model_3)
head(Model_3)

# Converting the data to a time series object
Model_3_ts <- ts(Model_3, start = c(2004, 1), frequency = 12)  

# Checking for missing values in each column
colSums(is.na(Model_3))

# Removing rows with missing values
Model_3 <- na.omit(Model_3)

# Defining endogenous variables for the VAR model
endog_vars <- cbind(Model_3$dlogWTI, 
                    Model_3$dlogUSFPCO, 
                    Model_3$dlogUSINDPROD, 
                    Model_3$ddlogUSINV, 
                    Model_3$dlogDXY, 
                    Model_3$dlogUSEPU)

# Defining exogenous variables for the VAR model
exog_vars <- cbind(Model_3$dlogVIX, Model_3$L1_dlogVIX, Model_3$L2_dlogVIX)

# Assigning column names to endogenous and exogenous variables
colnames(endog_vars) <- c("dlogWTI", "dlogUSFPCO", "dlogUSINDPROD", "ddlogUSINV", "dlogDXY", "dlogUSEPU")
colnames(exog_vars) <- c("dlogVIX", "L1_dlogVIX", "L2_dlogVIX")

# Estimating the VAR model with 4 lags and exogenous variables
VAR_model <- VAR(endog_vars, p = 4, exogen = exog_vars)

# Extracting residuals from the VAR model
var_residuals <- residuals(VAR_model)

# Performing CUSUM and CUSUM of Squares tests for all variables in the VAR model

# For dlogWTI
cusum_test_var1 <- efp(var_residuals[,1] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var1, main = "Rec-CUSUM for dlogWTI", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var1 <- efp(var_residuals[,1] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var1, main = "CUSUM of Squares for dlogWTI", xlab = "Observation", ylab = "Cumulative Sum of Squares")

# For dlogUSFPCO
cusum_test_var2 <- efp(var_residuals[,2] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var2, main = "Rec-CUSUM for dlogUSFPCO", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var2 <- efp(var_residuals[,2] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var2, main = "CUSUM of Squares for dlogUSFPCO", xlab = "Observation", ylab = "Cumulative Sum of Squares")

# For dlogUSINDPROD
cusum_test_var3 <- efp(var_residuals[,3] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var3, main = "Rec-CUSUM for dlogUSINDPROD", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var3 <- efp(var_residuals[,3] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var3, main = "CUSUM of Squares for dlogUSINDPROD", xlab = "Observation", ylab = "Cumulative Sum of Squares")

# For ddlogUSINV
cusum_test_var4 <- efp(var_residuals[,4] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var4, main = "Rec-CUSUM for ddlogUSINV", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var4 <- efp(var_residuals[,4] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var4, main = "CUSUM of Squares for ddlogUSINV", xlab = "Observation", ylab = "Cumulative Sum of Squares")

# For dlogDXY
cusum_test_var5 <- efp(var_residuals[,5] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var5, main = "Rec-CUSUM for dlogDXY", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var5 <- efp(var_residuals[,5] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var5, main = "CUSUM of Squares for dlogDXY", xlab = "Observation", ylab = "Cumulative Sum of Squares")

# For dlogUSEPU
cusum_test_var6 <- efp(var_residuals[,6] ~ 1, type = "Rec-CUSUM")
plot(cusum_test_var6, main = "Rec-CUSUM for dlogUSEPU", xlab = "Observation", ylab = "Cumulative Sum")
cusum_squares_test_var6 <- efp(var_residuals[,6] ~ 1, type = "OLS-CUSUM")
plot(cusum_squares_test_var6, main = "CUSUM of Squares for dlogUSEPU", xlab = "Observation", ylab = "Cumulative Sum of Squares")