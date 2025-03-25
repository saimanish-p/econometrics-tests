# Loading necessary libraries
library(tseries)
library(np)
library(vars)
library(haven)

# Loading Stata dataset
Model_3 <- read_dta("/Users/saimanishprabhakar/Desktop/Dissertation/Dissertation Data and Estimation/STATA Dissertation files/STATA data files/Model_3 (USEPU Endog & VIX Exog).dta")

# Checking for missing values and removing them if any
colSums(is.na(Model_3))
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

# Assigning column names to the endogenous and exogenous variables
colnames(endog_vars) <- c("dlogWTI", "dlogUSFPCO", "dlogUSINDPROD", "ddlogUSINV", "dlogDXY", "dlogUSEPU")
colnames(exog_vars) <- c("dlogVIX", "L1_dlogVIX", "L2_dlogVIX")

# Estimating the VAR model with 4 lags and exogenous variables
VAR_model <- VAR(endog_vars, p = 4, exogen = exog_vars)

# Extracting residuals from the VAR model
var_residuals <- residuals(VAR_model)

# Defining Diks-Panchenko test function
diks_panchenko_test <- function(X, Y, lag = 1) {
  # Ensuring the length of the series is the same
  n <- length(X)
  
  # Lagging the Y variable by 'lag'
  Y_lag <- c(rep(NA, lag), Y[1:(n - lag)])
  
  # Removing NAs caused by lagging
  X <- X[-(1:lag)]
  Y_lag <- Y_lag[-(1:lag)]
  
  # Removing any remaining NAs
  valid_indices <- complete.cases(X, Y_lag)
  X <- X[valid_indices]
  Y_lag <- Y_lag[valid_indices]
  
  # Ensuring data is not empty
  if (length(X) == 0 || length(Y_lag) == 0) {
    stop("No valid data after removing NAs")
  }
  
  # Creating data frames for bandwidth selection
  data_X <- data.frame(X = X)
  data_Y <- data.frame(Y_lag = Y_lag)
  
  # Performing nonparametric regression using npreg
  bw_X <- np::npregbw(Y_lag ~ X, data = data.frame(X = X, Y_lag = Y_lag), regtype = "lc", bwmethod = "cv.aic")
  bw_Y <- np::npregbw(X ~ Y_lag, data = data.frame(X = X, Y_lag = Y_lag), regtype = "lc", bwmethod = "cv.aic")
  
  np_X <- np::npreg(bws = bw_X, txdat = data_X, tydat = Y_lag)
  np_Y <- np::npreg(bws = bw_Y, txdat = data_Y, tydat = X)
  
  # Calculating the fitted values (predictions)
  pred_X <- fitted(np_X)
  pred_Y <- fitted(np_Y)
  
  # Calculating prediction errors
  diff_errors <- (X - pred_X)^2 - (Y_lag - pred_Y)^2
  
  # Applying a simple t-test to check the significance of the difference
  test_statistic <- mean(diff_errors) / (sd(diff_errors) / sqrt(length(diff_errors)))
  
  # Calculating p-value
  p_value <- 2 * (1 - pnorm(abs(test_statistic)))  # Two-tailed p-value
  
  return(list(statistic = test_statistic, p_value = p_value))
}

# Defining function to test all combinations
test_all_combinations <- function(residuals, lag = 1) {
  var_names <- colnames(residuals)
  n_vars <- ncol(residuals)
  results <- list()
  
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      X <- residuals[, i]
      Y <- residuals[, j]
      
      # Testing X -> Y
      result_XY <- diks_panchenko_test(X, Y, lag)
      results[[paste(var_names[i], "->", var_names[j])]] <- result_XY
      
      # Testing Y -> X
      result_YX <- diks_panchenko_test(Y, X, lag)
      results[[paste(var_names[j], "->", var_names[i])]] <- result_YX
    }
  }
  
  return(results)
}

# Performing the Diks-Panchenko test for all combinations
all_results <- test_all_combinations(var_residuals, lag = 1)

# Printing the results
for (combination in names(all_results)) {
  cat("\nTesting:", combination, "\n")
  cat("Diks-Panchenko Test Statistic:", all_results[[combination]]$statistic, "\n")
  cat("p-value:", all_results[[combination]]$p_value, "\n")
}