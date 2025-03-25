# Loading required libraries
library(haven)
library(vars)
library(zoo)
library(ggplot2)
library(reshape2)
library(gganimate)

# Loading Stata dataset
Model_3 <- read_dta("/Users/saimanishprabhakar/Desktop/Dissertation/Dissertation Data and Estimation/STATA Dissertation files/STATA data files/Model_3 (USEPU Endog & VIX Exog).dta")

# Checking for missing values and remove them if any
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

# Defining Function to calculate TVGC for all combinations and return a p-value matrix
tvgc_heatmap <- function(data, window_size = 60, max_lag = 1) {
  n_vars <- ncol(data) # Number of variables
  var_names <- colnames(data)
  
  # Initializing an empty matrix to store p-values
  p_values_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
  rownames(p_values_matrix) <- var_names
  colnames(p_values_matrix) <- var_names
  
  # Looping through all combinations of variables
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i != j) { # Skipping the case where i == j (essentially, a variable cannot Granger-cause itself)
        # Extracting the two variables for testing
        X <- data[, i]
        Y <- data[, j]
        
        # Performing rolling Granger causality (using rolling window)
        p_values <- numeric(nrow(data) - window_size + 1)
        
        for (k in 1:(nrow(data) - window_size + 1)) {
          # Selecting the window data
          X_window <- X[k:(k + window_size - 1)]
          Y_window <- Y[k:(k + window_size - 1)]
          
          # Fitting VAR model within the window
          var_data <- cbind(X_window, Y_window)
          colnames(var_data) <- c("X", "Y")
          var_model <- VAR(var_data, p = max_lag)
          
          # Performing Granger causality test
          gc_result <- causality(var_model, cause = "X")
          
          # Storing the p-value of the rolling window
          p_values[k] <- gc_result$Granger$p.value
        }
        
        # Storing the minimum p-value for each combination
        p_values_matrix[i, j] <- min(p_values)
      }
    }
  }
  
  # Printing the raw p-values matrix
  print("Raw p-values matrix for TVGC Heatmap:")
  print(p_values_matrix)
  
  return(p_values_matrix)
}

# Defining Function to calculate rolling window Granger causality
calculate_rolling_gc <- function(data, cause_vars, effect_var, dates, window_size = 60, max_lag = 1) {
  n_obs <- nrow(data)
  n_vars <- length(cause_vars)
  
  # Initializing the results matrix
  results <- matrix(NA, nrow = n_obs - window_size + 1, ncol = n_vars)
  colnames(results) <- cause_vars
  
  for (i in 1:(n_obs - window_size + 1)) {
    window_data <- data[i:(i + window_size - 1), c(cause_vars, effect_var)]
    
    for (j in 1:n_vars) {
      var_data <- window_data[, c(cause_vars[j], effect_var)]
      var_model <- VAR(var_data, p = max_lag)
      gc_result <- causality(var_model, cause = cause_vars[j])
      results[i, j] <- gc_result$Granger$p.value
    }
  }
  
  # Adding dates to the results
  results_df <- as.data.frame(results)
  results_df$Date <- dates[window_size:n_obs]
  
  # Printing the raw p-values only for WTI
  if (effect_var == "dlogWTI") {
    print("Raw p-values for rolling Granger causality (WTI):")
    print(results_df)
  }
  
  return(results_df)
}

# Modify the plot_rolling_gc function to include animation option
plot_rolling_gc <- function(gc_results, response_var, animate = FALSE) {
  gc_long <- reshape2::melt(gc_results, id.vars = "Date", variable.name = "Impulse", value.name = "p_value")
  
  subtitle <- paste0("Variables Granger causing ", response_var)
  
  p <- ggplot(gc_long, aes(x = Date, y = p_value, color = Impulse)) +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_y_continuous(trans = "log10") +
    ggtitle(paste0("Rolling Window Analysis of Time-Varying Granger Causality\n",
                   subtitle)) +
    theme_minimal() +
    scale_color_brewer(palette = "Set2") +  # Move scale_color_brewer after theme
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 12)) +
    annotate("text", x = max(gc_long$Date), y = 0.05, label = "α = 0.05",
             hjust = 1, vjust = -0.5, color = "red")
  
  if (animate) {
    p <- p +
      transition_reveal(Date) +
      ease_aes('linear')
    
    anim <- animate(p, nframes = 200, fps = 10, width = 1000, height = 600)
    return(anim)
  } else {
    return(p)
  }
}

# Process each variable
process_variable("dlogWTI", animate = TRUE)  
process_variable("dlogUSFPCO", animate = TRUE)
process_variable("dlogUSINDPROD")
process_variable("ddlogUSINV")
process_variable("dlogDXY")
process_variable("dlogUSEPU")