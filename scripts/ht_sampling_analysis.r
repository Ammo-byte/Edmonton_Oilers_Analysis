# Function to calculate Horvitz-Thompson estimate, standard error, and confidence interval
hte_ci <- function(sample, N) {
  n <- length(sample)  # Sample size
  pi <- rep(n / N, n)  # Inclusion probability for each element in sample
  ht_e <- sum(sample / pi) / N  # Horvitz-Thompson estimate for the population mean
  
  # Calculate the standard error of the Horvitz-Thompson estimate
  se_ht <- sqrt(sum(((1 - pi) / (pi^2 * N^2)) * (sample - ht_e)^2))
  
  # Construct a 95% confidence interval
  ci <- ht_e + 1.96 * c(-1, 1) * se_ht
  
  # Return estimate and confidence interval as a list
  return(list(e = ht_e, ci = ci))
}

# Load the Edmonton Oilers game data
edm_data <- read.csv("data/EDM.csv", 
                     header=TRUE, stringsAsFactors=TRUE)
new_edm <- subset(edm_data, situation == "all")  # Filter dataset to retain rows with 'situation' as 'all'

# Load watched games data and remove duplicates
gw <- read.csv("data/gamesWatched.csv", 
               header = TRUE, stringsAsFactors = TRUE)
unique_gw <- unique(gw)  # Unique entries in watched games data

# Extract goalsFor column from main dataset and calculate total games
gF <- new_edm$goalsFor
N <- length(gF)  # Total number of games (population size)
popmean <- mean(gF)  # Population mean for comparison

# Initialize a dataframe to store results across sample sizes
results <- data.frame(n = integer(), bias = numeric(), variance = numeric(),
                      mse = numeric(), coverage = numeric())

# Loop over sample sizes (from 100 to 1000 in increments of 100)
for (n in seq(100, 1000, by = 100)) {
  ht_e <- numeric(1000)  # Vector to store Horvitz-Thompson estimates
  ci <- matrix(NA, nrow = 1000, ncol = 2)  # Matrix to store confidence intervals
  
  # Perform 1000 simulations for each sample size
  for (i in 1:1000) {
    # Take a sample without replacement from the population
    sample <- sample(gF, n, replace = FALSE)
    
    # Calculate Horvitz-Thompson estimate and confidence interval
    list_all <- hte_ci(sample, N)
    ht_e[i] <- list_all$e  # Store the estimate
    ci[i, ] <- list_all$ci  # Store the confidence interval
  }
  
  # Calculate coverage: proportion of intervals containing the population mean
  coverage <- apply(X = ci, MARGIN = 1, FUN = function(u) {
    popmean >= u[1] & popmean <= u[2]
  })
  
  # Calculate bias, variance, and mean squared error for the current sample size
  bias <- mean(ht_e) - popmean
  variance <- var(ht_e)
  mse <- mean((ht_e - popmean)^2)
  coverage_r <- mean(coverage)
  
  # Append results to the dataframe
  results <- rbind(results, data.frame(n = n, bias = bias, variance = variance,
                                       mse = mse, coverage = coverage_r))
}

# Plot results for each metric across sample sizes
par(mfrow = c(1, 4))

# Plot 1: Bias vs. Sample Size
plot(results[, 1], results[, 2], type = "l", ylim = c(-0.01, 0.01), 
     xlab = "Sample Size (n)", ylab = "Bias", main = "Bias vs. n", col = 'lightblue')
abline(h = 0, lty = 2)  # Add horizontal line at 0 for reference

# Plot 2: Variance vs. Sample Size
plot(results[, 1], results[, 3], type = "l", ylim = c(0, 0.03), 
     xlab = "Sample Size (n)", ylab = "Variance", main = "Variance vs. n", col ='lightgreen')

# Plot 3: MSE vs. Sample Size
plot(results[, 1], results[, 4], type = "l", ylim = c(0, 0.03), 
     xlab = "Sample Size (n)", ylab = "MSE", main = "MSE vs. n", col = 'orange')

# Plot 4: Coverage vs. Sample Size
plot(results[, 1], results[, 5], type = "l", ylim = c(0, 1), 
     xlab = "Sample Size (n)", ylab = "Coverage", main = "Coverage vs. n", col = 'purple')
abline(h = 0.95, lty = 2)  # Add horizontal line at 0.95 for reference (ideal coverage)

# From the first plot we can see that as sample size gets larger, bias gets closer 
# to zero. Furthermore as sample size increases, variance and MSE decrease. 
# The coverage however stays around .95 the entire time. So larger sample sizes 
# are more accurate for HT estimator.