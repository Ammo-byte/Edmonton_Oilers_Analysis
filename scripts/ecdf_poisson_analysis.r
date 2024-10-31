# Function to check convergence between two theta values
# based on absolute or relative tolerance level
testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) {
  # Check if the difference is within the tolerance level, return TRUE if converged
  sum(abs(thetaNew - thetaOld)) < if (relative)
    tolerance * sum(abs(thetaOld)) else tolerance
}

# Function to compute the Empirical Cumulative Distribution Function (ECDF)
# for a given vector y
ecdf <- function(y) {
  n <- length(y)
  ecdf_val <- numeric(n)  # Initialize an empty vector for ECDF values
  # Calculate ECDF at each point by averaging indicator values
  for (i in 1:length(y)) {
    a <- ifelse(y <= y[i], 1, 0)
    mn <- mean(a)
    ecdf_val[i] <- mn
  }
  return(ecdf_val)
}

# Factory function to create psi and psi' (psi prime) functions
# used for minimum distance estimation (MDE) based on Cramér-von Mises criterion
createMDEPsiFns <- function(y) {
  n <- length(y)
  ecdf_v <- ecdf(y)  # Calculate ECDF values for y
  
  # Define psi function based on theta
  psi <- function(theta) {
    poisson_cdf <- ppois(y, theta)
    poisson_m1 <- ppois(y - 1, theta)  # CDF shifted by 1 for calculating derivative
    # Sum of squared differences between empirical and theoretical CDFs
    (2 / n) * sum((poisson_cdf - ecdf_v) * (poisson_m1 - poisson_cdf))
  }
  
  # Define psiPrime function, the derivative of psi, used in Newton's method
  psiPrime <- function(theta) {
    poisson_cdf <- ppois(y, theta)
    poisson_m1 <- ppois(y - 1, theta)
    poisson_m2 <- ppois(y - 2, theta)
    
    # Second-order difference to approximate the derivative
    (2 / n) * sum((poisson_m1 - poisson_cdf)^2 + 
                  (poisson_cdf - ecdf_v) * (poisson_m2 - 2 * poisson_m1 + poisson_cdf))
  }
  
  # Return a list with both psi and psiPrime functions
  list(psi = psi, psiPrime = psiPrime)
}

# Newton's method function to estimate theta for given psi functions
Newton <- function(theta = 0,
                   psiFn, 
                   psiPrimeFn, 
                   testConvergenceFn = testConvergence,
                   maxIterations = 100, 
                   tolerance = 1E-6,    
                   relative = FALSE) {

  converged <- FALSE
  i <- 0
  # Iterate until convergence or maximum iterations
  while (!converged & i <= maxIterations) {
    
    # Update theta value using Newton's update rule
    thetaNew <- theta - psiFn(theta) / psiPrimeFn(theta)
    
    # Check for convergence based on the testConvergence function
    converged <- testConvergenceFn(thetaNew, theta,
                                   tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
  }
  # Return results in a list
  list(theta = theta,
       converged = converged,
       iteration = i,
       fnValue = psiFn(theta)
  )
} 

# Load the dataset
edm_data <- read.csv("data/EDM.csv", 
                     header = TRUE, stringsAsFactors = TRUE)

# Filter data for "all" situations only
new_edm <- subset(edm_data, situation == "all")
gF <- new_edm$goalsFor  # Extract 'goalsFor' column for analysis

# Maximum Likelihood Estimation (MLE) for theta
theta_mle <- mean(gF)

# Generate psi and psiPrime functions using the createMDEPsiFns factory function
psis <- createMDEPsiFns(gF)

# Estimate theta using Newton's method, starting from MLE
result <- Newton(theta = theta_mle, psiFn = psis$psi, psiPrimeFn = psis$psiPrime)

# Sort 'goalsFor' data and calculate ECDF for plotting
sort_gF <- sort(gF)
ecdf_v <- ecdf(sort_gF)

# Plot ECDF of 'goalsFor' data
plot(sort_gF, ecdf_v, type = "s", 
     main = "Empirical Cumulative Distribution Function for goalsFor",
     xlab = "goalsFor", ylab = "Cumulative Probability",
     col = "blue", lwd = 2)

# Recalculate MLE and calculate MDE from the Newton's method result
theta_mle <- mean(gF)
theta_mde <- result$theta

# Calculate Poisson CDFs based on MLE and MDE theta values
p_cdf <- ppois(sort_gF, lambda = theta_mle)     # CDF using MLE
p_cdfmde <- ppois(sort_gF, lambda = theta_mde)  # CDF using MDE

# Overlay the Poisson CDFs on the ECDF plot
lines(sort_gF, p_cdf, type = "s", col = "red", lwd = 2, lty = 2)       # MLE line
lines(sort_gF, p_cdfmde, type = "s", col = "purple", lwd = 2, lty = 2) # MDE line

# Add legend to distinguish between ECDF, MLE, and MDE lines
legend("bottomright", legend = c("ECDF", "Poisson CDF with MLE", "Poisson CDF with MDE"), 
       col = c("blue", "red", "purple"), lwd = 2, lty = c(1, 2))

# Conclusion:
# The poisson models are all somewhat similar in shape with some minor differences. 
# The MLE and MDE both seem to fit data well but when you look a bit closer the MLE
# seems to be slightly closer to the ECDF shape. 