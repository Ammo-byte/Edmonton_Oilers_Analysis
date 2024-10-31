# Function to check convergence of theta values
testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) {
  sum(abs(thetaNew - thetaOld)) < if (relative)
    tolerance * sum(abs(thetaOld)) else tolerance
}

# Empirical Cumulative Distribution Function (ECDF) calculation
ecdf <- function(y) {
  n <- length(y)               # Number of elements in y
  ecdf_val <- numeric(n)       # Initialize vector to store ECDF values
  for (i in 1:length(y)) {
    a <- ifelse(y <= y[i], 1, 0)   # Indicator for values in y that are <= y[i]
    mn <- mean(a)                  # Mean gives the ECDF value at y[i]
    ecdf_val[i] <- mn              # Store ECDF value for current y[i]
  }
  return(ecdf_val)                 # Return vector of ECDF values
}

# Create Minimum Distance Estimation (MDE) psi and psiPrime functions
createMDEPsiFns <- function(y) {
  n <- length(y)               # Length of y
  ecdf_v <- ecdf(y)            # Calculate ECDF for y
  
  # Define the psi function to calculate distance between empirical and Poisson CDFs
  psi <- function(theta) {
    poisson_cdf <- ppois(y,theta)        # Poisson CDF for y with parameter theta
    poisson_m1 <- ppois(y-1, theta)      # Offset Poisson CDF for y - 1
    (2 / n) * sum((poisson_cdf - ecdf_v) * (poisson_m1 - poisson_cdf))  # Distance calculation
  }
  
  # Define psiPrime, the derivative of psi function, for Newton's method
  psiPrime <- function(theta) {
    poisson_cdf <- ppois(y, theta)         # Poisson CDF for y
    poisson_m1 <- ppois(y - 1, theta)      # Poisson CDF for y - 1
    poisson_m2 <- ppois(y - 2, theta)      # Poisson CDF for y - 2
    
    # Calculate the derivative based on differences between shifted Poisson CDF values
    (2 / n) * sum((poisson_m1 - poisson_cdf)^2 + 
                    (poisson_cdf - ecdf_v) * (poisson_m2 - 2 * poisson_m1 + poisson_cdf))
  }
  
  list(psi = psi, psiPrime = psiPrime)    # Return a list containing the psi and psiPrime functions
}

# Newton's method to find theta that minimizes psi function using psi and psiPrime functions
Newton <- function(theta = 0,
                   psiFn, 
                   psiPrimeFn, 
                   testConvergenceFn = testConvergence,
                   maxIterations = 100, 
                   tolerance = 1E-6,    
                   relative = FALSE) {

  converged <- FALSE           # Initialize convergence status as FALSE
  i <- 0                       # Initialize iteration counter
  
  # Iterate until convergence or reaching the maximum iterations
  while (!converged & i <= maxIterations) {
    thetaNew <- theta - psiFn(theta) / psiPrimeFn(theta)  # Update theta based on Newton's method
    
    # Check if theta has converged based on tolerance
    converged <- testConvergenceFn(thetaNew, theta,
                                   tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew          # Update theta for the next iteration
    i <- i + 1                 # Increment iteration counter
  }
  
  list(theta = theta,          # Return results as a list containing theta
       converged = converged,  # convergence status
       iteration = i,          # iteration count
       fnValue = psiFn(theta)) # final psi value
} 

# Load EDM dataset, filter for rows with 'all' situation, and select 'goalsFor' column
edm_data <- read.csv("data/EDM.csv", 
                     header=TRUE, stringsAsFactors=TRUE)
new_edm <- subset(edm_data, situation == "all")   # Filter rows where situation is 'all'
gF <-new_edm$goalsFor                             # Extract 'goalsFor' column
set.seed(341)                                     # Set seed for reproducibility
s <- sample(gF, size = 50, replace = FALSE)       # Draw a random sample of size 50
true_theta_mle <- mean(gF)                        # Calculate Maximum Likelihood Estimate (MLE)

psis <- createMDEPsiFns(gF)                       # Generate psi functions for MDE
result <- Newton(theta = true_theta_mle, psiFn = psis$psi, psiPrimeFn = psis$psiPrime) # Newton for MDE

n <- 1000                         # Number of random samples
thou_samples_mle <- numeric(n)    # Vector to store MLE estimates
thou_samples_mde <- numeric(n)    # Vector to store MDE estimates

# Loop to generate samples and calculate MLE and MDE for each sample
for (i in 1:n) {
  s <- sample(gF, size = 50, replace = FALSE)     # Draw a random sample of size 50
  psis_s <- createMDEPsiFns(s)                    # Generate psi functions for the sample
  theta_mle <- mean(s)                            # Calculate MLE as the sample mean
  theta_mde <- Newton(theta = theta_mle, psiFn = psis_s$psi, psiPrimeFn = psis_s$psiPrime)
  thou_samples_mle[i] <- theta_mle                # Store MLE for current sample
  thou_samples_mde[i] <- theta_mde$theta          # Store MDE for current sample
}

# Calculate bias, variance, and mean squared error (MSE) for MLE
bias_mle <- mean(thou_samples_mle)-true_theta_mle   # Bias of MLE
print(bias_mle)                                     # Print Bias
bias_mde <- mean(thou_samples_mde)-true_theta_mde   # Bias of MDE
print(bias_mde)                                     # Print Bias

variance_mle <- var(thou_samples_mle)               # Variance of MLE
print(variance_mle)                                 # Print Variance
variance_mde <- var(thou_samples_mde)               # Variance of MDE
print(variance_mde)                                 # Print Variance

mse_mle <- variance_mle + (bias_mle)^2              # Mean Squared Error (MSE) for MLE
print(mse_mle)                                      # Print MSE
mse_mde <- variance_mde + (bias_mde)^2              # Mean Squared Error (MSE) for MDE
print(mse_mde)                                      # Print MSE

# Conlusion: 
# The bias of the mle is closer to zero, the variance of mle is lower and mse is smaller for mle aswell. 
# Keeping these in mind and looking at the plots, mle seems to give a better estimation.