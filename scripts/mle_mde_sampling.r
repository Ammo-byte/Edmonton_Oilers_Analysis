# Function to test if the convergence criteria is met
testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) {
  # Calculate absolute difference and compare with tolerance
  sum(abs(thetaNew - thetaOld)) < if (relative)
    # For relative tolerance, scale by thetaOld's magnitude
    tolerance * sum(abs(thetaOld)) else tolerance
}

# Function to compute Empirical Cumulative Distribution Function (ECDF)
ecdf <- function(y) {
  n <- length(y)  # Get length of vector y
  ecdf_val <- numeric(n)  # Initialize an empty vector to store ECDF values
  
  # Calculate ECDF at each point in y
  for (i in 1:length(y)) {
    # Count all values <= current value y[i] and take the mean
    a <- ifelse(y <= y[i], 1, 0)
    mn <- mean(a)
    ecdf_val[i] <- mn  # Store ECDF value
  }
  return(ecdf_val)  # Return vector of ECDF values
}

# Factory function to create the psi and psiPrime functions used in minimum distance estimation (MDE)
createMDEPsiFns <- function(y) {
  n <- length(y)   # Length of input data y
  ecdf_v <- ecdf(y)  # Calculate the empirical CDF for y
  
  # Define the psi function based on input theta
  psi <- function(theta) {
    poisson_cdf <- ppois(y, theta)         # Calculate Poisson CDF for each y given theta
    poisson_m1 <- ppois(y - 1, theta)      # Poisson CDF shifted by 1 for calculating derivative
    (2 / n) * sum((poisson_cdf - ecdf_v) * (poisson_m1 - poisson_cdf))  # Compute psi
  }
  
  # Define the psiPrime function, the derivative of psi, for use in Newton's method
  psiPrime <- function(theta) {
    poisson_cdf <- ppois(y, theta)
    poisson_m1 <- ppois(y - 1, theta)
    poisson_m2 <- ppois(y - 2, theta)
    
    # Calculate the derivative using a second-order difference approximation
    (2 / n) * sum((poisson_m1 - poisson_cdf)^2 + 
                    (poisson_cdf - ecdf_v) * (poisson_m2 - 2 * poisson_m1 + poisson_cdf))
  }
  
  # Return both psi and psiPrime as a list
  list(psi = psi, psiPrime = psiPrime)
}

# Newton's method function to find the theta value that minimizes psi function
Newton <- function(theta = 0,
                   psiFn, 
                   psiPrimeFn, 
                   testConvergenceFn = testConvergence,
                   maxIterations = 100, 
                   tolerance = 1E-6,    
                   relative = FALSE) {

  converged <- FALSE  # Initialize convergence status as FALSE
  i <- 0              # Initialize iteration counter
  
  # Iterate until convergence criteria is met or maximum iterations reached
  while (!converged & i <= maxIterations) {
    # Update theta using Newton's update rule
    thetaNew <- theta - psiFn(theta) / psiPrimeFn(theta)
    
    # Check if updated theta has converged
    converged <- testConvergenceFn(thetaNew, theta,
                                   tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew  # Update theta for next iteration
    i <- i + 1         # Increment iteration counter
  }
  
  # Return theta value, convergence status, iteration count, and final psi value
  list(theta = theta,
       converged = converged,
       iteration = i,
       fnValue = psiFn(theta))
} 

# Load Edmonton Oilers game data and filter by "all" situations
edm_data <- read.csv("dataq/EDM.csv", 
                     header=TRUE, stringsAsFactors=TRUE)
new_edm <- subset(edm_data, situation == "all")  # Filter dataset for "all" situations
gF <- new_edm$goalsFor  # Extract 'goalsFor' column for analysis

# Calculate Maximum Likelihood Estimate (MLE) for theta
true_theta_mle <- mean(gF)

# Generate psi and psiPrime functions using createMDEPsiFns function
psis <- createMDEPsiFns(gF)

# Use Newton's method to find Minimum Distance Estimate (MDE), starting at MLE
result <- Newton(theta = true_theta_mle, psiFn = psis$psi, psiPrimeFn = psis$psiPrime)
true_theta_mde <- result$theta  # Store the MDE value for theta

# Set up for simulation study
n <- 1000                     # Number of samples to simulate
thou_samples_mle <- numeric(n) # Store MLE values
thou_samples_mde <- numeric(n) # Store MDE values

# Perform 1000 simulations to estimate sampling distributions of MLE and MDE
for (i in 1:n) {
  s <- sample(gF, size = 50, replace = FALSE)  # Take random sample of size 50 from data
  psis_s <- createMDEPsiFns(s)                 # Create psi functions for sample
  theta_mle <- mean(s)                         # Calculate MLE for sample
  theta_mde <- Newton(theta = theta_mle, psiFn = psis_s$psi, psiPrimeFn = psis_s$psiPrime)
  thou_samples_mle[i] <- theta_mle             # Store sample MLE
  thou_samples_mde[i] <- theta_mde$theta       # Store sample MDE
}

# Set up plot for side-by-side histograms of MLE and MDE distributions
par(mfrow = c(1, 2))

# Round true MLE and MDE for display in plots
round_mle <- round(true_theta_mle, 4)
round_mde <- round(true_theta_mde, 4)

# Plot histogram for the sampling distribution of MLE
hist(thou_samples_mle, breaks = 15, main = "Sampling Distribution of MLE",
     xlab = "MLE Estimates", col = "lightblue", border = "black", ylim= c(0,200), xlim = c(1.5, 4))
# Add vertical line for true MLE value
abline(v = true_theta_mle, col = "red", lwd = 2)
# Display true MLE value on plot
text(x = round_mle, y = 120, labels = paste(round_mle), col = "red", pos = 4)

# Plot histogram for the sampling distribution of MDE
hist(thou_samples_mde, breaks = 15, main = "Sampling Distribution of MDE",
     xlab = "MDE Estimates", col = "lightgreen", border = "black", ylim = c(0,200), xlim = c(1.5,4))
# Add vertical line for true MDE value
abline(v = true_theta_mde, col = "blue", lwd = 2)
# Display true MDE value on plot
text(x = round_mde, y = 120, labels = paste(round_mde), col = "blue", pos = 4)