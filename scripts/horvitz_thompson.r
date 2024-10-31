#Density of Goals Scored by Edmonton Oilers
# Load EDM dataset
edm_data <- read.csv("data/EDM.csv", 
                     header=TRUE, stringsAsFactors=TRUE)

# Filter dataset to include only rows where 'situation' is 'all'
new_edm <- subset(edm_data, situation == "all")

# Separate goals scored by location: 'HOME' and 'AWAY'
home <- new_edm$goalsFor[new_edm$home_or_away == "HOME"]
away <- new_edm$goalsFor[new_edm$home_or_away == "AWAY"]

# Calculate density of goals scored at home and away games, with bandwidth 0.5
d_home <- density(home, bw = 0.5)
d_away <- density(away, bw = 0.5)

# Calculate mean goals scored for home and away games
mean_home <- mean(home)
mean_away <- mean(away)

# Plot density of goals scored at home games
plot(d_home, col = "blue", lwd = 2, main = "Density of Goals Scored by Edmonton Oilers",
     xlab = "Goals Scored", ylab = "Density", xlim = range(c(d_home$x, d_away$x)), 
     ylim=range(c(d_home$y, d_away$y)))

# Add density line for goals scored at away games
lines(d_away, col = "red", lwd = 2)

# Add vertical line at mean goals scored for home games
abline(v = mean_home, col = "blue", lwd = 2, lty = 2)

# Display mean goals scored at home games as text on the plot
text(x = mean_home, y = 0, labels = paste(round(mean_home, 4)), col = "blue", pos = 4)

# Add vertical line at mean goals scored for away games
abline(v = mean_away, col = "red", lwd = 2, lty = 2)

# Display mean goals scored at away games as text on the plot
text(x = mean_away, y = 0, labels = paste(round(mean_away, 4)), col = "red", pos = 2)

# Add legend indicating colors and line types for home and away games
legend("topright", legend = c("Home Games", "Away Games"),
       col = c("blue", "red"), lwd = 2, lty = 1)


# Load the dataset containing games watched and remove duplicate entries
gw <- read.csv("data/gamesWatched.csv", 
               header = TRUE, stringsAsFactors = TRUE)
unique_gw <- unique(gw)

# Separate goals scored by location: 'HOME' and 'AWAY' games
home <- new_edm$goalsFor[new_edm$home_or_away == "HOME"]
away <- new_edm$goalsFor[new_edm$home_or_away == "AWAY"]

# Calculate total number of home and away games
N_h <- length(home)
N_a <- length(away)

# Calculate the mean goals scored for home and away games
mean_home <- mean(home)
mean_away <- mean(away)

# Merge the games watched data with the main dataset on 'gameDate' to filter watched games only
f_data <- merge(new_edm, gw, by = "gameDate", all.x = FALSE, all.y = TRUE)

# Remove duplicates from the merged data
u_f_data <- unique(f_data)

# Separate goals scored by Edmonton in watched home and away games (with duplicates)
home_f <- f_data$goalsFor[f_data$home_or_away == "HOME"]
away_f <- f_data$goalsFor[f_data$home_or_away == "AWAY"]

# Separate goals scored by Edmonton in watched home and away games (no duplicates)
home_fdup <- u_f_data$goalsFor[u_f_data$home_or_away == "HOME"]
away_fdup <- u_f_data$goalsFor[u_f_data$home_or_away == "AWAY"]

# Calculate number of home and away games in watched dataset
n_h <- length(home_f)
n_a <- length(away_f)

# Calculate number of unique home and away games in watched dataset
n_h_nodup <- length(home_fdup)
n_a_nodup <- length(away_fdup)

# Calculate inclusion probabilities for home and away games in the unique watched dataset
pi_h <- rep(1 - ((N_h - 1) / N_h) ^ n_h, n_h_nodup)
pi_a <- rep(1 - ((N_a - 1) / N_a) ^ n_a, n_a_nodup)

# Compute Horvitz-Thompson weights for unique home and away games
t_u_home <- home_fdup / N_h
t_u_away <- away_fdup / N_a

# Calculate Horvitz-Thompson estimates for mean goals scored at home and away games
ht_home <- sum(t_u_home / pi_h)
ht_away <- sum(t_u_away / pi_a)

# Display Horvitz-Thompson estimates
print(ht_home)
print(ht_away)

# Create joint inclusion probability matrix for home games
pi_uv_home <- matrix(1 - 2 * ((N_h - 1) / N_h)^n_h + ((N_h - 2) / N_h)^n_h, 
                     nrow = n_h_nodup, ncol = n_h_nodup)
diag(pi_uv_home) <- pi_h

# Create joint inclusion probability matrix for away games
pi_uv_away <- matrix(1 - 2 * ((N_a - 1) / N_a)^n_a + ((N_a - 2) / N_a)^n_a, 
                     nrow = n_a_nodup, ncol = n_a_nodup)
diag(pi_uv_away) <- pi_a

# Function to estimate variance using Horvitz-Thompson estimation
estVarHT <- function(t_u, pi_u, pi_uv){
  delta <- pi_uv - outer(pi_u, pi_u)
  estimateVar <- sum((delta / pi_uv) * outer(t_u / pi_u, t_u / pi_u))
  return(abs(estimateVar))
}

# Calculate variance of Horvitz-Thompson estimates for home and away games
var_ht_home <- estVarHT(t_u_home, pi_h, pi_uv_home)
var_ht_away <- estVarHT(t_u_away, pi_a, pi_uv_away)

# Calculate standard errors for home and away Horvitz-Thompson estimates
se_ht_home <- sqrt(var_ht_home)
se_ht_away <- sqrt(var_ht_away)

# Display standard errors
print(se_ht_home)
print(se_ht_away)

# Calculate the difference in Horvitz-Thompson estimates between home and away games
d_hat <- ht_home - ht_away

# Calculate standard error of the difference in estimates
se_diff <- sqrt(se_ht_home^2 + se_ht_away^2)

# Calculate 95% confidence interval for the difference in mean goals scored
ci_diff <- d_hat + 2 * c(-1, 1) * se_diff

# Display confidence interval for difference
print(ci_diff)

# Conclusion based on confidence interval
# Since 0 is part of this interval, there is significant difference between goals scored in home and away indicating home ice advantage is real.
