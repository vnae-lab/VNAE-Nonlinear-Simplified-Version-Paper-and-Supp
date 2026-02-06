# Load required package for root solving
if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}
library(rootSolve)

# --- PARAMETERS ---
omega_A <- 1.0     # Structural weight for Player A
theta_A <- 0.5     # Structural coefficient for Player A  
omega_B <- 1.0     # Structural weight for Player B
theta_B <- 1.0     # Structural coefficient for Player B

# --- SYSTEM DEFINITION ---
# Victoria-Nash field F(x,y;θ) = 0
system_eq <- function(vars) {
  x <- vars[1]
  y <- vars[2]
 
  eq1 <- (omega_A + theta_A) * x^2 + omega_A * y
  eq2 <- (omega_B + theta_B) * y^2 + omega_B * x
 
  c(eq1, eq2)
}

# --- FIND EQUILIBRIUM POINTS ---
equilibrium1 <- multiroot(f = system_eq, start = c(0, 0))
equilibrium2 <- multiroot(f = system_eq, start = c(-1, -1))

# --- VISUALIZE VNAE MANIFOLD ---
# From Example 2.2: y = -((ω_A+θ_A)/ω_A)x² defines the VNAE manifold
vnae_manifold <- function(x) {
  -((omega_A + theta_A)/omega_A) * x^2
}

# Create sequence for plotting
x_vals <- seq(-1.5, 1.5, length.out = 200)
y_vals <- vnae_manifold(x_vals)

# Plot the Victoria-Nash manifold
plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Strategy x (Player A)",
     ylab = "Strategy y (Player B)",
     main = "Victoria-Nash Manifold Γ_VNAE(θ)",
     xlim = c(-1.5, 1.5), ylim = c(-2.5, 0.5))
grid()

# Mark the equilibrium points
points(equilibrium1$root[1], equilibrium1$root[2], col = "red", pch = 19, cex = 1.5)
points(equilibrium2$root[1], equilibrium2$root[2], col = "red", pch = 19, cex = 1.5)

# Add legend
legend("topright",
       legend = c(expression(paste(Gamma["VNAE"], "(θ)")), "VNAE Equilibria"),
       col = c("blue", "red"),
       lty = c(1, NA),
       pch = c(NA, 19),
       bg = "white")

# --- DISPLAY RESULTS ---
cat("Victoria-Nash Asymmetric Equilibria:\n")
cat("=====================================\n")
cat("Parameters: ω_A =", omega_A, "θ_A =", theta_A, "| ω_B =", omega_B, "θ_B =", theta_B, "\n")
cat("Asymmetry coefficient: |θ_A - θ_B| =", abs(theta_A - theta_B), "\n\n")
cat("Equilibrium points found:\n")
cat("1. (x, y) = (", round(equilibrium1$root[1], 4), ",", round(equilibrium1$root[2], 4), ")\n")
cat("2. (x, y) = (", round(equilibrium2$root[1], 4), ",", round(equilibrium2$root[2], 4), ")\n")
cat("Both points lie on the curved VNAE manifold.\n")
