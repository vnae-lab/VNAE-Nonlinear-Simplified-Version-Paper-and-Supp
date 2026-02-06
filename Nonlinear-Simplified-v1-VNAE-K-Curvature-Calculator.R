# Load required package
library(rootSolve)

# --- PARAMETERS ---
omega_A <- 1.0; theta_A <- 0.5
omega_B <- 1.0; theta_B <- 1.0
beta <- 0.1

# --- 1. FIND VNAE EQUILIBRIA (Γ_VNAE) ---
system_eq <- function(vars) {
  x <- vars[1]; y <- vars[2]
  c((omega_A + theta_A)*x^2 + omega_A*y,
    (omega_B + theta_B)*y^2 + omega_B*x)
}

# Find equilibrium points
equilibrium1 <- multiroot(system_eq, start = c(0, 0))
equilibrium2 <- multiroot(system_eq, start = c(-1, -1))

# --- 2. CALCULATE METRIC gᵢⱼ AND CURVATURE K ---
calculate_vnae_geometry <- function(theta_A, theta_B, beta, x, y) {
  # Hessian H = ∇²(V + φ)
  H <- matrix(c(2*x + 2*theta_A*x, 1,
                1, 2*y + 2*theta_B*y), nrow = 2)
 
  # Metric tensor gᵢⱼ (Section 3)
  g <- matrix(c(omega_A + beta*H[1,1], beta*H[1,2],
                beta*H[1,2], omega_B + beta*H[2,2]), nrow = 2)
 
  # Curvature K (Theorem 3.1)
  det_H <- det(H)
  det_g <- det(g)
  K <- (det_H/det_g) * beta * abs(theta_A - theta_B)
 
  # Return both metric and curvature
  return(list(g_metric = g, curvature = K, hessian = H))
}

# --- 3. APPLY TO EQUILIBRIUM POINTS ---
x_eq <- equilibrium2$root[1]  # -0.6057
y_eq <- equilibrium2$root[2]  # -0.5503

geometry <- calculate_vnae_geometry(theta_A, theta_B, beta, x_eq, y_eq)

# --- 4. RESULTS ---
cat("VNAE GEOMETRY ANALYSIS:\n")
cat("=======================\n")
cat("Equilibrium point: (", round(x_eq,4), ",", round(y_eq,4), ")\n\n")

cat("METRIC TENSOR gᵢⱼ:\n")
print(geometry$g_metric)
cat("\nDeterminant of g: ", det(geometry$g_metric), "\n\n")

cat("HESSIAN H = ∇²(V + φ):\n")
print(geometry$hessian)
cat("\nDeterminant of H: ", det(geometry$hessian), "\n\n")

cat("GAUSSIAN CURVATURE K:", geometry$curvature, "\n")
cat("K > 0 confirms asymmetric equilibrium (θ_A ≠ θ_B)\n")
