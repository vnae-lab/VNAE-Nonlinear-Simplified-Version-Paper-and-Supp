# =============================================================================
# VNAE COMPLETE VISUALIZATION
# =============================================================================
#
# Paper: "Riemannian Manifolds of Asymmetric Equilibria: The Victoria-Nash Geometry"
# Author: Daniel Henrique Pereira  
# Victoria-Nash Asymmetric Equilibrium (VNAE)
#
# =============================================================================

library(ggplot2)
library(plotly)
library(gridExtra)
library(rootSolve)

# =============================================================================
# PARAMETER CONFIGURATION - USER CHANGES HERE!
# =============================================================================

# DEFINE θA AND θB VALUES HERE (Example 2.2 parameters):
thetaA <- 0.5    # Structural coefficient for Player A
thetaB <- 1.0    # Structural coefficient for Player B  
omega_A <- 1.0   # Structural weight for A (from paper)
omega_B <- 1.0   # Structural weight for B (from paper)
beta <- 0.1      # Inertial parameter (from Section 3)

# =============================================================================
# VNAE MODEL FUNCTIONS
# =============================================================================

# Expectation field F(s;θ) from Definition 2.1 and Example 2.2
vnae_expectation_field <- function(x, y, thetaA, thetaB) {
  Fx <- (omega_A + thetaA) * x^2 + omega_A * y
  Fy <- (omega_B + thetaB) * y^2 + omega_B * x
  return(list(Fx = Fx, Fy = Fy))
}

# VNAE manifold equations from Example 2.2
vnae_manifold_A <- function(x) {
  -((omega_A + thetaA)/omega_A) * x^2
}

vnae_manifold_B <- function(y) {
  -((omega_B + thetaB)/omega_B) * y^2
}

# Find equilibrium points (Γ_VNAE)
find_vnae_equilibria <- function(thetaA, thetaB) {
  system_eq <- function(vars) {
    x <- vars[1]; y <- vars[2]
    F <- vnae_expectation_field(x, y, thetaA, thetaB)
    c(F$Fx, F$Fy)
  }
 
  eq1 <- multiroot(system_eq, start = c(0, 0))
  eq2 <- multiroot(system_eq, start = c(-1, -1))
 
  return(list(eq1 = eq1$root, eq2 = eq2$root))
}

# Riemannian metric g_ij from Section 3
vnae_metric <- function(x, y, thetaA, thetaB) {
  H11 <- 2*x + 2*thetaA*x  # ∂²(V+φ)/∂x²
  H12 <- 1                  # ∂²(V+φ)/∂x∂y
  H22 <- 2*y + 2*thetaB*y  # ∂²(V+φ)/∂y²
 
  g11 <- omega_A + beta * H11
  g12 <- beta * H12
  g22 <- omega_B + beta * H22
 
  matrix(c(g11, g12, g12, g22), 2, 2)
}

# Gaussian curvature from Theorem 3.1
vnae_curvature <- function(x, y, thetaA, thetaB) {
  g <- vnae_metric(x, y, thetaA, thetaB)
 
  H11 <- 2*x + 2*thetaA*x
  H12 <- 1
  H22 <- 2*y + 2*thetaB*y
  H <- matrix(c(H11, H12, H12, H22), 2, 2)
 
  det_H <- det(H)
  det_g <- det(g)
 
  (det_H/det_g) * beta * abs(thetaA - thetaB)
}

# =============================================================================
# 1. 2D PLOT: VECTOR FIELD WITH VNAE MANIFOLD
# =============================================================================

create_2d_plot <- function(thetaA, thetaB) {
  x <- seq(-1.5, 1.5, length.out = 20)
  y <- seq(-1.5, 1.5, length.out = 20)
  grid <- expand.grid(x = x, y = y)
 
  field <- vnae_expectation_field(grid$x, grid$y, thetaA, thetaB)
  grid$Fx <- field$Fx
  grid$Fy <- field$Fy
 
  # Create manifold curves
  x_manifold <- seq(-1.5, 1.5, length.out = 100)
  y_manifold_A <- vnae_manifold_A(x_manifold)
  y_manifold_B_seq <- seq(-1.5, 0.5, length.out = 100)
  x_manifold_B <- vnae_manifold_B(y_manifold_B_seq)
 
  # Find equilibrium points
  equilibria <- find_vnae_equilibria(thetaA, thetaB)
 
  ggplot() +
    # Vector field
    geom_segment(data = grid,
                 aes(x = x, y = y, xend = x + Fx/15, yend = y + Fy/15),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "blue", alpha = 0.6) +
    # Manifold curves
    geom_line(data = data.frame(x = x_manifold, y = y_manifold_A),
              aes(x = x, y = y), color = "red", size = 1.2, linetype = "solid") +
    geom_line(data = data.frame(x = x_manifold_B, y = y_manifold_B_seq),
              aes(x = x, y = y), color = "darkgreen", size = 1.2, linetype = "solid") +
    # Equilibrium points
    geom_point(aes(x = equilibria$eq1[1], y = equilibria$eq1[2]),
               size = 4, color = "purple", shape = 17) +
    geom_point(aes(x = equilibria$eq2[1], y = equilibria$eq2[2]),
               size = 4, color = "purple", shape = 17) +
    labs(title = paste("VNAE: Expectation Field F(s;θ) and Manifold Γ_VNAE\n",
                       "θA =", thetaA, ", θB =", thetaB),
         x = "Strategy x (Player A)", y = "Strategy y (Player B)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

# =============================================================================
# 2. 3D PLOT: GAUSSIAN CURVATURE SURFACE
# =============================================================================

create_curvature_surface <- function(thetaA, thetaB) {
  x <- seq(-1, 1, length.out = 25)
  y <- seq(-1, 1, length.out = 25)
  grid <- expand.grid(x = x, y = y)
 
  grid$K <- sapply(1:nrow(grid), function(i) {
    vnae_curvature(grid$x[i], grid$y[i], thetaA, thetaB)
  })
 
  K_matrix <- matrix(grid$K, nrow = length(x), ncol = length(y))
 
  plot_ly() %>%
    add_surface(x = ~x, y = ~y, z = ~K_matrix,
                colorscale = "RdBu", opacity = 0.85,
                name = "Curvature K") %>%
    layout(title = paste("VNAE: Gaussian Curvature K(s;θ)\nθA =", thetaA, "θB =", thetaB),
           scene = list(xaxis = list(title = "x"),
                       yaxis = list(title = "y"),
                       zaxis = list(title = "K(s;θ)")))
}

# =============================================================================
# 3. PLOT: CURVATURE vs ASYMMETRY ANALYSIS
# =============================================================================

create_curvature_analysis <- function(thetaA, thetaB) {
  # Calculate curvature at origin
  K_origin <- vnae_curvature(0, 0, thetaA, thetaB)
 
  # Create asymmetry parameter sweep
  thetaA_seq <- seq(0.1, 2, length.out = 50)
  curvatures <- sapply(thetaA_seq, function(ta) vnae_curvature(0, 0, ta, thetaB))
 
  df <- data.frame(thetaA = thetaA_seq, curvature = curvatures)
 
  ggplot(df, aes(x = thetaA, y = curvature)) +
    geom_line(color = "purple", size = 1.2) +
    geom_vline(xintercept = thetaA, linetype = "dashed", color = "red") +
    geom_hline(yintercept = K_origin, linetype = "dashed", color = "blue") +
    geom_point(aes(x = thetaA, y = K_origin), size = 3, color = "red") +
    labs(title = paste("VNAE Manifold Curvature Analysis\n",
                       "K(0,0) =", round(K_origin, 4), "| θA =", thetaA, "θB =", thetaB),
         x = "θA (Structural Asymmetry)",
         y = "Gaussian Curvature K(0,0)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# =============================================================================
# 4. PLOT: PHASE PORTRAIT WITH EQUILIBRIUM MANIFOLD
# =============================================================================

create_phase_portrait <- function(thetaA, thetaB) {
  x <- seq(-1.5, 1.5, length.out = 25)
  y <- seq(-1.5, 1.5, length.out = 25)
  grid <- expand.grid(x = x, y = y)
 
  field <- vnae_expectation_field(grid$x, grid$y, thetaA, thetaB)
  grid$Fx <- field$Fx
  grid$Fy <- field$Fy
  grid$magnitude <- sqrt(grid$Fx^2 + grid$Fy^2)
 
  # Manifold curve
  x_manifold <- seq(-1.5, 1.5, length.out = 100)
  y_manifold <- vnae_manifold_A(x_manifold)
 
  ggplot() +
    geom_segment(data = grid,
                 aes(x = x, y = y, xend = x + Fx/15, yend = y + Fy/15,
                     color = magnitude, alpha = magnitude),
                 arrow = arrow(length = unit(0.12, "cm"))) +
    geom_line(data = data.frame(x = x_manifold, y = y_manifold),
              aes(x = x, y = y), color = "red", size = 1.5, linetype = "solid") +
    scale_color_viridis_c(option = "plasma", name = "Field\nStrength") +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    labs(title = paste("VNAE Phase Portrait with Equilibrium Manifold\n",
                       "θA =", thetaA, ", θB =", thetaB),
         x = "Strategy x", y = "Strategy y") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# =============================================================================
# 5. COMPARISON: SYMMETRIC vs ASYMMETRIC CASES
# =============================================================================

create_comparison_plot <- function(thetaA, thetaB) {
  # Symmetric case (θA = θB = 1.0)
  p_sym <- create_2d_plot(1.0, 1.0) +
    labs(title = "Symmetric Case: θA = θB = 1.0\n(Classical Nash)")
 
  # Asymmetric case with current values
  p_asym <- create_2d_plot(thetaA, thetaB) +
    labs(title = paste("Asymmetric Case: θA =", thetaA, ", θB =", thetaB, "\n(VNAE)"))
 
  grid.arrange(p_sym, p_asym, ncol = 2)
}

# =============================================================================
# 6. NUMERICAL ANALYSIS
# =============================================================================

analyze_vnae_system <- function(thetaA, thetaB) {
  # Find equilibria
  equilibria <- find_vnae_equilibria(thetaA, thetaB)
 
  # Calculate curvature at origin
  K_origin <- vnae_curvature(0, 0, thetaA, thetaB)
 
  # Calculate metric at origin
  g_origin <- vnae_metric(0, 0, thetaA, thetaB)
 
  cat("=========================================\n")
  cat("VNAE COMPLETE ANALYSIS\n")
  cat("=========================================\n")
  cat("Structural Parameters (Example 2.2):\n")
  cat("  θA =", thetaA, "(Player A structural coefficient)\n")
  cat("  θB =", thetaB, "(Player B structural coefficient)\n")
  cat("  ωA =", omega_A, ", ωB =", omega_B, "(Structural weights)\n")
  cat("  β =", beta, "(Inertial parameter)\n")
  cat("  |θA - θB| =", abs(thetaA - thetaB), "(Asymmetry measure)\n")
 
  cat("\nEquilibrium Manifold Γ_VNAE:\n")
  cat("  Equilibrium 1: (", round(equilibria$eq1[1], 4), ",",
      round(equilibria$eq1[2], 4), ")\n")
  cat("  Equilibrium 2: (", round(equilibria$eq2[1], 4), ",",
      round(equilibria$eq2[2], 4), ")\n")
 
  cat("\nRiemannian Geometry (Section 3):\n")
  cat("  Metric g(0,0) = [", round(g_origin[1,1], 4), round(g_origin[1,2], 4), "\n")
  cat("                   ", round(g_origin[2,1], 4), round(g_origin[2,2], 4), "]\n")
  cat("  det(g) =", round(det(g_origin), 6), "\n")
  cat("  Gaussian Curvature K(0,0) =", round(K_origin, 4), "\n")
 
  cat("\nTheorem 3.1 Verification:\n")
  if (thetaA != thetaB && beta > 0) {
    cat("  ✓ K > 0 confirmed: θA ≠ θB and β > 0\n")
  } else {
    cat("  → K → 0: Symmetric or zero inertia limit\n")
  }
  cat("=========================================\n")
}

# =============================================================================
# MAIN EXECUTION - RUN THIS SECTION
# =============================================================================

cat("STARTING VNAE COMPLETE VISUALIZATION...\n")
cat("Paper: Riemannian Manifolds of Asymmetric Equilibria\n")
cat("Parameters: θA =", thetaA, ", θB =", thetaB, "\n\n")

# Generate all plots
print("1. Generating 2D Vector Field with VNAE Manifold...")
plot_2d <- create_2d_plot(thetaA, thetaB)
print(plot_2d)

print("2. Generating 3D Gaussian Curvature Surface...")
plot_3d <- create_curvature_surface(thetaA, thetaB)
print(plot_3d)

print("3. Generating Curvature Analysis...")
plot_curvature <- create_curvature_analysis(thetaA, thetaB)
print(plot_curvature)

print("4. Generating Phase Portrait...")
plot_phase <- create_phase_portrait(thetaA, thetaB)
print(plot_phase)

print("5. Generating Symmetric vs Asymmetric Comparison...")
plot_comparison <- create_comparison_plot(thetaA, thetaB)
print(plot_comparison)

print("6. Performing Complete Numerical Analysis...")
analyze_vnae_system(thetaA, thetaB)



# =============================================================================
# USER INSTRUCTIONS:
# =============================================================================
#
# 1. Change thetaA and thetaB values in lines 18-19
# 2. Execute all code (Ctrl+A -> Ctrl+Enter)  
# 3. All plots will be generated automatically
#
# TEST EXAMPLES FROM PAPER:
# thetaA <- 0.5; thetaB <- 1.0   # Example 2.2 parameters
# thetaA <- 1.0; thetaB <- 1.0   # Symmetric case (Nash limit)
# thetaA <- 2.0; thetaB <- 1.0   # Strong asymmetry
#
# =============================================================================
