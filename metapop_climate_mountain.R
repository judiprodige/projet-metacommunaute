###############################################################################
#  Métapopulation spatialement explicite en paysage montagneux               #
#  sous changement climatique                                                 #
#                                                                             #
#  Équation de Hanski (SPOM) :                                               #
#  dp_x/dt = c_x * sum_{y≠x} K(x,y) * p_y * A_y * (1-p_x) - (e_x/A_x)*p_x #
#                                                                             #
#  Scénarios : baseline, réchauffement modéré (+2°C), sévère (+4°C)          #
###############################################################################

rm(list = ls())

# ===========================================================================
# 0. Packages
# ===========================================================================
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",     quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",     quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(42)

# ===========================================================================
# 1. Paramètres du paysage montagneux
# ===========================================================================

# Nombre de patchs d'habitat
N <- 30

# Coordonnées spatiales (km) – distribués dans un massif montagneux
coords <- data.frame(
  x = runif(N, 0, 30),  # étendue E-O
  y = runif(N, 0, 30)   # étendue N-S
)

# Altitude (m) : modèle simplifié en cône avec bruit
centre <- c(15, 15)
dist_centre <- sqrt((coords$x - centre[1])^2 + (coords$y - centre[2])^2)
altitude <- 2500 - 80 * dist_centre + rnorm(N, 0, 100)
altitude <- pmax(altitude, 400)  # minimum 400 m

# Surface des patchs (ha) – décroît avec l'altitude (effet conique)
A_base <- pmax(0.5, 10 - 0.003 * altitude + rnorm(N, 0, 1))

# Température locale (°C) – gradient altitudinal : -0.6°C / 100 m
T_base <- 15 - 0.006 * altitude

# Niche thermique de l'espèce
T_opt  <- 6    # température optimale (°C)
T_sd   <- 3    # tolérance thermique (écart-type)

# Assembler les patchs
patches <- data.frame(
  id       = 1:N,
  x_coord  = coords$x,
  y_coord  = coords$y,
  altitude = altitude,
  A_base   = A_base,
  T_base   = T_base
)

cat("── Patchs générés ──\n")
print(summary(patches[, c("altitude", "A_base", "T_base")]))

# ===========================================================================
# 2. Fonctions du modèle
# ===========================================================================

# Matrice de distances inter-patchs (km)
dist_mat <- as.matrix(dist(coords))

# Noyau de dispersion exponentiel négatif
#   K(x,y) = exp(-alpha * d(x,y))
make_kernel <- function(alpha = 0.3) {
  K <- exp(-alpha * dist_mat)
  diag(K) <- 0
  K
}

# Aptitude climatique d'un patch (gaussienne centrée sur T_opt)
climate_suitability <- function(T_local, T_opt, T_sd) {
  exp(-0.5 * ((T_local - T_opt) / T_sd)^2)
}

# Surface effective sous climat modifié
#   A_eff = A_base * suitability(T_local)
# + réduction de surface en altitude (montée de la limite arborée, etc.)

# Taux d'extinction local : e_x / A_x  (rescue effect via aire)
# Taux de colonisation modulé par la capacité climatique

# ===========================================================================
# 3. Simulation : dp_x/dt (Euler)
# ===========================================================================

simulate_metapop <- function(patches, K, delta_T = 0,
                              c0 = 0.5, e0 = 0.2, alpha = 0.3,
                              T_opt = 6, T_sd = 3,
                              dt = 0.1, t_max = 200,
                              p_init = NULL) {

  N <- nrow(patches)

  # Température locale modifiée par le réchauffement
  T_local <- patches$T_base + delta_T

  # Aptitude climatique par patch
  suit <- climate_suitability(T_local, T_opt, T_sd)

  # Surface effective (réduite si le climat devient inadéquat)
  A_eff <- patches$A_base * suit
  # Seuil : en dessous de 0.1 ha effectif, le patch est inhabitable
  A_eff <- pmax(A_eff, 0.01)

  # Taux d'extinction : inversement proportionnel à la surface effective
  e_x <- e0 / A_eff

  # Taux de colonisation modulé par l'aptitude locale
  c_x <- c0 * suit

  # Occupancy initiale
  if (is.null(p_init)) {
    p <- ifelse(suit > 0.3, runif(N, 0.3, 0.8), runif(N, 0, 0.1))
  } else {
    p <- p_init
  }

  # Stockage des résultats
  n_steps <- ceiling(t_max / dt)
  save_every <- max(1, floor(n_steps / 500))
  t_save <- numeric()
  p_save <- list()
  metapop_p <- numeric()

  for (step in 1:n_steps) {

    # Taux de colonisation pour chaque patch x
    colonisation <- c_x * ((K %*% (p * A_eff)) * (1 - p))[, 1]

    # Taux d'extinction pour chaque patch x
    extinction <- e_x * p

    # Mise à jour (Euler explicite)
    dp <- colonisation - extinction
    p  <- p + dt * dp
    p  <- pmin(pmax(p, 0), 1)  # borner entre 0 et 1

    # Sauvegarder
    if (step %% save_every == 0) {
      t_save    <- c(t_save, step * dt)
      p_save    <- c(p_save, list(p))
      metapop_p <- c(metapop_p, mean(p))
    }
  }

  list(
    t       = t_save,
    p_mat   = do.call(rbind, p_save),     # matrice temps × patchs
    p_mean  = metapop_p,                  # occupancy moyenne
    patches = mutate(patches,
                     T_local  = T_local,
                     suit     = suit,
                     A_eff    = A_eff,
                     p_final  = p),
    delta_T = delta_T
  )
}

# ===========================================================================
# 4. Lancer les scénarios
# ===========================================================================

K <- make_kernel(alpha = 0.2)

cat("\n── Simulation : Baseline (ΔT = 0°C) ──\n")
res_base <- simulate_metapop(patches, K, delta_T = 0)

cat("── Simulation : Réchauffement modéré (ΔT = +2°C) ──\n")
res_mod  <- simulate_metapop(patches, K, delta_T = 2)

cat("── Simulation : Réchauffement sévère (ΔT = +4°C) ──\n")
res_sev  <- simulate_metapop(patches, K, delta_T = 4)

cat("── Simulation : Réchauffement extrême (ΔT = +6°C) ──\n")
res_ext  <- simulate_metapop(patches, K, delta_T = 6)

# ===========================================================================
# 5. Figures
# ===========================================================================

theme_metapop <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

# ── 5a. Dynamique temporelle de l'occupancy moyenne ──

df_dyn <- bind_rows(
  data.frame(t = res_base$t, p = res_base$p_mean, scenario = "Baseline (ΔT=0°C)"),
  data.frame(t = res_mod$t,  p = res_mod$p_mean,  scenario = "Modéré (ΔT=+2°C)"),
  data.frame(t = res_sev$t,  p = res_sev$p_mean,  scenario = "Sévère (ΔT=+4°C)"),
  data.frame(t = res_ext$t,  p = res_ext$p_mean,  scenario = "Extrême (ΔT=+6°C)")
)

p1 <- ggplot(df_dyn, aes(t, p, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Baseline (ΔT=0°C)"  = "#2166AC",
                                 "Modéré (ΔT=+2°C)"   = "#FDAE61",
                                 "Sévère (ΔT=+4°C)"   = "#D73027",
                                 "Extrême (ΔT=+6°C)"  = "#67001F")) +
  labs(title    = "Dynamique de la métapopulation sous changement climatique",
       subtitle = "Occupancy moyenne (p̄) au cours du temps",
       x = "Temps", y = "Occupancy moyenne (p̄)",
       color = "Scénario") +
  ylim(0, 1) +
  theme_metapop

p1

# ── 5b. Carte spatiale : occupancy finale par scénario ──

map_data <- bind_rows(
  mutate(res_base$patches, scenario = "Baseline (ΔT=0°C)"),
  mutate(res_mod$patches,  scenario = "Modéré (ΔT=+2°C)"),
  mutate(res_sev$patches,  scenario = "Sévère (ΔT=+4°C)"),
  mutate(res_ext$patches,  scenario = "Extrême (ΔT=+6°C)")
)

p2 <- ggplot(map_data, aes(x_coord, y_coord)) +
  geom_point(aes(size = A_eff, fill = p_final),
             shape = 21, alpha = 0.85, stroke = 0.3) +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF",
                       high = "#1A9850", midpoint = 0.5,
                       limits = c(0, 1)) +
  scale_size_continuous(range = c(1, 8), name = "Surface eff. (ha)") +
  facet_wrap(~scenario, ncol = 2) +
  labs(title    = "Occupancy finale par patch et scénario",
       subtitle = "Taille = surface effective | Couleur = occupancy",
       x = "X (km)", y = "Y (km)", fill = "p_x final") +
  coord_equal() +
  theme_metapop
p2

# ── 5c. Occupancy finale vs altitude ──

p3 <- ggplot(map_data, aes(altitude, p_final, color = scenario)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c("Baseline (ΔT=0°C)"  = "#2166AC",
                                 "Modéré (ΔT=+2°C)"   = "#FDAE61",
                                 "Sévère (ΔT=+4°C)"   = "#D73027",
                                 "Extrême (ΔT=+6°C)"  = "#67001F")) +
  labs(title    = "Occupancy finale selon l'altitude",
       subtitle = "Déplacement altitudinal de la zone de persistance",
       x = "Altitude (m)", y = "Occupancy finale (p_x)",
       color = "Scénario") +
  ylim(0, 1) +
  theme_metapop

p3

# ── 5d. Aptitude climatique vs altitude ──

p4 <- ggplot(map_data, aes(altitude, suit, color = scenario)) +
  geom_line(stat = "smooth", method = "loess", linewidth = 1.2, se = FALSE) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("Baseline (ΔT=0°C)"  = "#2166AC",
                                 "Modéré (ΔT=+2°C)"   = "#FDAE61",
                                 "Sévère (ΔT=+4°C)"   = "#D73027",
                                 "Extrême (ΔT=+6°C)"  = "#67001F")) +
  labs(title    = "Aptitude climatique selon l'altitude",
       subtitle = "Déplacement de la niche vers les sommets",
       x = "Altitude (m)", y = "Aptitude climatique",
       color = "Scénario") +
  theme_metapop

p4

# ── 5e. Analyse de sensibilité : capacité de dispersion ──

cat("\n── Analyse de sensibilité : alpha (dispersion) ──\n")
alpha_vals <- c(0.1, 0.2, 0.3, 0.5, 0.8)
sens_results <- list()

for (a in alpha_vals) {
  K_a <- make_kernel(alpha = a)
  res <- simulate_metapop(patches, K_a, delta_T = 4)
  sens_results[[as.character(a)]] <- data.frame(
    t     = res$t,
    p     = res$p_mean,
    alpha = paste0("α = ", a)
  )
}

df_sens <- bind_rows(sens_results)

p5 <- ggplot(df_sens, aes(t, p, color = alpha)) +
  geom_line(linewidth = 1) +
  labs(title    = "Sensibilité à la capacité de dispersion (ΔT = +4°C)",
       subtitle = "α petit = grande distance de dispersion",
       x = "Temps", y = "Occupancy moyenne (p̄)",
       color = "Paramètre") +
  ylim(0, 1) +
  theme_metapop

p5
# ── 5f. Seuil d'extinction : balayage de ΔT ──

cat("── Balayage du seuil de réchauffement ──\n")
delta_T_seq <- seq(0, 8, by = 0.5)
p_equil <- numeric(length(delta_T_seq))

for (i in seq_along(delta_T_seq)) {
  res_i <- simulate_metapop(patches, K, delta_T = delta_T_seq[i], t_max = 300)
  p_equil[i] <- tail(res_i$p_mean, 1)
}

df_thresh <- data.frame(delta_T = delta_T_seq, p_equil = p_equil)

p6 <- ggplot(df_thresh, aes(delta_T, p_equil)) +
  geom_line(linewidth = 1.2, color = "#D73027") +
  geom_point(size = 2.5, color = "#67001F") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 6, y = 0.15, label = "Seuil quasi-extinction",
           color = "grey40", size = 3.5) +
  labs(title    = "Seuil d'extinction de la métapopulation",
       subtitle = "Occupancy à l'équilibre en fonction du réchauffement",
       x = "Réchauffement ΔT (°C)", y = "Occupancy à l'équilibre (p̄)") +
  ylim(0, 1) +
  theme_metapop

p6

# ===========================================================================
# 6. Composer et sauvegarder
# ===========================================================================

# Figure composite
fig_main <- (p1 | p6) / (p3 | p4) +
  plot_annotation(
    title    = "Métapopulation en paysage montagneux sous changement climatique",
    subtitle = expression(dp[x]/dt == c[x] * Sigma * K(x,y) * p[y] * A[y] * (1-p[x]) - e[x]/A[x] * p[x]),
    theme    = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("metapop_climate_composite.png", fig_main,
       width = 14, height = 10, dpi = 200, bg = "white")
cat("✔ Figure composite sauvegardée : metapop_climate_composite.png\n")

ggsave("metapop_carte_scenarios.png", p2,
       width = 12, height = 10, dpi = 200, bg = "white")
cat("✔ Carte spatiale sauvegardée : metapop_carte_scenarios.png\n")

ggsave("metapop_sensibilite_dispersion.png", p5,
       width = 10, height = 6, dpi = 200, bg = "white")
cat("✔ Analyse de sensibilité sauvegardée : metapop_sensibilite_dispersion.png\n")

# ===========================================================================
# 7. Résumé numérique
# ===========================================================================

cat("\n══════════════════════════════════════════════════════\n")
cat("  RÉSUMÉ : OCCUPANCY À L'ÉQUILIBRE PAR SCÉNARIO\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  Baseline (ΔT=0°C)  : p̄ = %.3f\n", tail(res_base$p_mean, 1)))
cat(sprintf("  Modéré   (ΔT=+2°C) : p̄ = %.3f\n", tail(res_mod$p_mean, 1)))
cat(sprintf("  Sévère   (ΔT=+4°C) : p̄ = %.3f\n", tail(res_sev$p_mean, 1)))
cat(sprintf("  Extrême  (ΔT=+6°C) : p̄ = %.3f\n", tail(res_ext$p_mean, 1)))
cat("══════════════════════════════════════════════════════\n")

# Nombre de patchs quasi-éteints (p < 0.05) par scénario
cat("\n  Patchs quasi-éteints (p < 0.05) :\n")
cat(sprintf("  Baseline : %d / %d\n", sum(res_base$patches$p_final < 0.05), N))
cat(sprintf("  Modéré   : %d / %d\n", sum(res_mod$patches$p_final  < 0.05), N))
cat(sprintf("  Sévère   : %d / %d\n", sum(res_sev$patches$p_final  < 0.05), N))
cat(sprintf("  Extrême  : %d / %d\n", sum(res_ext$patches$p_final  < 0.05), N))
cat("══════════════════════════════════════════════════════\n")

cat("\n✔ Script terminé. Bonne exploration !\n")
