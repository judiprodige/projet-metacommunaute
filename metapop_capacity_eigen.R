###############################################################################
#  Analyse spectrale de la capacité de métapopulation                         #
#                                                                             #
#  Capacité (Hanski & Ovaskainen) : lambda_M = plus grande valeur propre de   #
#      M_xy = (A_x * s_x^1.5) * K(x,y) * (A_y * s_y^1.5)                      #
#  Persistance <=> lambda_M > e0/c0                                           #
#                                                                             #
#  Le vecteur propre associé w donne la contribution de chaque patch :        #
#      contribution_x = w_x^2   (somme = 1)                                   #
#                                                                             #
#  Prérequis : metapop_climate_mountain.R (paysage, K, paramètres)            #
###############################################################################

source("metapop_climate_mountain.R")

# ===========================================================================
# 1. Capacité de métapopulation et décomposition spectrale
# ===========================================================================

c0 <- 0.5; e0 <- 0.2

# Aptitude climatique et surface effective sous un réchauffement donné
suitability_dT <- function(delta_T) {
  climate_suitability(patches$T_base + delta_T, T_opt = 6, T_sd = 3)
}
A_eff_dT <- function(delta_T) pmax(patches$A_base * suitability_dT(delta_T), 0.01)

# Matrice de capacité : symétrique, donc valeurs propres réelles
capacity_matrix <- function(delta_T) {
  s <- suitability_dT(delta_T)
  v <- sqrt(s) * A_eff_dT(delta_T)   # = A_x * s_x^1.5 (hors plancher)
  outer(v, v) * K
}

# Décomposition : lambda_M, vecteur propre normé, contribution par patch
capacity_eigen <- function(delta_T) {
  M <- capacity_matrix(delta_T)
  e <- eigen(M, symmetric = TRUE)
  w <- abs(e$vectors[, 1])           # Perron-Frobenius : w > 0
  w <- w / sqrt(sum(w^2))
  list(lambda   = e$values[1],
       lambda2  = e$values[2],
       w        = w,
       contrib  = w^2,               # somme = 1
       persist  = e$values[1] > e0 / c0)
}

# Équilibre analytique : p* = S/(S + delta), point fixe unique
equilibrium <- function(delta_T, tol = 1e-12, iter = 1e4) {
  s <- suitability_dT(delta_T); Ae <- A_eff_dT(delta_T)
  delta <- e0 / (c0 * s * Ae)
  p <- rep(0.5, nrow(patches))
  for (k in seq_len(iter)) {
    p_new <- as.numeric(K %*% (p * Ae)); p_new <- p_new / (p_new + delta)
    if (max(abs(p_new - p)) < tol) return(p_new)
    p <- p_new
  }
  p
}

# ===========================================================================
# 2. Résumé numérique par scénario
# ===========================================================================

scenarios <- c(0, 2, 4, 6)
lab_dT <- function(dT) ifelse(dT == 0, "ΔT = 0 °C", sprintf("ΔT = +%d °C", dT))

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  CAPACITÉ DE MÉTAPOPULATION ET CONTRIBUTION DES PATCHS\n")
cat("  seuil de persistance : lambda_M > e0/c0 =", e0 / c0, "\n")
cat("══════════════════════════════════════════════════════════════\n")

for (dT in scenarios) {
  E <- capacity_eigen(dT); ctr <- E$contrib
  cat(sprintf("\n%s\n", lab_dT(dT)))
  cat(sprintf("  lambda_1 = %8.3f   lambda_2 = %7.3f   lambda_2/lambda_1 = %.2f\n",
              E$lambda, E$lambda2, E$lambda2 / E$lambda))
  cat(sprintf("  persistance : %s   (marge = %.0f x le seuil)\n",
              ifelse(E$persist, "OUI", "NON"), E$lambda / (e0 / c0)))
  cat(sprintf("  concentration : top 5 = %.0f %% de la capacité, top 10 = %.0f %%\n",
              100 * sum(sort(ctr, decreasing = TRUE)[1:5]),
              100 * sum(sort(ctr, decreasing = TRUE)[1:10])))
  cat(sprintf("  altitude moyenne pondérée par w² : %.0f m (paysage : %.0f m)\n",
              sum(ctr * patches$altitude), mean(patches$altitude)))
  top <- order(ctr, decreasing = TRUE)[1:3]
  cat("  patchs pivots :",
      paste(sprintf("#%d (%.0f m, %.0f %%)", top, patches$altitude[top], 100 * ctr[top]),
            collapse = " · "), "\n")
}

# Sensibilité au climat décomposée sur le vecteur propre :
#   d ln(lambda_M)/d(ΔT) = 3 * sum_x w_x^2 * d ln(s_x)/d(ΔT)
cat("\n  Sensibilité au réchauffement (élasticité) :\n")
for (dT in scenarios) {
  E <- capacity_eigen(dT)
  dln_s <- -(patches$T_base + dT - 6) / 3^2
  cat(sprintf("  %s : d ln(lambda_M)/dΔT = %+.3f par °C\n",
              lab_dT(dT), sum(3 * E$contrib * dln_s)))
}

# Validation : retrait effectif du patch le plus contributeur
cat("\n  Retrait du patch pivot (ΔT = +4 °C) :\n")
E4 <- capacity_eigen(4); M4 <- capacity_matrix(4)
piv <- which.max(E4$contrib)
lam_sans <- max(eigen(M4[-piv, -piv], symmetric = TRUE, only.values = TRUE)$values)
cat(sprintf("  patch #%d : lambda_M %.2f -> %.2f  (perte réelle %.0f %%, w² = %.0f %%)\n",
            piv, E4$lambda, lam_sans, 100 * (E4$lambda - lam_sans) / E4$lambda,
            100 * E4$contrib[piv]))
cat("══════════════════════════════════════════════════════════════\n")

# ===========================================================================
# 3. Figure : contribution de chaque patch à la capacité
# ===========================================================================

df <- bind_rows(lapply(scenarios, function(dT) {
  E <- capacity_eigen(dT)
  data.frame(
    id       = patches$id,
    x_coord  = patches$x_coord,
    y_coord  = patches$y_coord,
    altitude = patches$altitude,
    A_eff    = A_eff_dT(dT),
    suit     = suitability_dT(dT),
    contrib  = E$contrib,
    p_star   = equilibrium(dT),
    lambda   = E$lambda,
    scenario = factor(lab_dT(dT), levels = lab_dT(scenarios))
  )
}))

lab_lambda <- df |>
  group_by(scenario) |>
  summarise(lambda = first(lambda), .groups = "drop") |>
  mutate(txt = sprintf("λ_M = %.1f", lambda))

theme_cap <- theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 10),
        strip.text    = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# ── (a) carte : où se trouve le cœur de la capacité ──
pa <- ggplot(df, aes(x_coord, y_coord)) +
  geom_point(aes(size = A_eff, fill = 100 * contrib),
             shape = 21, stroke = 0.25, colour = "white") +
  geom_text(data = lab_lambda, aes(x = 0.5, y = 29.5, label = txt),
            hjust = 0, vjust = 1, size = 3.2, colour = "grey30", inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Contribution w² (%)") +
  scale_size_continuous(range = c(1, 7), name = "A_eff (ha)") +
  facet_wrap(~scenario, nrow = 1) +
  coord_equal(xlim = c(0, 30), ylim = c(0, 30)) +
  labs(title    = "Où se concentre la capacité de métapopulation",
       subtitle = "Contribution w² de chaque patch à λ_M · taille = surface effective",
       x = "X (km)", y = "Y (km)") +
  theme_cap
pa

# ── (b) contribution selon l'altitude ──
pb <- ggplot(df, aes(altitude, 100 * contrib, colour = scenario)) +
  geom_point(aes(size = A_eff), alpha = 0.75) +
  scale_colour_manual(values = c("#2166AC", "#7A9E4F", "#D9822B", "#8A3324"),
                      name = NULL) +
  scale_size_continuous(range = c(0.8, 4), guide = "none") +
  labs(title    = "La capacité remonte le versant",
       subtitle = "Les patchs pivots ne sont pas les mêmes selon le scénario",
       x = "Altitude (m)", y = "Contribution à λ_M (%)") +
  theme_cap
pb

# ── (c) concentration : courbe cumulée par rang ──
df_rank <- df |>
  group_by(scenario) |>
  arrange(desc(contrib), .by_group = TRUE) |>
  mutate(rang = row_number(), cum = 100 * cumsum(contrib)) |>
  ungroup()

pc <- ggplot(df_rank, aes(rang, cum, colour = scenario)) +
  geom_abline(slope = 100 / 40, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.9) +
  annotate("text", x = 27, y = 60, label = "contribution\nuniforme",
           colour = "grey50", size = 3, hjust = 0) +
  scale_colour_manual(values = c("#2166AC", "#7A9E4F", "#D9822B", "#8A3324"),
                      name = NULL) +
  labs(title    = "Quelques patchs portent tout le réseau",
       subtitle = "Contribution cumulée des patchs classés par importance",
       x = "Patchs classés par contribution décroissante", y = "Capacité cumulée (%)") +
  theme_cap
pc

fig <- pa / (pb | pc) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title = "Analyse spectrale : quels patchs font la capacité de métapopulation ?",
    subtitle = "λ_M = plus grande valeur propre de M ; w = vecteur propre associé ; contribution = w²",
    theme = theme(plot.title    = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(colour = "grey40")))

ggsave("metapop_contribution_patchs.png", fig,
       width = 13, height = 9.5, dpi = 200, bg = "white", type = "cairo")
cat("✔ Figure sauvegardée : metapop_contribution_patchs.png\n")
