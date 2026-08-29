###############################################################################
#  Vélocité climatique vs. vélocité de colonisation                           #
#                                                                             #
#  Le climat impose une vitesse : l'optimum thermique remonte de              #
#      v_clim = r / 0.006 = 167 m par °C de réchauffement                     #
#  La métapopulation répond avec un temps de relaxation tau, d'où une         #
#  vitesse de suivi :                                                         #
#      v_col ≈ 1 / (0.006 * tau)   mètres par unité de temps (pour 1 °C)      #
#                                                                             #
#  Le retard accumulé sous réchauffement progressif vaut ≈ r * tau degrés.    #
#  Il est toujours NÉGATIF ici : la métapopulation observée est plus riche    #
#  que son équilibre — c'est une dette d'extinction, pas un déficit de        #
#  colonisation.                                                              #
#                                                                             #
#  Prérequis : metapop_climate_mountain.R                                     #
###############################################################################

source("metapop_climate_mountain.R")

c0 <- 0.5; e0 <- 0.2
lapse <- 0.006          # °C par mètre (gradient altitudinal du script)
alt   <- patches$altitude

suitability_dT <- function(delta_T) {
  climate_suitability(patches$T_base + delta_T, T_opt = 6, T_sd = 3)
}
A_eff_dT <- function(delta_T) pmax(patches$A_base * suitability_dT(delta_T), 0.01)

# ===========================================================================
# 1. Équilibre, relaxation, vélocités
# ===========================================================================

equilibrium <- function(delta_T, tol = 1e-13, iter = 2e4) {
  s <- suitability_dT(delta_T); Ae <- A_eff_dT(delta_T)
  delta <- e0 / (c0 * s * Ae)
  p <- rep(0.5, length(alt))
  for (k in seq_len(iter)) {
    q <- as.numeric(K %*% (p * Ae)); q <- q / (q + delta)
    if (max(abs(q - p)) < tol) return(q)
    p <- q
  }
  p
}

# Temps de relaxation : inverse de la valeur propre dominante de la jacobienne
#   dp_i/dt = c_i (sum_j K_ij A_j p_j)(1 - p_i) - e_i p_i
relaxation_time <- function(delta_T) {
  p <- equilibrium(delta_T)
  s <- suitability_dT(delta_T); Ae <- A_eff_dT(delta_T)
  C <- as.numeric(K %*% (p * Ae))
  J <- outer(c0 * s * (1 - p), rep(1, length(p))) * t(t(K) * Ae)
  diag(J) <- -c0 * s * C - e0 / Ae
  -1 / max(Re(eigen(J, only.values = TRUE)$values))
}

# Vitesse de suivi : mètres d'altitude par unité de temps, pour 1 °C de marge
colonisation_velocity <- function(delta_T) 1 / (lapse * relaxation_time(delta_T))

# Vitesse imposée par le climat, pour un réchauffement de r °C / unité de temps
climate_velocity <- function(r) r / lapse

# ===========================================================================
# 2. Simulation transitoire : ΔT(t) = r * t
# ===========================================================================

simulate_transient <- function(r, delta_T_max = 6, dt = 0.02, p_init = NULL) {
  p <- if (is.null(p_init)) equilibrium(0) else p_init
  t_cur <- 0; save_every <- max(1, floor((delta_T_max / r / dt) / 400))
  tr <- list()
  step <- 0
  while (r * t_cur < delta_T_max) {
    dT <- r * t_cur
    s <- suitability_dT(dT); Ae <- A_eff_dT(dT)
    dp <- c0 * s * (as.numeric(K %*% (p * Ae)) * (1 - p)) - (e0 / Ae) * p
    p <- pmin(pmax(p + dt * dp, 0), 1)
    t_cur <- t_cur + dt; step <- step + 1
    if (step %% save_every == 0) {
      tr[[length(tr) + 1]] <- data.frame(
        t = t_cur, delta_T = dT, p_bar = mean(p),
        centroid = if (sum(p) > 1e-9) sum(p * alt) / sum(p) else NA_real_)
    }
  }
  list(p = p, delta_T = r * t_cur, r = r, trace = bind_rows(tr))
}

# ===========================================================================
# 3. Résumé numérique
# ===========================================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  VÉLOCITÉ CLIMATIQUE vs VÉLOCITÉ DE COLONISATION\n")
cat("  v_clim = 167 m par °C · v_col = 1 / (0,006 · tau)\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("%-6s %-8s %-8s %-11s %-11s %-10s\n",
            "ΔT", "p̄*", "tau", "v_col", "centroïde", "optimum"))
for (dT in c(0, 2, 4, 6, 7, 8, 8.3)) {
  p <- equilibrium(dT); tau <- relaxation_time(dT)
  cen <- if (sum(p) > 1e-9) sum(p * alt) / sum(p) else NA
  cat(sprintf("%-6.1f %-8.3f %-8.2f %-11.0f %-11.0f %-10.0f\n",
              dT, mean(p), tau, 1 / (lapse * tau), cen, (15 + dT - 6) / lapse))
}
cat("\n  Le ralentissement critique : tau est multiplié par ~38 entre ΔT = 0\n")
cat("  et le seuil. La métapopulation ralentit quand il faudrait accélérer.\n")

# Retard accumulé sous réchauffement progressif
grid_dT <- seq(0.5, 9, by = 0.02)                       # branche décroissante
grid_p  <- sapply(grid_dT, function(d) mean(equilibrium(d)))
lag_in_degrees <- function(p_bar) approx(grid_p, grid_dT, xout = p_bar)$y

rates <- c(0.005, 0.02, 0.05, 0.1, 0.3, 1)
cat("\n  Retard sous réchauffement progressif (ΔT final = +6 °C) :\n")
cat(sprintf("  %-9s %-13s %-9s %-11s %-9s\n",
            "r (°C/ut)", "v_clim (m/ut)", "p̄", "ΔT équiv.", "retard"))
debt <- lapply(rates, function(r) {
  R <- simulate_transient(r, 6); pb <- mean(R$p); de <- lag_in_degrees(pb)
  cat(sprintf("  %-9.3f %-13.1f %-9.3f %-11.2f %+.2f °C\n",
              r, climate_velocity(r), pb, de, de - 6))
  data.frame(r = r, v_clim = climate_velocity(r), p_bar = pb, lag = de - 6)
})
debt <- bind_rows(debt)
cat("\n  Retard toujours négatif = DETTE D'EXTINCTION : l'occupation observée\n")
cat("  surestime la viabilité. Près du seuil, d'un facteur 5.\n")
cat("══════════════════════════════════════════════════════════════\n")

# ===========================================================================
# 4. Figure
# ===========================================================================

dT_grid <- seq(0, 8.3, by = 0.1)
df_vel <- data.frame(
  delta_T = dT_grid,
  v_col   = sapply(dT_grid, colonisation_velocity)
)
rates_shown <- c(0.01, 0.05, 0.2)
df_clim <- data.frame(r = rates_shown, v = climate_velocity(rates_shown),
                      lab = sprintf("r = %.2f °C/ut", rates_shown))

theme_vel <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey40", size = 10),
        panel.grid.minor = element_blank(), legend.position = "bottom")

# ── (a) les deux vélocités se croisent ──
pa <- ggplot(df_vel, aes(delta_T, v_col)) +
  geom_hline(data = df_clim, aes(yintercept = v), colour = "#B4552F",
             linetype = "dashed", linewidth = 0.6) +
  geom_text(data = df_clim, aes(x = 0.1, y = v * 1.25, label = lab),
            inherit.aes = FALSE, hjust = 0, size = 3, colour = "#B4552F") +
  geom_line(linewidth = 1.2, colour = "#2C6E8F") +
  scale_y_log10() +
  labs(title    = "Les deux vélocités se croisent",
       subtitle = "v_col chute avec le réchauffement ; v_clim reste constante",
       x = "Réchauffement ΔT (°C)",
       y = "Vitesse (m d'altitude par unité de temps, échelle log)") +
  theme_vel

# ── (b) ralentissement critique ──
pb <- ggplot(data.frame(delta_T = dT_grid,
                        tau = sapply(dT_grid, relaxation_time)),
             aes(delta_T, tau)) +
  geom_line(linewidth = 1.2, colour = "#8A3324") +
  labs(title    = "Ralentissement critique",
       subtitle = "Temps de relaxation à l'approche du seuil (ΔT* = 8,4 °C)",
       x = "Réchauffement ΔT (°C)", y = "Temps de relaxation τ") +
  theme_vel

# ── (c) dette d'extinction : transitoire vs équilibre ──
df_eq <- data.frame(delta_T = grid_dT, p = grid_p, type = "Équilibre")
df_tr <- bind_rows(lapply(c(0.05, 0.3, 1), function(r) {
  tr <- simulate_transient(r, 6)$trace
  data.frame(delta_T = tr$delta_T, p = tr$p_bar,
             type = sprintf("Transitoire r = %.2f", r))
}))

pc <- ggplot(bind_rows(df_eq, df_tr), aes(delta_T, p, colour = type)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Équilibre" = "grey30",
                                 "Transitoire r = 0.05" = "#7FAECC",
                                 "Transitoire r = 0.30" = "#3874A0",
                                 "Transitoire r = 1.00" = "#0A304C"),
                      name = NULL) +
  labs(title    = "La dette d'extinction",
       subtitle = "Sous réchauffement progressif, l'occupation reste au-dessus de son équilibre",
       x = "Réchauffement ΔT (°C)", y = "Occupancy moyenne (p̄)") +
  theme_vel

# ── (d) l'occupation ne rejoint jamais son optimum ──
df_track <- data.frame(
  delta_T  = dT_grid,
  optimum  = (15 + dT_grid - 6) / lapse,
  occupe   = sapply(dT_grid, function(d) {
    p <- equilibrium(d); if (sum(p) > 1e-9) sum(p * alt) / sum(p) else NA }))

pd <- ggplot(df_track, aes(delta_T)) +
  geom_ribbon(aes(ymin = occupe, ymax = pmin(optimum, 3000)),
              fill = "#B4552F", alpha = 0.12) +
  geom_line(aes(y = optimum), linewidth = 1, colour = "#B4552F") +
  geom_line(aes(y = occupe),  linewidth = 1.2, colour = "#2C6E8F") +
  geom_hline(yintercept = max(alt), linetype = "dotted", colour = "grey40") +
  annotate("text", x = 0.2, y = max(alt) + 60, hjust = 0, size = 3,
           colour = "grey40", label = sprintf("sommet du massif (%.0f m)", max(alt))) +
  annotate("text", x = 6.2, y = 2750, hjust = 1, size = 3, colour = "#B4552F",
           label = "optimum thermique") +
  annotate("text", x = 6.2, y = 1700, hjust = 1, size = 3, colour = "#2C6E8F",
           label = "centroïde occupé") +
  labs(title    = "L'écart n'est pas un retard, c'est une impasse",
       subtitle = "L'optimum sort du massif ; l'occupation plafonne sous le sommet",
       x = "Réchauffement ΔT (°C)", y = "Altitude (m)") +
  theme_vel

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title    = "Vélocité climatique vs. vélocité de colonisation",
    subtitle = "v_clim = 167 m par °C · v_col = 1/(0,006·τ) · retard ≈ r·τ",
    theme    = theme(plot.title = element_text(face = "bold", size = 15),
                     plot.subtitle = element_text(colour = "grey40")))

ggsave("metapop_velocite_climatique.png", fig,
       width = 13, height = 9.5, dpi = 200, bg = "white")
cat("✔ Figure sauvegardée : metapop_velocite_climatique.png\n")
