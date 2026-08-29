###############################################################################
#  Métacommunauté trophique : une chaîne alimentaire par patch                #
#                                                                             #
#  Contrainte trophique (théorie trophique de la biogéographie insulaire,     #
#  Gravel et al. 2011 ; cf. Cameron 2022, ressources/) :                      #
#                                                                             #
#    • le niveau 1 (ressource basale) suit le SPOM classique                  #
#    • le niveau i > 1 ne peut coloniser un patch QUE si le niveau i−1 y est  #
#      déjà présent, et il DISPARAÎT avec lui                                 #
#                                                                             #
#    dp_i,x/dt = c·s_x·(Σ_y K_xy A_y p_i,y)·(p_{i−1,x} − p_i,x)              #
#                − (i · e / A_x^eff) · p_i,x                                  #
#                                                                             #
#  Le facteur i sur l'extinction vient de la cascade : le niveau i s'éteint   #
#  de lui-même ET avec chacun des niveaux dont il dépend.                     #
#                                                                             #
#  DÉMOGRAPHIE ET NICHE IDENTIQUES POUR LES TROIS NIVEAUX : toute différence  #
#  observée est donc l'effet de la seule contrainte trophique.                #
#                                                                             #
#  Prérequis : metapop_climate_mountain.R                                     #
###############################################################################

source("metapop_climate_mountain.R")

c0 <- 0.5; e0 <- 0.2; L <- 3          # trois niveaux : plante, herbivore, prédateur
niveaux <- c("1 · Ressource basale", "2 · Consommateur", "3 · Prédateur")
alt <- patches$altitude

suitability_dT <- function(delta_T) {
  climate_suitability(patches$T_base + delta_T, T_opt = 6, T_sd = 3)
}
A_eff_dT <- function(delta_T) pmax(patches$A_base * suitability_dT(delta_T), 0.01)

# ===========================================================================
# 1. Équilibre de la chaîne — résolu niveau par niveau
# ===========================================================================
# Le niveau i ne dépend que de i−1 : la chaîne se résout du bas vers le haut.
#   p*_i,x = c·s_x·S_i,x·p_{i−1,x} / (c·s_x·S_i,x + i·e/A_x)

equilibrium_chain <- function(delta_T, tol = 1e-12, iter = 2e4) {
  s <- suitability_dT(delta_T); Ae <- A_eff_dT(delta_T)
  P <- matrix(0, nrow = length(alt), ncol = L)
  support <- rep(1, length(alt))            # niveau 0 : le patch lui-même
  for (i in seq_len(L)) {
    loss <- i * e0 / Ae                     # cascade d'extinction
    p <- support * 0.5
    for (k in seq_len(iter)) {
      col <- c0 * s * as.numeric(K %*% (p * Ae))
      q <- col * support / (col + loss)
      if (max(abs(q - p)) < tol) { p <- q; break }
      p <- q
    }
    P[, i] <- p; support <- p
  }
  P
}

# Capacité de métapopulation du niveau i, sachant le niveau inférieur.
#   persistance <=> lambda > 1, avec v_x = sqrt(c·s_x·support_x / loss_x) · A_x^eff
capacity_level <- function(delta_T, i, support) {
  s <- suitability_dT(delta_T); Ae <- A_eff_dT(delta_T)
  loss <- i * e0 / Ae
  v <- sqrt(c0 * s * pmax(support, 1e-12) / loss) * Ae
  max(eigen(outer(v, v) * K, symmetric = TRUE, only.values = TRUE)$values)
}

capacities <- function(delta_T) {
  P <- equilibrium_chain(delta_T)
  support <- cbind(1, P[, -L, drop = FALSE])
  sapply(seq_len(L), function(i) capacity_level(delta_T, i, support[, i]))
}

# Seuil d'extinction de chaque niveau
threshold_level <- function(i, interval = c(0, 9)) {
  f <- function(d) capacities(d)[i] - 1
  if (f(interval[1]) < 0) return(NA_real_)
  if (f(interval[2]) > 0) return(NA_real_)
  uniroot(f, interval, tol = 1e-3)$root
}

# ===========================================================================
# 2. Résumé numérique
# ===========================================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  MÉTACOMMUNAUTÉ TROPHIQUE — CHAÎNE À", L, "NIVEAUX\n")
cat("  démographie et niche thermique identiques aux trois niveaux\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("%-7s %-24s %-24s %-14s\n", "ΔT",
            "occupancy moyenne", "patchs occupés / 30", "longueur chaîne"))
for (dT in c(0, 2, 4, 6)) {
  P <- equilibrium_chain(dT)
  cat(sprintf("%-7s %-24s %-24s %-14.2f\n",
      sprintf("+%d °C", dT),
      paste(sprintf("%.3f", colMeans(P)), collapse = " / "),
      paste(sprintf("%d", colSums(P > 0.05)), collapse = " / "),
      mean(rowSums(P))))
}

cat("\n  Seuils d'extinction par niveau trophique :\n")
seuils <- sapply(seq_len(L), threshold_level)
for (i in seq_len(L)) {
  cat(sprintf("  %-22s ΔT* = %s\n", niveaux[i],
              ifelse(is.na(seuils[i]), "hors de [0 ; 9]", sprintf("%.2f °C", seuils[i]))))
}
cat(sprintf("\n  Le prédateur disparaît %.1f °C avant sa ressource :\n",
            seuils[1] - seuils[L]))
cat("  la contrainte trophique, à démographie égale, suffit à décaler le seuil.\n")

cat("\n  Longueur de chaîne selon l'altitude (ΔT = +4 °C) :\n")
P4 <- equilibrium_chain(4); len4 <- rowSums(P4)
q <- cut(alt, breaks = c(0, 1300, 1700, 2100, Inf),
         labels = c("< 1300 m", "1300–1700 m", "1700–2100 m", "> 2100 m"))
for (lv in levels(q)) if (any(q == lv))
  cat(sprintf("  %-14s %.2f niveaux (n = %d patchs)\n", lv, mean(len4[q == lv]), sum(q == lv)))
cat("══════════════════════════════════════════════════════════════\n")

# ===========================================================================
# 3. Figure
# ===========================================================================

scen <- c(0, 2, 4, 6)
dT_grid <- seq(0, 9, by = 0.15)

df_occ <- bind_rows(lapply(dT_grid, function(d) {
  P <- equilibrium_chain(d)
  data.frame(delta_T = d, niveau = niveaux, p = colMeans(P))
}))

df_cap <- bind_rows(lapply(dT_grid, function(d)
  data.frame(delta_T = d, niveau = niveaux, lambda = capacities(d))))

df_patch <- bind_rows(lapply(scen, function(d) {
  P <- equilibrium_chain(d)
  data.frame(x_coord = patches$x_coord, y_coord = patches$y_coord,
             altitude = alt, A_eff = A_eff_dT(d), longueur = rowSums(P),
             scenario = factor(ifelse(d == 0, "ΔT = 0 °C", sprintf("ΔT = +%d °C", d)),
                               levels = c("ΔT = 0 °C", sprintf("ΔT = +%d °C", scen[-1]))))
}))

pal <- c("#3A8167", "#C4693F", "#2C6E8F")
names(pal) <- niveaux

theme_tro <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey40", size = 10),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "bottom")

# ── (a) occupancy par niveau ──
pa <- ggplot(df_occ, aes(delta_T, p, colour = niveau)) +
  geom_vline(xintercept = seuils[!is.na(seuils)], linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(title    = "Le sommet de la chaîne cède en premier",
       subtitle = "Occupancy à l'équilibre par niveau trophique ; pointillés = seuils",
       x = "Réchauffement ΔT (°C)", y = "Occupancy moyenne (p̄)") +
  theme_tro
pa
# ── (b) capacité par niveau ──
pb <- ggplot(df_cap, aes(delta_T, lambda, colour = niveau)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1.1) +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.2, 400)) +          # sous le seuil, λ plonge vers 0
  annotate("text", x = 0.1, y = 1.45, hjust = 0, size = 3, colour = "grey40",
           label = "seuil de persistance λ = 1") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(title    = "Capacité de métapopulation par niveau",
       subtitle = "Chaque niveau hérite de la fragilité de tous ceux du dessous",
       x = "Réchauffement ΔT (°C)", y = "λ (échelle log)") +
  theme_tro
pb
# ── (c) longueur de chaîne dans l'espace ──
pc <- ggplot(df_patch, aes(x_coord, y_coord)) +
  geom_point(aes(size = A_eff, fill = longueur), shape = 21,
             stroke = 0.25, colour = "white") +
  scale_fill_gradientn(colours = c("#EDEDE8", "#9CC4B4", "#3A8167", "#14402F"),
                       limits = c(0, 3), name = "Longueur de chaîne") +
  scale_size_continuous(range = c(1, 6), name = "A_eff (ha)") +
  facet_wrap(~scenario, nrow = 1) + coord_equal() +
  labs(title    = "La chaîne raccourcit avant de disparaître",
       subtitle = "Nombre attendu de niveaux présents par patch",
       x = "X (km)", y = "Y (km)") +
  theme_tro
pc
# ── (d) longueur de chaîne selon l'altitude ──
pd <- ggplot(df_patch, aes(altitude, longueur, colour = scenario)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_colour_manual(values = c("#2C6E8F", "#7A9E4F", "#D9822B", "#8A3324"), name = NULL) +
  labs(title    = "Les patchs marginaux perdent leurs niveaux supérieurs",
       subtitle = "Longueur de chaîne selon l'altitude",
       x = "Altitude (m)", y = "Nombre de niveaux présents") +
  theme_tro

pd

fig <- (pa | pb) / pc / pd +
  plot_layout(heights = c(1, 0.95, 0.8)) +
  plot_annotation(
    title    = "Métacommunauté trophique : une chaîne alimentaire par patch",
    subtitle = "Le niveau i ne colonise que si i−1 est présent, et s'éteint avec lui",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(colour = "grey40")))

ggsave("metapop_reseau_trophique.png", fig,
       width = 13, height = 13, dpi = 200, bg = "white")
cat("✔ Figure sauvegardée : metapop_reseau_trophique.png\n")
