# projet-metacommunaute

Projet test dans le cadre d'une école d'été sur l'écologie spatiale.

Ce dépôt contient un **modèle de métapopulation spatialement explicite** (SPOM,
*Stochastic Patch Occupancy Model* de Hanski) appliqué à un paysage **montagneux**
soumis au **changement climatique**, ainsi qu'un document pédagogique expliquant
le modèle.

## Contenu du dépôt

| Fichier | Description |
|---|---|
| [metapop_climate_mountain.R](metapop_climate_mountain.R) | Script R complet : génération du paysage, simulation, figures, résumé numérique |
| [Modèle de métapopulation spatialement explicite.html](Modèle%20de%20métapopulation%20spatialement%20explicite.html) | Document explicatif (théorie + interprétation des résultats + code commenté) |
| [Modèle de métapopulation spatialement explicite.pdf](Modèle%20de%20métapopulation%20spatialement%20explicite.pdf) | Même document, version PDF |
| [ressources/](ressources/) | Littérature de référence (`cameron_victor_MSc_2022.pdf`) |

## Le modèle

Équation de Hanski en temps continu, intégrée par Euler explicite :

```
dp_x/dt = c_x * Σ_{y≠x} K(x,y) · p_y · A_y · (1 − p_x)  −  (e_x / A_x) · p_x
```

- `p_x` : probabilité d'occupation du patch *x*
- `K(x,y) = exp(−α · d(x,y))` : noyau de dispersion exponentiel négatif
- `A_x` : surface (effective) du patch — plus grand = extinction plus faible (*rescue effect*)
- `c_x`, `e_x` : taux de colonisation et d'extinction

### Couplage au climat

1. **Gradient altitudinal** : `T = 15 − 0.006 · altitude` (−0,6 °C / 100 m).
2. **Niche thermique gaussienne** : aptitude `suit = exp(−0.5·((T − T_opt)/T_sd)²)`,
   avec `T_opt = 6 °C` et `T_sd = 3 °C`.
3. **Modulation des paramètres** : `A_eff = A_base · suit`, `c_x = c0 · suit`,
   `e_x = e0 / A_eff`. Un réchauffement `ΔT` décale la température locale et pousse
   la zone favorable vers les sommets (*piège sommital* : moins de surface, moins
   de connectivité en altitude).

### Paysage simulé

40 patchs tirés aléatoirement (`set.seed(42)`) sur 30 × 30 km, altitude en cône
(sommet ≈ 2500 m au centre, plancher à 400 m), surface des patchs décroissante
avec l'altitude.

## Scénarios et analyses

- Quatre scénarios de réchauffement : **ΔT = 0, +2, +4, +6 °C**.
- **Analyse de sensibilité** à la capacité de dispersion : `α ∈ {0.1, 0.2, 0.3, 0.5, 0.8}` à ΔT = +4 °C.
- **Balayage du seuil d'extinction** : ΔT de 0 à 8 °C par pas de 0,5 °C, occupancy à l'équilibre.

## Exécution

```r
source("metapop_climate_mountain.R")
```

Dépendances (installées automatiquement si absentes) : `ggplot2`, `dplyr`,
`tidyr`, `patchwork`.

### Sorties produites

- `metapop_climate_composite.png` — figure composite (dynamique temporelle, seuil
  d'extinction, occupancy vs altitude, aptitude climatique vs altitude)
- `metapop_carte_scenarios.png` — cartes spatiales de l'occupancy finale par scénario
- `metapop_sensibilite_dispersion.png` — sensibilité à la dispersion
- Un résumé console : occupancy à l'équilibre et nombre de patchs quasi-éteints
  (`p < 0.05`) par scénario

## Document d'accompagnement

Le document HTML/PDF suit la structure suivante : contexte des métapopulations,
équation SPOM, sens écologique des paramètres, entrée du climat dans le modèle,
spécificités du milieu montagneux (piège sommital, connectivité altitudinale,
microrefuges), résultats de simulation, sensibilité à la dispersion, implications
pour la conservation, et le code R commenté.
