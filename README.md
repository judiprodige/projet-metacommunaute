# projet-metacommunaute

Projet test dans le cadre d'une école d'été sur l'écologie spatiale.

Ce dépôt contient un **modèle de métapopulation spatialement explicite** (SPOM,
*Stochastic Patch Occupancy Model* de Hanski) appliqué à un paysage **montagneux**
soumis au **changement climatique**, son **analyse analytique** (équilibre, capacité
de métapopulation, contribution des patchs), et un document pédagogique.

## Contenu du dépôt

| Fichier | Description |
|---|---|
| [metapop_climate_mountain.R](metapop_climate_mountain.R) | Paysage, simulation par Euler, scénarios climatiques, sensibilité, figures |
| [metapop_capacity_eigen.R](metapop_capacity_eigen.R) | Équilibre analytique, capacité λ_M, analyse spectrale et contribution des patchs |
| [metapop_velocite_climatique.R](metapop_velocite_climatique.R) | Réchauffement transitoire, temps de relaxation, vélocités climat/colonisation, dette d'extinction |
| [Modèle de métapopulation spatialement explicite.html](Modèle%20de%20métapopulation%20spatialement%20explicite.html) | Document explicatif en 12 sections (théorie, résultats, équilibre, analyse spectrale, vélocités, code) |
| [Modèle de métapopulation spatialement explicite.pdf](Modèle%20de%20métapopulation%20spatialement%20explicite.pdf) | Version PDF — **antérieure**, ne contient pas les sections 8 à 10 |
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
2. **Niche thermique gaussienne** : aptitude `s = exp(−0.5·((T − T_opt)/T_sd)²)`,
   avec `T_opt = 6 °C` et `T_sd = 3 °C`.
3. **Modulation des paramètres** : `A_eff = A_base · s`, `c_x = c0 · s`,
   `e_x = e0 / A_eff`. Un réchauffement `ΔT` décale la température locale et pousse
   la zone favorable vers les sommets (*piège sommital* : moins de surface, moins
   de connectivité en altitude).

### Paysage simulé

30 patchs tirés aléatoirement (`set.seed(42)`) sur 30 × 30 km, altitude en cône
(sur ce tirage : 941 – 2246 m), surface des patchs décroissante avec l'altitude,
noyau de dispersion `α = 0,2`.

## Résultats analytiques

Le modèle admet une solution d'équilibre explicite, vérifiée contre la simulation
d'Euler (écart < 10⁻¹³ par patch) :

```
p*_x = S_x / (S_x + δ_x)      S_x = Σ_y K(x,y)·A_eff,y·p_y      δ_x = e0 / (c0·A_x·s_x²)
```

- Le climat entre **au carré** dans `δ_x` : il ralentit la colonisation *et* accélère
  l'extinction via la surface effective.
- L'équilibre positif est **unique et globalement attractif** (l'itération est concave),
  donc `p_init` n'influence pas le résultat. Le passage à l'extinction est une
  bifurcation **transcritique** — pas de bistabilité ni d'hystérésis.

Persistance ⟺ `λ_M > e0/c0`, où `λ_M` est la plus grande valeur propre de
`M_xy = (A_x·s_x^1.5) · K(x,y) · (A_y·s_y^1.5)`. Son vecteur propre `w` donne la
contribution `w²_x` de chaque patch (somme = 1), et toutes les élasticités :
`∂lnλ_M/∂lnA_x = 2w²_x`, `∂lnλ_M/∂ln s_x = 3w²_x`.

| ΔT | p̄* | patchs `p<0,05` | λ_M | marge / seuil | top 5 patchs | altitude pondérée w² |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0,981 | 0 / 30 | 76,3 | 191 × | 71 % | 1391 m |
| +2 | 0,965 | 0 / 30 | 43,9 | 110 × | 61 % | 1681 m |
| +4 | 0,774 | 0 / 30 | 21,6 | 54 × | 59 % | 1949 m |
| +6 | 0,416 | 8 / 30 | 6,5 | 16 × | 74 % | 2092 m |
| **+8,42** | **0** | 30 / 30 | **0,40** | **1 ×** | — | — |

Seuil d'extinction du paysage : **ΔT\* = 8,42 °C**. À +6 °C, un seul patch de crête
(#15, 2176 m) porte 22 % de la capacité du réseau, et dix patchs en portent 99 %.
Le cœur fonctionnel remonte de 1391 m à 2092 m, soit ≈ 117 m par degré.

## Analyses incluses

- Quatre scénarios de réchauffement : **ΔT = 0, +2, +4, +6 °C**.
- **Sensibilité à la dispersion** : `α ∈ {0.1, 0.2, 0.3, 0.5, 0.8}` à ΔT = +4 °C
  (p̄ final de 0,88 à 0,26 — la connectivité pèse autant que le climat).
- **Balayage du seuil** : ΔT de 0 à 8 °C par pas de 0,5 °C.
- **Analyse spectrale** : λ_M, écart spectral λ₂/λ₁ (0,48 → 0,31 : le réseau se réduit
  à un cœur unique), contribution par patch, validation par retrait effectif du patch pivot
  (−18 % sur λ_M, là où l'élasticité `2w²` prédit −32 %).
- **Vélocités** : réchauffement transitoire `ΔT(t) = r·t`, temps de relaxation τ, comparaison
  `v_clim = 167·r` m/ut contre `v_col = 1/(0,006·τ)`, et retard `≈ r·τ`.

## Exécution

```r
source("metapop_climate_mountain.R")     # simulation + figures
source("metapop_capacity_eigen.R")       # équilibre, capacité, contributions
source("metapop_velocite_climatique.R")  # vélocités, transitoire, dette d'extinction
```

Le second script source le premier ; les deux s'exécutent depuis la racine du dépôt.
Dépendances (installées automatiquement si absentes) : `ggplot2`, `dplyr`, `tidyr`,
`patchwork`.

### Sorties produites

`metapop_climate_mountain.R` :

- `metapop_climate_composite.png` — dynamique temporelle, seuil d'extinction,
  occupancy vs altitude, aptitude climatique vs altitude
- `metapop_carte_scenarios.png` — cartes spatiales de l'occupancy finale par scénario
- `metapop_sensibilite_dispersion.png` — sensibilité à la dispersion
- Console : occupancy à l'équilibre et patchs quasi-éteints (`p < 0.05`) par scénario

`metapop_capacity_eigen.R` :

- `metapop_contribution_patchs.png` — cartes de contribution `w²` par scénario,
  contribution selon l'altitude, courbe de contribution cumulée
- Console : λ_M, λ₂/λ₁, concentration, patchs pivots et élasticités climatiques

`metapop_velocite_climatique.R` :

- `metapop_velocite_climatique.png` — croisement des vélocités, ralentissement critique,
  dette d'extinction, écart optimum / occupation
- Console : τ et v_col par scénario, retard accumulé selon le rythme de réchauffement
- Fonction réutilisable `simulate_transient(r)` pour un réchauffement progressif

## Document d'accompagnement

Le document HTML est structuré en 11 sections : contexte des métapopulations,
équation SPOM, sens écologique des paramètres, entrée du climat dans le modèle,
spécificités du milieu montagneux (piège sommital, connectivité altitudinale,
microrefuges), résultats de simulation, sensibilité à la dispersion,
**conditions d'équilibre**, **capacité de métapopulation et analyse spectrale**,
**vélocité climatique vs. vélocité de colonisation**, implications pour la conservation,
et le code R commenté.

Une **présentation en 8 diapositives** reprend l'ensemble sous forme visuelle
(contexte, espèces cibles illustrées, modèle, scénarios, cartes et dispersion, patchs pivots,
vélocités, pistes de validation) — voir l'artéfact publié « Archipel Vertical ».
