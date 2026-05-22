# Statistical methods
### Biosurfactant Physiology Analysis — detailed methods for manuscript and thesis

---

## 1. Data preparation (Module 1)

Raw surface tension (ST, mN/m) and optical density (OD) readings were imported from 10 independent experiments (exp2–exp11). Six readings per experimental unit were collected on Days 1, 3, 5, 7, 9, and 11. Three derived metrics were computed per strain × temperature (× salinity) combination:

**Δγ (surface tension reduction, mN/m)**
Δγ at time t = ST(Day 1) − ST(Day t)

Baseline ST was defined as the Day 1 reading per strain × temperature (× salinity) combination. Baseline homogeneity across experiments was verified prior to pooling (expected range: 68.5–71.0 mN/m).

**Δγ/OD (biomass-normalised biosurfactant efficiency)**
Δγ/OD = Δγ / ODAverage, computed for all time points where ODAverage > 0.01. Values at ODAverage ≤ 0.01 were set to NA to avoid division by near-zero denominators at the inoculation stage.

ODAverage was the mean of three technical OD replicates (OD1, OD2, OD3) measured from the same culture flask at each time point.

**AUC(Δγ/OD) — area under the curve**
AUC was computed using the trapezoidal rule over the 11-day measurement window:

AUC = Σ [(Δγ/OD_i + Δγ/OD_{i+1}) / 2] × Δt

where Δt = 2 days between consecutive readings. AUC captures sustained production across the full time course, complementing the single-point Max(Δγ/OD) metric.

---

## 2. Producer classification (Module 2)

Strains were classified into production tiers using a Gaussian mixture model (GMM) fitted to the distribution of Max(Δγ/OD) pooled across all strains and conditions (R package `mclust`, BIC-selected number of components). Classification thresholds were set at the intersection points of fitted Gaussian components, yielding data-driven rather than arbitrary quantile cutoffs.

Three tiers were defined: High Producer, Moderate Producer, Low/Non-Producer. The GMM approach was preferred over percentile-based thresholds to avoid circular dependency between the classification variable and the response variable used in downstream mixed models.

---

## 3. Temperature effects — Module 3

**Primary model (H1: cold stimulation main effect)**

```
Max(Δγ/OD) ~ Temperature + (1 | Strain) + (1 | Experiment)
```

Fitted using restricted maximum likelihood (REML) via `lme4::lmer()`. Temperature was treated as a three-level fixed effect (4°C, 15°C, 30°C). Strain was included as a random intercept to account for the inherent between-strain variance in baseline production capacity and the repeated-measures structure (each strain appears at all three temperatures). Experiment was included as a random block to absorb batch-to-batch variation arising from sequential experimental runs.

The model was run independently for three response variables: Max(Δγ/OD), AUC(Δγ/OD), and final ST (Day 11). Convergence of conclusions across all three metrics was required for a finding to be reported as robust.

**Assumption checks**
- Normality of residuals: Shapiro-Wilk test (W > 0.95) + Q-Q plot
- Homoscedasticity: Levene test per temperature group + residuals vs fitted plot
- If Max(Δγ/OD) residuals were right-skewed: log10(Max(Δγ/OD) + 1) transformation applied

**Post-hoc comparisons**
Pairwise temperature comparisons (4°C vs 15°C, 4°C vs 30°C, 15°C vs 30°C) were performed using `emmeans::emmeans()` with Tukey adjustment. Effect sizes were reported as Cohen's d computed from the estimated marginal means and pooled standard errors.

**Variance explained**
Marginal R² (fixed effects only) and conditional R² (full model) were computed using `MuMIn::r.squaredGLMM()`.

**Non-parametric backup**
If normality assumptions were violated after transformation: Kruskal-Wallis test followed by Dunn test with Bonferroni correction (`coin` package).

---

## 4. Temperature × salinity interaction — Module 4

**Primary model (H2: salinity modulation)**

```
Max(Δγ/OD) ~ Temperature * Salinity + (1 | Strain) + (1 | Experiment)
```

The interaction term Temperature:Salinity tests whether the effect of temperature on biosurfactant production is modulated by salinity concentration. A significant positive interaction indicates synergistic combined stress; a significant negative interaction indicates antagonism.

For the 11 strains with data in both temperature-only and temperature × salinity experiments, the temperature-only data served as a "zero salinity" reference, enabling the full 4-level salinity comparison (0M, 0.5M, 1.0M, 1.5M) for those strains.

**Bliss Independence model for synergy testing**
Synergy between cold temperature and high salinity was formally tested using the Bliss Independence model. Each stress factor's effect was normalised to a fractional inhibition/activation scale (0–1). The expected combined effect under independence is:

E(A+B) = E(A) + E(B) − E(A) × E(B)

Synergy is defined as observed effect > E(A+B); antagonism as observed effect < E(A+B). This approach is preferred over custom synergy indices because it has an established statistical framework and is directly comparable to published combinatorial stress studies.

**Simple effects decomposition**
Significant interactions were decomposed using simple effects analysis: the effect of temperature was estimated separately at each salinity level, and vice versa.

---

## 5. Strain classification — Module 5

**Cold Stimulation Index (CSI)**
CSI = AUC(Δγ/OD) at 4°C / AUC(Δγ/OD) at 30°C

Strains were classified as Cold-stimulated (CSI > 1.5), Neutral (0.7 ≤ CSI ≤ 1.5), or Warm-preferred (CSI < 0.7).

**Hierarchical clustering**
Strains were clustered on a matrix of Max(Δγ/OD) values across all temperature × salinity combinations using Ward's D2 linkage on Euclidean distances (`pheatmap::pheatmap()`). Cluster stability was assessed by bootstrapping (1000 iterations, `pvclust` package).

**OD vs production decoupling**
Spearman correlation between Max(OD) and Max(Δγ/OD) per temperature was computed to test whether biosurfactant production is independent of growth rate. A non-significant or weak correlation formally rules out the confound that high production merely reflects high biomass.

---

## 6. Genomic integration — Module 6

**Lipopeptide BGC vs production**
Strains were coded as BGC-positive (Lichenysin or Surfactin cluster present and intact, from antiSMASH 7.1.0) or BGC-negative. A Wilcoxon rank-sum test compared Max(Δγ/OD) at 4°C between the two groups. Point-biserial correlation (equivalent to rank-biserial correlation for Wilcoxon) provided the effect size.

**Spearman correlation: BGC count vs production**
The total number of BGCs per strain (from WGS pipeline) was correlated with Max(Δγ/OD) at 4°C using Spearman's rho.

**Stress gene linkage**
Biosurfactant regulation gene counts (from stress gene analysis in WGS pipeline) were correlated with Max(Δγ/OD) at 4°C and AUC(Δγ/OD) at 4°C using Spearman's rho.

**Within-species comparisons**
For *B. zhangzhouensis* strains (n = 7–8 depending on salinity data availability), phylogenetic position (from IQ-TREE2 core genome tree) was plotted alongside CSI values to test whether cold-inducibility is a derived trait within the clade.

---

## Software versions

| Software | Version | Purpose |
|---|---|---|
| R | ≥ 4.3.0 | All statistical analysis |
| lme4 | ≥ 1.1-35 | Linear mixed-effects models |
| lmerTest | ≥ 3.1-3 | F-tests and p-values for lme4 |
| emmeans | ≥ 1.10 | Post-hoc pairwise comparisons |
| MuMIn | ≥ 1.47 | R² for mixed models |
| mclust | ≥ 6.0 | Gaussian mixture model classification |
| coin | ≥ 1.4 | Non-parametric tests |
| ggplot2 | ≥ 3.5 | Visualisation |
| patchwork | ≥ 1.2 | Multi-panel figures |
| pheatmap | ≥ 1.0.12 | Heatmaps |
| ComplexHeatmap | ≥ 2.18 | Advanced heatmaps |

