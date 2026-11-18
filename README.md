# pbiasr

`pbiasr` is a lightweight source package that wraps the 21 publication-bias methods described in the uploaded paper and code.

## Main design

- **Required input:** `yi` and `sei`
- **Optional input:** `obs`
- **Main function:** `pbias_table()`
- **Sample data:** `pbias_sample_data()`
- **Default sorting:** `Logit Distance` ranking for the detected paper case

The package returns publication-bias test results only:

1. `bias_summary`: raw one-row-per-method publication-bias results
2. `publication_table`: formatted four-column display table sorted by the selected metric
3. `case`: paper case based on number of studies and I-squared
4. `interpretation`: short interpretation text for the publication-bias tests
5. `meta`: run metadata

Estimated effects and model-parameter estimates are not reported in the returned object or console output.

## Paper case and ranking

`pbias_table()` fits a baseline random-effects model internally to estimate I-squared for case assignment only. It does not report the estimated effect.

The paper case is defined by:

| Code | Number of studies |
|---|---|
| S | `K <= 10` |
| M | `10 < K <= 100` |
| L | `K > 100` |

and:

| Code | I-squared |
|---|---|
| S | `I2 < 0.29` |
| M | `0.29 <= I2 < 0.60` |
| L | `I2 >= 0.60` |

For example, `ML` means medium-sized meta-analysis with large heterogeneity. The performance values are bundled with the package and come from **OUR PUBLICATION BIAS** paper.

## Important note about input requirements

Your uploaded code implements EGGER4 / EGGER5 via `MAIVE` and `FATIV`. Those methods use `obs` (a sample-size proxy). Because of that, **it is not possible to run all 21 methods from effect size and standard error alone**. In this package:

- methods that need only `yi` and `sei` run directly
- EGGER4 / EGGER5 run only when `obs` is supplied
- otherwise they return `NA` with a clear status message

## Installation

Install the current release from GitHub:

```r
remotes::install_github("sho-125/BiasDetectionMethods")
```

From an unpacked source folder:

```r
remotes::install_local("pbiasr")
```

Or from a local tarball:

```r
install.packages("pbiasr_0.4.1.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(pbiasr)

sample_data1 <- pbias_sample_data(1)
str(sample_data1)

out1 <- pbias_table(
  yi = sample_data1$eff,
  sei = sample_data1$se,
  obs = sample_data1$obs
)

print(out1)
```

Second sample dataset:

```r
sample_data2 <- pbias_sample_data(2)
str(sample_data2)

out2 <- pbias_table(
  yi = sample_data2$eff,
  sei = sample_data2$se,
  obs = sample_data2$obs
)

print(out2)
```

Or load both sample datasets at once:

```r
samples <- pbias_sample_data("all")
names(samples)
```


The printed publication-bias table follows this interpretation:

- **Null hypothesis:** no publication bias
- **Rejection:** the method detects the existence of publication bias
- **Estimated p-value:** calculated from the user's meta-analysis data
- **Balanced Accuracy:** based on the corresponding results from **OUR PUBLICATION BIAS** paper
- If a method cannot be run, `P-value` shows a short status such as `requires obs`, `n/a`, or `failed`.
- `Unranked*` in the Balanced Accuracy column means the estimator did not meet the applicability/convergence criteria for that paper case.

The printed note also points users to the interactive Shiny app:

```text
https://pbiasmethods.shinyapps.io/main/
```

## Returned object

`pbias_table()` returns an object of class `pbias_result` with:

- `bias_summary`
- `publication_table` (Method / P-value / Balanced Accuracy / Logit Distance)
- `case`
- `interpretation`
- `meta`

`pbias_sample_data(1)` loads `SampleData1.RData`, `pbias_sample_data(2)` loads `SampleData2.RData`, and `pbias_sample_data("all")` loads both bundled sample datasets from `inst/extdata`.

## Method mapping

| Paper label | Code label in package |
|---|---|
| BEGG | `Begg` |
| EGGER1 | `FE` |
| EGGER2 | `RE` |
| EGGER3 | `WLS` |
| EGGER4 | `MAIVE Type I` |
| EGGER5 | `MAIVE Type II` |
| EK | `EK` |
| SKEW1 | `SKEWNESS` |
| SKEW2 | `SKEWCombined` |
| PSM3 | `3PSM` |
| PSM4 | `4PSM` |
| AK1 | `AK1` |
| AK2 | `AK2` |
| PUNIF | `Puniform` |
| TES | `TES` |
| PSST | `PSST` |
| TESS | `TESS` |
| CALI05 | `Caliper05` |
| CALI10 | `Caliper10` |
| CALI15 | `Caliper15` |
| CALI20 | `Caliper20` |

## Numerical corrections in 0.2.9

- AK1 and AK2 convert the minimized negative log-likelihood back to a log-likelihood before computing `2 * (logLik_selected - logLik_no_selection)`. AK1 uses 1 degree of freedom and AK2 uses 3, matching the supplied simulation implementation.
- AK coefficient p-values use Student t tail probabilities. Hessian standard errors use the square roots of the covariance diagonal, avoiding square roots of negative off-diagonal entries. Invalid fits and missing selection intervals are reported explicitly.
- Constant standard errors do not yield valid Egger slope tests. Constant sample sizes are marked unavailable for MAIVE/FATIV. Invalid p-values are reported as unavailable rather than successful results.
- Internal selection-model coefficient p-values use two-sided normal tails. This does not change their publication-bias likelihood-ratio tests.

The default sort is Logit Distance. Use `ranking_metric = "balanced_accuracy"` to sort by Balanced Accuracy instead. The reference I-squared decision rule is based on the Size x I2 (Extended) simulation results.

TES uses the reference right-tail alpha of 0.025 when `alpha = 0.05`. Default caliper windows reproduce the literal 1.96-centered boundaries in the simulation code (for example CALI05 uses 1.86-1.96 and 1.96-2.06), including the five-observation applicability rule.


## Closest individual simulation conditions (0.4.1)

```r
near <- pbias_closest(k = 65, tau = 0.2)
near                       # a readable summary of the matched conditions
near$summary               # mean Balanced Accuracy across matches, sorted
near$performance           # all 21 methods for each match, as numeric data
near$matches               # the full experimental settings

# If your heterogeneity estimate is a variance, convert it to tau:
pbias_closest(k = 65, tau = sqrt(0.04), censoring = "all")

# Restrict to the original design and medium/high p-hacking:
pbias_closest(k = 30, tau = 0.2, design = "original", qrp = "with")
```

The defaults are **Extended / QRP All / Censoring All**. Matching searches all **864 individual conditions**, or 432 under Original. Each full K / true-effect / population-tau / QRP / Censoring condition retains its own mean estimated tau. There is no pooling across Censoring or QRP levels before matching.

K is matched first: exact K, or nearest lower/upper K. Within each selected K, the immediately lower and upper mean estimated tau are selected. Exact matches take priority; ties are retained. For K=65 and tau=0.2, this searches 108 conditions at K=60 and 108 at K=100, then selects four individual conditions.

QRP and Censoring filters independently accept `without` (none), `with` (medium/high) and `all`. Mean tau and mean tau-squared average valid simulation replications within each condition. The latter is the mean estimated variance, not squared mean tau. `matches$n_tau` records the valid count.

```r
near <- pbias_closest(65, 0.2, qrp = "all", censoring = "all")
near$matches     # four separate neighbours, with explicit QRP/Censoring labels and full settings
near$performance # 21 methods per condition: TPR or FPR and applicability
near$summary     # descriptive performance AFTER selecting neighbours
```

A condition is biased when QRP or Censoring is medium/high. A biased condition supplies TPR; an unbiased condition supplies FPR. Individual-condition BA is unavailable. After selecting neighbours, summary mean TPR averages the biased conditions' rates and mean FPR averages the unbiased conditions' rates. **BA = (mean TPR + 1 - mean FPR) / 2**. Each class gets half the weight, and conditions within a class receive equal weight. Both classes and a valid rate in every selected condition are required; missing results are not filled with zero. This describes the selected simulation neighbours and is not an interpolated prediction for the user's meta-analysis. Applicability and missing-result coverage are retained.

If only one class is selected, BA remains `NA` and mean rejection rates remain available. Hover and CSV exports identify every condition's QRP/Censoring levels and its TPR/FPR contribution. Low applicability is flagged separately. The Shiny app uses this same package function.

`pooling = "none"` is the only supported selection mode; older requests for `pooling = "publication_bias"` fail clearly rather than silently changing the search. The legacy `universe` filter is still supported. Raw test corrections from 0.2.9 remain in place.
