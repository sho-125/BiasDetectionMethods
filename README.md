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
install.packages("pbiasr_0.3.1.tar.gz", repos = NULL, type = "source")
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


## Closest simulation conditions (0.3.0)

```r
near <- pbias_closest(k = 65, tau = 0.2)
near                       # a readable summary of the matched conditions
near$summary               # mean Balanced Accuracy across matches, sorted
near$performance           # all 21 methods for each match, as numeric data
near$matches               # the full experimental settings

# If your heterogeneity estimate is a variance, convert it to tau:
pbias_closest(k = 65, tau = sqrt(0.04), universe = "all")

# Restrict to the original design and medium/high p-hacking:
pbias_closest(k = 30, tau = 0.2, design = "original", qrp = "with")
```

The defaults are **Extended / all QRP levels / All conditions / pooled censoring levels**. K is matched first: an exact K, or the nearest simulated K below and above. Within each selected K, the function returns the immediately lower and upper mean estimated tau. Exact matches take priority; ties are retained, and boundary messages explain when only one side is available.

**Grouped results are the default.** Each group holds K, true effect, population tau and QRP fixed. `universe = "all"` averages all three censoring levels. `universe = "pb_only"` includes every biased condition: censoring medium/high OR QRP medium/high. Thus PB-only groups contain two censoring levels at QRP none and all three at positive QRP. Both universes have the same K / effect / population-tau / QRP groups.

Use `pooling = "none"` for individual conditions, `qrp = "without"` for QRP none, or `qrp = "with"` for medium/high QRP. QRP levels remain separate in all modes.

Within a full condition, mean tau averages valid estimates among the 1,000 replications. `matches$n_tau` records the valid count (667 of 864,000 source tau estimates are missing). Group means give equal weight to each eligible censoring condition. `mean_tau` is the mean estimated standard deviation; `mean_tau2` is the mean estimated variance, and differs from `mean_tau^2`.

`near$summary` gives one row per method, sorted by mean Balanced Accuracy. Each matched group receives equal weight. An average is available only when every selected match has a finite value; `ba_available` and `n_matches` show coverage. Per-match values remain in `near$performance`. Low applicability is flagged separately, including the number of matches below 90%.

Rejection rates exclude missing tests. Applicability reports the fraction of all replications with valid results. TPR requires biased conditions; FPR requires unbiased conditions; Balanced Accuracy requires both. Missing metrics remain `NA`, including for biased-only groups and individual conditions. No additional control conditions are invented. All numeric rates are proportions. Pooled rates give equal weight to available component-condition rates.

The summaries include all 864 Extended conditions. Original restricts K to 10, 30, 60 and 100 (432 conditions based on Carter et al., 2019). The Shiny app uses this same package function for matching and its average-performance chart; hovering a bar shows each match's value.
