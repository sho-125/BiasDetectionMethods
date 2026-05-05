# pbiasr

`pbiasr` is a lightweight source package that wraps the 21 publication-bias methods described in the uploaded paper and code.

## Main design

- **Required input:** `yi` and `sei`
- **Optional input:** `obs`
- **Main function:** `pbias_table()`
- **Sample data:** `pbias_sample_data()`
- **Default sorting:** `Balanced Accuracy` ranking for the detected paper case

The package returns publication-bias test results only:

1. `bias_summary`: raw one-row-per-method publication-bias results
2. `publication_table`: formatted three-column display table sorted by Balanced Accuracy
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

From an unpacked source folder:

```r
remotes::install_local("pbiasr")
```

Or from a local tarball:

```r
install.packages("pbiasr_0.2.7.tar.gz", repos = NULL, type = "source")
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
- `not recommended` in the Balanced Accuracy column means the estimator did not meet the applicability/convergence criteria for that paper case.

The printed note also points users to the interactive Shiny app:

```text
https://biasdetectionmethods.shinyapps.io/main/
```

## Returned object

`pbias_table()` returns an object of class `pbias_result` with:

- `bias_summary`
- `publication_table` (Estimator / P-value / Balanced Accuracy)
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
