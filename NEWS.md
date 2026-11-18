# pbiasr 0.4.0

- Corrected closest-condition selection: all eligible individual conditions remain separate (864 Extended or 432 Original with All filters). K is matched first, then within-condition mean estimated tau.
- Removed pooling before matching. Explicit legacy pooling requests now produce a migration error.
- Method performance is summarised after selection. Mean TPR uses selected biased conditions and mean FPR selected unbiased conditions; descriptive BA gives those classes equal weight. Both classes and complete test rates are required.
- Individual-condition BA remains unavailable. Console, app, hover and exports distinguish individual rates from the cross-condition BA comparison.

# pbiasr 0.3.2

- Added an independent `censoring` filter with `without`, `with` and `all`, alongside the same three QRP options.
- The standard workflow always groups the selected Censoring levels. Shiny no longer asks users to choose a condition-summary mode.
- QRP-only conditions continue to count as biased. Censoring filters select actual Censoring levels, independently of that classification.
- Retained legacy `universe` and `pooling` arguments for existing R scripts. Console output and app exports use QRP and Censoring terminology.

# pbiasr 0.3.1

- Closest lookup now defaults to Extended / all QRP / All conditions / pooled censoring levels.
- The biased-only universe includes QRP-only conditions (censoring or QRP defines bias), consistently for matching and performance.
- Matched mean estimated tau-squared is shown alongside mean tau. These are distinct averages.
- Added a sorted method summary: equal-weight mean Balanced Accuracy across matched groups, per-metric coverage, and applicability. Complete averages require a valid value in every match; detailed match results are retained.

# pbiasr 0.3.0

- Added `pbias_closest(k, tau)` to find and present neighbouring simulation conditions and performance of all 21 methods.
- Individual conditions are the default: mean estimated tau is computed across replications within each condition. Neighbouring matches stay separate. Optional pooling across publication-bias severity is explicit.
- Bundled summaries of the 864-condition Extended design and 432-condition Original subset, with data provenance, applicability and missing control metrics retained.
- Retained the 0.2.9 AK likelihood-ratio, coefficient p-value, TES, caliper and applicability corrections.

# pbiasr 0.2.9

- Corrected AK1/AK2 likelihood-ratio signs and degrees of freedom, coefficient tail probabilities, Hessian standard errors and invalid-fit handling.
- Matched TES and caliper defaults to the supplied simulation implementation.
- Reported unavailable regression slopes and invalid p-values explicitly.
