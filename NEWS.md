# pbiasr 0.3.0

- Added `pbias_closest(k, tau)` to find and present neighbouring simulation conditions and performance of all 21 methods.
- Individual conditions are the default: mean estimated tau is computed across replications within each condition. Neighbouring matches stay separate. Optional pooling across publication-bias severity is explicit.
- Bundled summaries of the 864-condition Extended design and 432-condition Original subset, with data provenance, applicability and missing control metrics retained.
- Retained the 0.2.9 AK likelihood-ratio, coefficient p-value, TES, caliper and applicability corrections.

# pbiasr 0.2.9

- Corrected AK1/AK2 likelihood-ratio signs and degrees of freedom, coefficient tail probabilities, Hessian standard errors and invalid-fit handling.
- Matched TES and caliper defaults to the supplied simulation implementation.
- Reported unavailable regression slopes and invalid p-values explicitly.
