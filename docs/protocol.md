# Simulation protocol

## Question and estimands

Compare frequentist performance of two-independent-sample procedures under small sample sizes, heteroscedasticity, skewness, and heavy tails. The main estimand is the arithmetic mean difference, E[X] − E[Y]. Mann–Whitney is included to illustrate both distribution sensitivity and the consequences of choosing a procedure that addresses a different question.

## Data-generating mechanisms

The main experiment crosses:

| Factor | Values |
| --- | --- |
| Smaller group size n | 10, 20, 30, 50, 100 |
| Baseline distribution | Normal; standardized lognormal; standardized t with 3 df |
| Design | Equal sizes and variances; equal sizes with variance ratio 4:1; sizes n and 2n with variance ratio 4:1 |
| Standardized mean shift δ | 0, 0.5, 0.8 |

There are 135 scenarios, with 10,000 independent datasets in each. Samples are independent within and between groups. All seven methods receive identical data within a repetition.

Baseline variables have mean zero and variance one:

- Z, with Z standard normal.
- (exp(Z) − exp(0.5)) / sqrt((e − 1)e).
- T / sqrt(3), with T distributed as t with 3 df.

Group scales are (1, 1) or (2, 1). In the unbalanced design, the smaller group has the larger variance. The opposite pairing of size and variance was not studied. The first group's mean is shifted by δ × sqrt((σ₁² + σ₂²)/2). All δ = 0 scenarios therefore have equal population means, including heteroscedastic lognormal scenarios.

## Procedures

All tests are two-sided with p < 0.05. The Mann–Whitney implementation is asymptotic with continuity correction and no ties in the simulated continuous data.

Permutation statistics are the raw difference of means and the difference divided by its Welch standard error. For each dataset, 1,999 random partitions preserve the original group sizes. Each partition is sampled independently; a partition may recur. The p-value is (b + 1)/(B + 1), with b counting absolute permuted statistics at least as extreme as the observed statistic. An absolute tolerance of 10⁻¹² is used in comparisons.

Bootstrap observations are sampled separately with replacement inside each group, using 1,999 resamples. The percentile interval uses the 0.025 and 0.975 quantiles of resampled differences. The bootstrap-t interval uses quantiles of (d* − d)/SE*, with the separate-variance standard error, and transforms them back as [d − q₀.₉₇₅SE, d − q₀.₀₂₅SE]. Quantiles use NumPy's linear interpolation. Bootstrap rejection is defined by exclusion of zero from the 95% interval; a separate bootstrap p-value is not computed. A nonfinite bootstrap pivot fails the run rather than silently dropping a resample; no such failure occurred in the published run.

## Measures

- Rejection frequency; type I error relative to equal means when δ = 0 for mean-targeting procedures, power when δ > 0.
- Coverage and mean width of 95% Student, Welch, percentile-bootstrap, and bootstrap-t intervals for the mean difference.
- Bias of the common arithmetic mean-difference estimator. This is not a separate estimator for each test.
- Mean effect estimate after selection by p < 0.05, calculated for Student and Welch. This is conditional selection bias and must not be confused with the unconditional bias of the estimator.

Mann–Whitney rejection under differing distributions is not automatically a false positive for an equal-means null. In the heteroscedastic lognormal δ = 0 scenario, means are equal but P(X < Y) ≈ 0.63624. Its usual asymptotic implementation is also not a heteroscedasticity-robust test of the weak null P(X < Y) = 0.5. Confidence intervals for the mean difference were not constructed for Mann–Whitney or permutation procedures in this experiment; corresponding CSV fields are empty.

Monte Carlo standard errors are reported. At a rejection probability of 0.05 and 10,000 repetitions, MCSE is approximately 0.00218, or 0.22 percentage points; at probability 0.5 it is 0.005, or 0.50 percentage points. Finite-resampling randomness is included in the empirically measured repeated procedure. Paired standard errors of differences between methods were not calculated. Small numerical differences do not by themselves establish superiority.

## Multiple testing

A separate experiment crosses normal/lognormal distributions, n = 10/30, equal groups with equal variances/heteroscedastic groups n and 2n, and 0/4 nonnull endpoints. Each of 16 scenarios has 10,000 independent families of 20 independent endpoints. Nonnull endpoints have δ = 0.8. Student and Welch p-values receive no adjustment, Bonferroni, Holm, or Benjamini–Hochberg adjustment.

FWER is P(V > 0); FDR is E[V/max(R,1)], with V false rejections and R total rejections. Power is the average fraction of the four nonnull endpoints detected. Under the global null FWER equals FDR. Dependence between endpoints was not simulated.

## Reproduction

The master seed is 20261024. Main scenario indices follow n → distribution → design → δ in simulate.py. A NumPy SeedSequence([20261024, scenario_index]) supplies independent child streams for data and native resampling seeds. Data generation uses NumPy's default PCG64 generator. The C++ kernel uses xorshift64* with unbiased bounded-integer draws and a separate seed per scenario. Family simulations use SeedSequence([20261024, 999, family_scenario_index]). Threads do not share random states.

The supplementary selection analysis regenerates exactly the same main datasets; it adds no independent datasets. Published run versions and settings are retained in results/published/run_metadata.json. The public CLI protects this directory from overwriting. Compiled binaries are built locally and are not distributed.

## Scope

These are synthetic continuous outcomes. No patient data were used. Paired data, clustering, missingness, discrete outcomes, alternative distributions, BCa intervals, robust estimators, and distribution-specific models were not studied. The t₃ distribution has finite variance but no finite fourth moment. Results describe these simulation conditions and do not rank methods universally. This is an exploratory reproducible analysis, not a preregistered or peer-reviewed validation.

## References

1. Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Statistics in Medicine. 2019;38:2074–2102. https://doi.org/10.1002/sim.8086
2. Fay MP, Proschan MA. Wilcoxon–Mann–Whitney or t-test? On assumptions for hypothesis tests and multiple interpretations of decision rules. Statistics Surveys. 2010;4:1–39. https://doi.org/10.1214/09-SS051
3. Chung E, Romano JP. Exact and asymptotically robust permutation tests. Annals of Statistics. 2013;41:484–507. https://doi.org/10.1214/13-AOS1090
4. Phipson B, Smyth GK. Permutation P-values should never be zero. Statistical Applications in Genetics and Molecular Biology. 2010;9:Article39. https://doi.org/10.2202/1544-6115.1585
5. Efron B, Tibshirani R. Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. Statistical Science. 1986;1:54–75. https://doi.org/10.1214/ss/1177013815
6. Holm S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics. 1979;6:65–70.
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. JRSS B. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
