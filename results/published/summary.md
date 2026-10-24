# Results summary

Every value below is read from the CSV results.

| Scenario | Procedure | Quantity | Result |
| --- | --- | --- | --- |
| normal; n₁=10, n₂=20; variances 4:1 | Student | Type I error | 11.4% |
| normal; n₁=10, n₂=20; variances 4:1 | Welch | Type I error | 5.1% |
| normal; n₁=10, n₂=20; variances 4:1 | Permutation_studentized | Type I error | 5.9% |
| normal; n₁=10, n₂=20; variances 4:1 | Bootstrap_t | Type I error | 5.4% |
| lognormal; n₁=10, n₂=20; variances 4:1 | Student | Type I error | 12.6% |
| lognormal; n₁=10, n₂=20; variances 4:1 | Welch | Type I error | 11.8% |
| lognormal; n₁=10, n₂=20; variances 4:1 | Permutation_studentized | Type I error | 12.7% |
| lognormal; n₁=10, n₂=20; variances 4:1 | Bootstrap_t | Type I error | 13.1% |
| Normal; 20 independent null endpoints; n₁=n₂=30 | Welch + None | FWER | 64.6% |
| Normal; 20 independent null endpoints; n₁=n₂=30 | Welch + Holm | FWER | 5.1% |

With normal equal-variance groups of 10 and true mean difference 0.5, the mean estimate was 0.501 before selection and 1.097 among the 1833 Welch-significant results (p < 0.05).

For unequal designs, the smaller group has the larger variance. These results do not establish a universal ranking. Mann–Whitney rejection under differing distributions is not automatically type I error for equal means.
