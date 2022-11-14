
# DBCS

An R package for constructing design-based confidence sequences.
<https://arxiv.org/abs/2210.08639>.

All code was written and edited by Dae Woong Ham
(<daewoongham@g.harvard.edu>)

## Installation

You can install the development version from
[GitHub](https://github.com/daewoongham97/DBCS) with:

``` r
devtools::install_github("daewoongham97/DBCS")
library(DBCS)
```

## Example

This is a basic example for constructing design-based confidence
sequences following Example 1 in the above paper.

``` r
signup_sim = function(n, baseline = 0.15, trt_effect = -0.10, seed = 1) {
  set.seed(seed)
  Y0 = sample(c(0,1), size = n, replace = TRUE, prob = c(1 - baseline, baseline))
  Y1 = sample(c(0,1), size = n, replace = TRUE,
              prob = c(1 - baseline - trt_effect, baseline + trt_effect))

  W = sample(c(0,1), size = n, replace = TRUE)
  y = Y0
  y[W == 1] = Y1[W == 1]
  df = data.frame(W, y)
  return(df)
}

trt_effect = -0.1
baseline = 0.15
n = 10000
seed = 6
df = signup_sim(n = n, baseline = baseline, trt_effect = trt_effect, seed = seed)

treatment = "W"
response = "y"

# Exact
exact_CS = DB_CS(df, treatment, response, nonasymp = TRUE, M = 1)
# Asymptotic
asymp_CS = DB_CS(df, treatment, response)

tail(asymp_CS$lower); tail(asymp_CS$upper)
```

    ## [1] -0.1184137 -0.1184019 -0.1183900 -0.1183782 -0.1183663 -0.1183545

    ## [1] -0.07168133 -0.07167416 -0.07166699 -0.07165982 -0.07165266 -0.07164549
