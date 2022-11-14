
# DBCS

An R package for doing fun sequntial/anytime-valid inference in linear
models.

<https://arxiv.org/abs/2210.08639>.

All code was written and edited by Dae Woong Ham
(<daewoongham@g.harvard.edu>)

## Installation

You can install the development version from
[GitHub](https://github.com/daewoongham97/DBCS) with:

## Example

This is a basic example for constructing design-based confidence
sequences following Example 1 in the above paper.

``` r
devtools::install_github("daewoongham97/DBCS")
```

    ## Downloading GitHub repo daewoongham97/DBCS@HEAD

    ## 
    ## * checking for file ‘/private/var/folders/t7/mgw35wp13yq5hzlx18v6qcyh0000gn/T/RtmpU1qR3e/remotes193b32650a34/daewoongham97-DBCS-9c7a5f8/DESCRIPTION’ ... OK
    ## * preparing ‘DBCS’:
    ## * checking DESCRIPTION meta-information ... OK
    ## * checking for LF line-endings in source and make files and shell scripts
    ## * checking for empty or unneeded directories
    ## * building ‘DBCS_0.0.0.9000.tar.gz’

``` r
library(DBCS)

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
```
