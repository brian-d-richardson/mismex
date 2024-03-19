
# mismex: causal inference with a mismeasured exposure <img id="mismex_hex" src="man/figures/mismex_hex.png" align="right" width="125"/>

Brian D. Richardson

## Installation

Installation of `mismex` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done with the following code:

``` r
# install the package
devtools::install_github(repo = "brian-d-richardson/mismex", 
                         ref = "main")
```

``` r
# load the package
library(mismex)

# load additional packages
library(MASS)
library(dplyr)
```

The `mismex` package contains functions to estimate causal effects in
the presence of confounding and a mismeasured exposure. The methods
implemented are introduced in the paper, “Addressing confounding and
continuous exposure measurement error using corrected score functions,”
which is currently in progress.

## Example

An example of the 3 proposed estimators (g-formula, IPW, and doubly
robust) used on a simulated data set is provided below.

### Data Generation

Data are generated below according to the data generation process
described in the third simulaton study.

``` r
## define parameters
 
n = 2000                                  # sample size
B = 30                                    # number of MC replicates
seed = 1                                  # random number seed
mc.seed <- 123                            # MCCS seed
cov.e <- 0.16                             # var(epsilon)
inv.link <- inv.ident                     # inverse link
d.inv.link <- d.inv.ident                 # derivative of inverse link
g <- c(1.5, 0.7, 0.9, -0.6, -0.7, 0.4)    # outcome model parameters
formula <- "~A*L1 + A*L2"                 # outcome model formula
ps.formula <- "~L1 + L2"                  # propensity score model formula
ipw.formula <- "~A"                       # MSM formula
```

``` r
## generate data

set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
L2 <- rnorm(n, 1, sqrt(0.5))                             # confounder 2
A <- rnorm(n, 2 + 0.9*L1 - 0.6*L2, sqrt(1.1))            # exposure
a <- seq(min(A), max(A), length = 4)                     # grid of exposures
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
Y <- rnorm(n, EY, sqrt(0.16))                            # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                    # mismeasured A
dat0 <- data.frame(Y, A, L1, L2)                         # oracle data
datstar <- data.frame(Y, A = Astar, L1, L2)              # mismeasured data
args <- list(formula = formula,                          # arguments for fitting
             ps.formula = ps.formula,
             inv.link = inv.link,
             d.inv.link = d.inv.link)
```

### G-Formula Estimation

The dose response curve at chosen values of the exposure $a$ are
estimated here using the MCCS g-formula method.

``` r
## G-formula
gfmla.res <- fit.gfmla.mccs(
  data = datstar,
  a = a,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = list(formula = formula,   
              inv.link = inv.link,
              d.inv.link = d.inv.link))

cbind(est = round(gfmla.res$est, 2),
      stde = round(sqrt(diag(gfmla.res$var)), 2))
```

    ##         est stde
    ## g.0    1.52 0.05
    ## g.1    0.69 0.02
    ## g.2    0.88 0.05
    ## g.3   -0.58 0.03
    ## g.4   -0.70 0.02
    ## g.5    0.39 0.01
    ## EYa.1 -0.48 0.06
    ## EYa.2  1.57 0.02
    ## EYa.3  3.62 0.03
    ## EYa.4  5.67 0.07

### IPW Estimation

Parameters in the marginal structural model are estimated here using the
MCCS IPW method.

``` r
## IPW
ipw.res <- fit.ipw.mccs(
  data = datstar,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = list(formula = ipw.formula,   
              ps.formula = ps.formula,
              inv.link = inv.link,
              d.inv.link = d.inv.link))

cbind(est = round(ipw.res$est, 2),
      stde = round(sqrt(diag(ipw.res$var)), 2))
```

    ##                est stde
    ## g.0           1.78 0.19
    ## g.1           0.55 0.09
    ## coef.a.l.1    2.06 0.05
    ## coef.a.l.2    0.88 0.05
    ## coef.a.l.3   -0.66 0.04
    ## log.var.a.l1  0.13 0.04

### Double Robust Estimation

The dose response curve is again estimated here using the double robust
method.

``` r
## Double Robust
dr.res <- fit.dr.mccs(
  data = datstar,
  a = a,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = list(formula = formula,   
              ps.formula = ps.formula,
              inv.link = inv.link,
              d.inv.link = d.inv.link))

cbind(est = round(dr.res$est, 2),
      stde = round(sqrt(diag(dr.res$var)), 2))
```

    ##                est stde
    ## g.0           1.51 0.05
    ## g.1           0.70 0.02
    ## g.2           0.89 0.05
    ## g.3          -0.59 0.03
    ## g.4          -0.72 0.02
    ## g.5           0.40 0.02
    ## coef.a.l.1    2.06 0.05
    ## coef.a.l.2    0.88 0.05
    ## coef.a.l.3   -0.66 0.04
    ## log.var.a.l1  0.13 0.04
    ## EYa.1        -0.52 0.07
    ## EYa.2         1.56 0.02
    ## EYa.3         3.64 0.03
    ## EYa.4         5.72 0.07
