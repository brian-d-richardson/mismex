
- [mismex: causal inference with a mismeasured exposure
  <img id="mismex_hex" src="man/figures/mismex_hex.png" align="right" width="125"/>](#mismex-causal-inference-with-a-mismeasured-exposure-)
  - [Installation](#installation)
  - [Examples](#examples)
    - [Data Generation](#data-generation)
    - [G-Formula Estimation](#g-formula-estimation)
    - [IPW Estimation](#ipw-estimation)
    - [Double Robust Estimation](#double-robust-estimation)
    - [Nonlinear Marginal Structural
      Model](#nonlinear-marginal-structural-model)
    - [Case Cohort Sampling](#case-cohort-sampling)
  - [Simulation Study](#simulation-study)
  - [Application to HVTN 505 Trial
    Data](#application-to-hvtn-505-trial-data)

# mismex: causal inference with a mismeasured exposure <img id="mismex_hex" src="man/figures/mismex_hex.png" align="right" width="125"/>

Brian D. Richardson

## Installation

Installation of `mismex` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done with the following code:

``` r
## install the package
devtools::install_github(repo = "brian-d-richardson/mismex", 
                         ref = "main")
```

``` r
## load the package
library(mismex)
library(devtools)

## load additional packages
library(MASS)
library(dplyr)
library(tidyverse)
library(ggplot2)
```

The `mismex` package contains functions to estimate causal effects in
the presence of confounding and a mismeasured exposure. The methods
implemented are introduced in the paper, “Addressing confounding and
continuous exposure measurement error using corrected score functions,”
which is currently in progress.

## Examples

An example of the three proposed estimators (g-formula, IPW, and doubly
robust) used on a simulated data set is provided below.

### Data Generation

Data are generated below according to the data generation process
described in the third simulation study:

- sample size $n=2000$,
- two independent confounders
  $L_1 \sim \textrm{Bernoulli}(0.5), L_2 \sim N(0, 0.16)$
- univariate exposure $A$, where
  $A|\pmb{L} \sim N(0.1 - 0.1L_1 + 0.3L_2, 0.04)$
- binary outcome $Y$ with
  $Y|A,\pmb{L} \sim \textrm{Bernoulli}(\beta_0 + \beta_1 A + \beta_2 L_1 + \beta_3 L_2 + \beta_4 AL_1 + \beta_5 AL_2)$,
  where
  $\pmb{\beta} = (\beta_0, \beta_1, \beta_2, \beta_3, \beta_4, \beta_5) = (0.35, 0.15, 0.25, 0.2, 0.05, 0.1)$
- resulting MSM $\eta(a;\pmb{\gamma}) = \gamma_0 + \gamma_1 a$, where
  $\pmb{\gamma} = (\gamma_0, \gamma_1) = (0.475, 0.175)$
- mismeasured exposure $A^* = A + \epsilon$, where
  $\epsilon \sim N(0, 0.02)$

``` r
## define parameters
n = 2000                                  # sample size
seed = 1                                  # random number seed
mc.seed <- 123                            # MCCS seed
cov.e <- 0.02                             # var(epsilon)
inv.link <- inv.ident                     # inverse link
d.inv.link <- d.inv.ident                 # derivative of inverse link
g <- c(0.35, 0.15, 0.25, 0.2, 0.05, 0.1)  # outcome model parameters
formula <- "~A*L1 + A*L2"                 # outcome model formula
ps.formula <- "~L1 + L2"                  # propensity score model formula
ipw.formula <- "~A"                       # MSM formula
```

``` r
## generate data

set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
L2 <- rnorm(n, 0, 0.16)                                  # confounder 2
EA <- 0.1 - 0.1*L1 + 0.3*L2                              # E(A|L)
A <- rnorm(n, EA, sqrt(0.04))                            # exposure
a <- seq(-0.5, 0.5, length.out = 10)                     # grid of exposures
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
EY[EY < 0] <- 0; EY[EY > 1] <- 1                         # constrain EY
Y <- rbinom(n, 1, EY)                                    # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                    # mismeasured A
datstar <- data.frame(Y, A = Astar, L1, L2)              # mismeasured data

head(datstar, 5)
```

    ##   Y           A L1          L2
    ## 1 0  0.21478594  0  0.18159441
    ## 2 1  0.07379852  0  0.17790910
    ## 3 0 -0.08944964  1 -0.13932442
    ## 4 1 -0.15502348  1  0.03371705
    ## 5 0 -0.25402005  0  0.01110330

### G-Formula Estimation

The dose response curve at chosen values of the exposure $a$ are
estimated here using the MCCS g-formula method.

Before fitting the model, we determine an appropriate number of
Monte-Carlo replicates $B$ in order for the G-formula MCCS function to
approximate the CS function. We can do this by evaluating the MCCS
function for a sequence of $B$ values at a particular parameter value,
say the naive g-formula estimator (ignoring measurement error).

``` r
## g-formula arguments
gfmla.args <- list(formula = formula,   
                   inv.link = inv.link,
                   d.inv.link = d.inv.link)

## naive estimator
gfmla.naive <- fit.glm(data = datstar,
                       args = gfmla.args,
                       return.var = F)$est

## assess MCCS estimating function over grid of B values
gfmla.B.tuning <- tune.B(
  get.psi = get.psi.glm,
  data = datstar,
  cov.e = cov.e,
  BB = seq(1, 50, by = 2),
  args = gfmla.args,
  mc.seed = 123)

gfmla.B.tuning$plot
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Based on these plots, the G-formula MCCS with $B=30$ seems sufficient.
Using this value, we proceed with estimating the G-formula parameters.
These include the parameters $\pmb{\beta} = (\beta_0, \dots, \beta_5)$
in the outcome model $Y|\pmb{L},A$, and the dose response curve
$\textrm{E}\{Y(a)\}$ evaluated at four points (-0.5, -0.39, -0.28,
-0.17, -0.06, 0.06, 0.17, 0.28, 0.39, 0.5) in the support of $A$.

The function `fit.gfmla.mccs` returns a list with three items: the
parameter estimates, their estimated covariance matrix, and a
bias-corrected covariance matrix estimate. Below a table is shown of
estimates and (crude and bias-corrected) standard errors.

``` r
## number of MC replicates
B <- 30

## G-formula
gfmla.res <- fit.gfmla.mccs(
  data = datstar,
  a = a,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = gfmla.args)

cbind(est = round(gfmla.res$est, 2),
      stde = round(sqrt(diag(gfmla.res$var)), 2),
      bc.stde = round(sqrt(diag(gfmla.res$bc.var)), 2))
```

    ##         est stde bc.stde
    ## g.0    0.37 0.02    0.02
    ## g.1    0.02 0.09    0.09
    ## g.2    0.23 0.02    0.02
    ## g.3    0.11 0.07    0.07
    ## g.4    0.23 0.13    0.13
    ## g.5    0.75 0.35    0.35
    ## EYa.1  0.41 0.04    0.04
    ## EYa.2  0.43 0.03    0.03
    ## EYa.3  0.44 0.02    0.02
    ## EYa.4  0.46 0.02    0.02
    ## EYa.5  0.47 0.01    0.01
    ## EYa.6  0.49 0.01    0.01
    ## EYa.7  0.50 0.01    0.01
    ## EYa.8  0.51 0.02    0.02
    ## EYa.9  0.53 0.02    0.02
    ## EYa.10 0.54 0.03    0.03

### IPW Estimation

Parameters $\pmb{\gamma} = (\gamma_0, \gamma_1)$ in the marginal
structural model, as well as coefficients and variance in the propensity
score models for $A|\pmb{L}$ and $A$, are estimated here using the MCCS
IPW method. We use the same number of MC replicates $B=30$ as for IPW,
but a similar strategy as with g-formula could be used to tune this to
an appropriate number.

``` r
## IPW arguments
ipw.args <- list(formula = ipw.formula,   
                 ps.formula = ps.formula,
                 inv.link = inv.link,
                 d.inv.link = d.inv.link)

## IPW estimation
ipw.res <- fit.ipw.mccs(
  data = datstar,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = ipw.args)

cbind(est = round(ipw.res$est, 2),
      stde = round(sqrt(diag(ipw.res$var)), 2),
      bc.stde = round(sqrt(diag(ipw.res$bc.var)), 2))
```

    ##                est stde bc.stde
    ## g.0           0.48 0.01    0.01
    ## g.1           0.14 0.08    0.08
    ## coef.a.l.1    0.10 0.01    0.01
    ## coef.a.l.2   -0.10 0.01    0.01
    ## coef.a.l.3    0.27 0.03    0.03
    ## log.var.a.l1 -3.17 0.05    0.05
    ## mean.a.1      0.05 0.01    0.01
    ## cov.a.1       0.05 0.00    0.00

### Double Robust Estimation

Finally, we estimate outcome model parameters, propensity model
parameters, and the dose response curve using the doubly robust MCCS
method, again with $B=30$ replicates.

``` r
## arguments for double robust estimation
dr.args <- list(formula = formula, 
                ps.formula = ps.formula,
                inv.link = inv.link,
                d.inv.link = d.inv.link)

## Double Robust
dr.res <- fit.dr.mccs(
  data = datstar,
  a = a,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = dr.args)

cbind(est = round(dr.res$est, 2),
      stde = round(sqrt(diag(dr.res$var)), 2),
      bc.stde = round(sqrt(diag(dr.res$bc.var)), 2))
```

    ##                est stde bc.stde
    ## g.0           0.37 0.02    0.02
    ## g.1           0.02 0.11    0.11
    ## g.2           0.23 0.02    0.02
    ## g.3           0.13 0.08    0.08
    ## g.4           0.23 0.16    0.16
    ## g.5           0.73 0.54    0.54
    ## coef.a.l.1    0.10 0.01    0.01
    ## coef.a.l.2   -0.10 0.01    0.01
    ## coef.a.l.3    0.27 0.03    0.03
    ## log.var.a.l1 -3.17 0.05    0.05
    ## mean.a.1      0.05 0.01    0.01
    ## cov.a.1       0.05 0.00    0.00
    ## EYa.1         0.41 0.05    0.05
    ## EYa.2         0.42 0.04    0.04
    ## EYa.3         0.44 0.03    0.03
    ## EYa.4         0.46 0.02    0.02
    ## EYa.5         0.47 0.01    0.01
    ## EYa.6         0.49 0.01    0.01
    ## EYa.7         0.50 0.02    0.02
    ## EYa.8         0.52 0.02    0.02
    ## EYa.9         0.53 0.03    0.03
    ## EYa.10        0.55 0.04    0.04

### Nonlinear Marginal Structural Model

Suppose the data come arise from a nonlinear outcome model. In
particular, we generate:

- single confounder $L_1 \sim U(0,1)$,
- univariate continuous exposure $A$ with $A|L \sim N(L,0.25)$,
- mismeasured exposure $A^* = A + \epsilon$, where
  $\epsilon \sim N(0, 0.05)$,
- normal outcome $Y$ with
  $Y|A,L \sim N(0.25A + 0.5A^2 - 0.5A^3 + L, 0.16)$.

This data generating process leads to a mean potential outcome at $a$ of
$\textrm{E}\{Y(a)\} = 0.5 + 0.25a + 0.5a^2 - 0.5a^3$.

These data are generated and plotted below.

``` r
# simulate data
cov.e <- 0.05                                   # measurement error variance
inv.link <- inv.ident                           # inverse link
d.inv.link <- d.inv.ident                       # deriv of inv link
g <- c(0, 0.25, 0.5, -0.5, 1)                   # outcome model parameters
formula <- "~A + I(A^2) + I(A^3) + L"           # outcome model formula
args <- list(formula = formula,                 # model fitting arguments
             inv.link = inv.link,
             d.inv.link = d.inv.link)
set.seed(seed)
L <- runif(n)                                                  # confounder
A <- rnorm(n, L, sqrt(0.25))
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
Y <- rnorm(n, EY, 0.16)                                        # outcome
Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
dat0 <- data.frame(Y, A, L)                 # oracle data
datstar <- data.frame(Y, Astar, L)          # mismeasured data
colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L")
a <- seq(-1, 2, length = 10)                # grid of exposure values

# plot data
ggplot(NULL, aes(x = A, y = Y)) + 
  geom_point() +
  ggtitle("Cubic MSM Data",
          subtitle = "Using True Exposure Values")
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggplot(NULL, aes(x = Astar, y = Y)) + 
  geom_point() +
  ggtitle("Cubic MSM Data",
          subtitle = "Using Measured Exposure Values")
```

![](README_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

The dose response curve is then estimated below using the g-formula MCCS
method.

``` r
## G-formula for cubic MSM
gfmla.mccs.cubic <- fit.gfmla.mccs(data = datstar, a = a,
                                   args = args, cov.e = cov.e,
                                   B = B, mc.seed = mc.seed)

cbind(est = round(gfmla.mccs.cubic$est, 2),
      stde = round(sqrt(diag(gfmla.mccs.cubic$var)), 2),
      bc.stde = round(sqrt(diag(gfmla.mccs.cubic$bc.var)), 2))
```

    ##          est stde bc.stde
    ## g.0    -0.03 0.02    0.02
    ## g.1     0.28 0.03    0.04
    ## g.2     0.66 0.07    0.08
    ## g.3    -0.64 0.05    0.06
    ## g.4     1.01 0.02    0.02
    ## EYa.1   1.49 0.11    0.12
    ## EYa.2   0.77 0.04    0.04
    ## EYa.3   0.47 0.02    0.02
    ## EYa.4   0.47 0.02    0.02
    ## EYa.5   0.61 0.01    0.01
    ## EYa.6   0.76 0.01    0.02
    ## EYa.7   0.77 0.02    0.02
    ## EYa.8   0.50 0.02    0.02
    ## EYa.9  -0.20 0.07    0.08
    ## EYa.10 -1.45 0.17    0.19

### Case Cohort Sampling

Consider the same data generating process as in the above example for
the DR estimator, but now using case cohort sampling. That is, a
proportion proportion $\pi=0.25$ of cases are selected to have exposure
$\pmb{A}$ and covariates $\pmb{L}$ measured, as indicated by having
`R=1`. For other cases (with `R=0`), no exposure or covariates are
observed.

``` r
## set parameters
n = 2000                                      # sample size
pi.cc <- 0.25                                 # sampling proportion
seed = 1                                      # random number seed
mc.seed <- 123                                # MCCS seed
cov.e <- 0.02                                 # var(epsilon)
inv.link <- inv.ident                         # inverse link
d.inv.link <- d.inv.ident                     # deriv of inv link
g <- c(0.35, 0.15, 0.25, 0.2, 0.05, 0.1)      # outcome model parameters
formula <- "~A*L1 + A*L2"                     # outcome model formula
ps.formula <- "~L1 + L2"                      # propensity score model formula
ipw.formula <- "~A"                           # ipw.formula

## simulate data
set.seed(seed)
L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
L2 <- rnorm(n, 0, 0.16)                                  # confounder 2
EA <- 0.1 - 0.1*L1 + 0.3*L2                              # E(A|L)
A <- rnorm(n, EA, sqrt(0.04))                            # exposure
Astar <- A + rnorm(n, 0, sqrt(cov.e))                    # mismeasured A
a <- seq(-0.5, 0.5, length.out = 10)                     # grid of exposures
EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
EY[EY < 0] <- 0; EY[EY > 1] <- 1                         # constrain EY
Y <- rbinom(n, 1, EY)                                    # outcome
R <- rbinom(n, 1, pi.cc)                                 # c-c sampling
A[R == 0 & Y == 0] <-
  Astar[R == 0 & Y == 0] <- NA
dat0 <- data.frame(Y, A, L1, L2, R)                      # oracle data
datstar <- data.frame(Y, A = Astar, L1, L2, R)           # measured data
args <- list(formula = formula,                          # arguments for fitting
             ps.formula = ps.formula,
             inv.link = inv.link,
             d.inv.link = d.inv.link)
```

To account for case cohort sampling, we add a sampling weight variable
`cc.wts` to the data set, given by
$\frac{(1-Y)R}{\widehat{pi}_{cc}} + Y$, where
$\widehat{pi}_{cc} = \sum_i \frac{R_i(1-Y_i)}{1-Y_i}$.

``` r
## estimate case-cohort weights
pi.cc.hat <- mean(datstar$R[datstar$Y == 0])
datstar$cc.wts <- (1 - datstar$Y) * datstar$R / pi.cc.hat + Y

## corrected doubly robust estimator
dr.cc.mccs <- fit.dr.mccs(data = datstar, args = args, a = a,
                          cov.e = cov.e, B = B, mc.seed = mc.seed)

cbind(est = round(dr.cc.mccs$est, 2),
      stde = round(sqrt(diag(dr.cc.mccs$var)), 2),
      bc.stde = round(sqrt(diag(dr.cc.mccs$bc.var)), 2))
```

    ##                est stde bc.stde
    ## g.0           0.35 0.03    0.03
    ## g.1           0.13 0.14    0.15
    ## g.2           0.24 0.04    0.04
    ## g.3           0.19 0.14    0.15
    ## g.4           0.03 0.27    0.31
    ## g.5          -0.02 0.90    0.94
    ## coef.a.l.1    0.11 0.01    0.01
    ## coef.a.l.2   -0.12 0.02    0.02
    ## coef.a.l.3    0.26 0.06    0.06
    ## log.var.a.l1 -3.24 0.09    0.09
    ## mean.a.1      0.05 0.01    0.01
    ## cov.a.1       0.04 0.00    0.00
    ## EYa.1         0.40 0.08    0.08
    ## EYa.2         0.42 0.06    0.07
    ## EYa.3         0.44 0.05    0.05
    ## EYa.4         0.45 0.03    0.04
    ## EYa.5         0.47 0.02    0.02
    ## EYa.6         0.49 0.02    0.02
    ## EYa.7         0.50 0.03    0.03
    ## EYa.8         0.52 0.04    0.04
    ## EYa.9         0.53 0.06    0.06
    ## EYa.10        0.55 0.07    0.07

## Simulation Study

Code used for the simulation study in the accompanying paper can be
found in the `simulation` folder. Within this folder are the following
sub directories:

- `sim_scripts`: R scripts to run simulations in parallel on a computing
  cluster,
- `sim_data`: simulation results produced by scripts in sim_scripts,
- `sim_analysis`: Rmd files to analyze the results in `sim_data`,
- `sim_figures`: figures produced in `sim_analysis`.

## Application to HVTN 505 Trial Data

Code for the data application in the accompanying paper can be found in
the `application` folder. The application code uses publicly available
files from the HVTN 505 trial. To get the data files (called at the top
of “application_revised.R”), first go to Start Page: /HVTN Public
Data/HVTN 505 (scharp.org). Next, navigate to the folder “correlates
analysis” and download the tar.gz located there. After download, this
can be installed locally by running “R CMD INSTALL
HVTN505_2019-08-08.tar.gz” in the Terminal or by running
‘devtools::install_local(“HVTN505_2019-4-25.tar.gz”)’. This will include
the necessary data object “dat.505.rda”. Then navigate back to the
parent folder (same url as before) and into the Fong et al. folder,
where the other data file “primary505_for_sharing_upd.csv” can be
downloaded directly.
