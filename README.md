
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

## Example

An example of the three proposed estimators (g-formula, IPW, and doubly
robust) used on a simulated data set is provided below.

### Data Generation

Data are generated below according to the data generation process
described in the third simulation study:

- sample size $n=2000$,
- two independent confounders
  $L_1 \sim \textrm{Bernoulli}(0.5), L_2 \sim N(1, 0.5)$
- univariate exposure $A$, where
  $A|\pmb{L} \sim N(2 + 0.9L_1 - 0.6L2, 1.1)$
- continuous outcome $Y$, with
  $Y|\pmb{L},A \sim N(1.5 + 0.7A + 0.9L_1 - 0.7L_2 + - 0.6AL_1 + 0.4AL_2, 0.16)$
- resulting MSM $\textrm{E}\{Y(a)\} = \gamma_0 + \gamma_1a$, where
  $\pmb{\gamma} = (\gamma_0, \gamma_1) = (1.35, 0.75)$
- mismeasured exposure $A^* = A + \epsilon$, where
  $\epsilon \sim N(0, 0.16)$

``` r
## define parameters
 
n = 2000                                  # sample size
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
datstar <- data.frame(Y, A = Astar, L1, L2)              # mismeasured data

head(datstar, 5)
```

    ##           Y          A L1        L2
    ## 1 2.2186524  1.1634983  0 1.8025415
    ## 2 2.5664615  1.7145178  0 1.7862545
    ## 3 2.1174644  4.3731142  1 0.3842672
    ## 4 2.8111597  1.7921312  1 1.1490097
    ## 5 0.9702103 -0.4505236  0 1.0490701

### G-Formula Estimation

The dose response curve at chosen values of the exposure $a$ are
estimated here using the MCCS g-formula method.

Before fitting the model, we determine an appropriate number of
Monte-Carlo replicates $B$ in order for the G-formula MCCS function to
approximate the CS function. We can do this by evaluating the MCCS
function for a sequence of $B$ values, and at a particular parameter
value, say the naive g-formula estimator (ignoring measurement error).

``` r
## g-formula arguments
gfmla.args <- list(formula = formula,   
                   inv.link = inv.link,
                   d.inv.link = d.inv.link)

## naive estimator
gfmla.naive <- fit.glm(data = datstar,
                       args = gfmla.args,
                       return.var = F)$est

## grid of possible B values
B.grid <- seq(1, 50, by = 2)

## store psi and computation time B
B.search <- vapply(
  X = 1:length(B.grid),
  FUN.VALUE = numeric(8),
  FUN = function(ii) {

    st <- Sys.time()
    get.psi.glm.mccs <- make.mccs(
      get.psi = get.psi.glm, data = datstar, args = gfmla.args,
      cov.e = cov.e, B = B.grid[ii], mc.seed = mc.seed)
    psi <- get.psi.glm.mccs(x = gfmla.naive)
    et <- Sys.time()

    return(c(B = B.grid[ii],
             Time = et - st,
             psi = psi))
  }) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(c("B", "Time", paste0("psi", 0:5))) %>% 
  pivot_longer(cols = !B)

## plot results
ggplot(data = B.search,
       aes(x = B,
           y = value)) +
  geom_line() +
  facet_wrap(~ name,
             scales = "free") +
  labs(y = "") +
  ggtitle("Score Values and Computation Time by Number of MC Replicates B")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Based on these plots, the G-formula MCCS seems to stabilize around
$B=30$. Using this value, we proceed with estimating the G-formula
parameters. These include the parameters
$\pmb{\beta} = (\beta_0, \dots, \beta_5)$ in the outcome model
$Y|\pmb{L},A$, and the dose response curve $\textrm{E}\{Y(a)\}$
evaluated at four points (-2.5, 0.3, 3, 5.8) in the support of $A$.

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

Parameters $\pmb{\gamma} = (\gamma_0, \gamma_1)$ in the marginal
structural model, as well as coefficients and variance in the propensity
score model $A|\pmb{L}$, are estimated here using the MCCS IPW method.
We use the same number of MC replicates $B=30$ as for IPW, but a similar
strategy as with g-formula could be used to tune this to an appropriate
number.

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
      stde = round(sqrt(diag(ipw.res$var)), 2))
```

    ##                est stde
    ## g.0           1.77 0.19
    ## g.1           0.55 0.09
    ## coef.a.l.1    2.06 0.05
    ## coef.a.l.2    0.88 0.05
    ## coef.a.l.3   -0.66 0.04
    ## log.var.a.l1  0.13 0.04

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

### Using Known Weights

If the propensity weights are known, or if they are estimated
externally, then they can be supplied using the `ps.wts` argument. In
this case, the sandwich variance estimator treats weights as known.

``` r
## Known PS weights
ps.wts <- dnorm(A, mean(A), sd(A)) /
  dnorm(A, 2 + 0.9*L1 - 0.6*L2, sqrt(1.1))
  
## Double Robust (known weights)
dr.res.known.wts <- fit.dr.mccs(
  data = datstar,
  a = a,
  cov.e = cov.e,
  B = B,
  mc.seed = mc.seed,
  return.var = TRUE,
  args = dr.args,
  ps.wts = ps.wts)

cbind(est = round(dr.res.known.wts$est, 2),
      stde = round(sqrt(diag(dr.res.known.wts$var)), 2))
```

    ##         est stde
    ## g.0    1.53 0.04
    ## g.1    0.69 0.02
    ## g.2    0.87 0.05
    ## g.3   -0.58 0.03
    ## g.4   -0.70 0.02
    ## g.5    0.39 0.01
    ## EYa.1 -0.47 0.07
    ## EYa.2  1.58 0.02
    ## EYa.3  3.62 0.03
    ## EYa.4  5.66 0.07
