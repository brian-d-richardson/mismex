
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
```

The `mismex` package contains functions to estimate causal effects in
the presence of confounding and a mismeasured exposure. The methods
implemented are introduced in the paper, “Addressing confounding and
continuous exposure measurement error using corrected score functions,”
which is currently in progress.
