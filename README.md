# PORTO: Preference and ORdinal Threshold Orientations

**Developer and maintainer:** Thyago Nepomuceno, Department of Statistics,
Federal University of Pernambuco, Brazil. Email: nepomuceno@de.ufpe.be.

PORTO implements ordinal weight restrictions for radial DEA and Benefit-of-the-Doubt (BoD)
composite indicators, as well as directional BoD for desirable and undesirable indicators.
Version 0.1.0 is an initial source-package release, not a CRAN publication.

## Install from GitHub

Create a repository (e.g. `PORTO`) and upload **the contents of this folder** to
its root, including `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, `tests/`, and the
hidden `.github/` workflow. Then, in R:

```r
install.packages("remotes")
remotes::install_github("usuario/PORTO") # replace usuario with your GitHub login
library(PORTO)
help(package = "PORTO")
```

The sole non-base dependency is `lpSolve`. For local installation, use
`remotes::install_local("path/to/PORTO")` or `R CMD INSTALL PORTO`.

## Functions

| Function | Data | Rank ties |
|---|---|---|
| `dea.order()` | Radial DEA or unit-input/unit-output BoD | No |
| `dea.order.ties()` | Radial DEA or unit-input/unit-output BoD | Yes |
| `dir.dea.order()` | Directional BoD (good and bad indicators) | No |
| `dir.dea.order.ties()` | Directional BoD (good and bad indicators) | Yes |

All functions use `Weak = TRUE` by default and `Weak = FALSE` for strict
preferences. Optionally use `Weak = TRUE, Gamma = 0.5` for moderate preferences.
`Gamma=0` is weak, `Gamma=1` is strict. The `*.ties()` routines use `TIES="pooled"`.
All ranks are oriented with smaller numbers indicating higher priority; rank
matrices must correspond exactly to the data matrices by row and column.

```r
Y <- matrix(c(5, 8, 8, 5, 6, 6), ncol = 2, byrow = TRUE)
OY <- matrix(rep(c(1, 2), 3), ncol = 2, byrow = TRUE)
weak <- dea.order(Y = Y, OY = OY)
strict <- dea.order(Y = Y, OY = OY, Weak = FALSE)
moderate <- dea.order(Y = Y, OY = OY, Gamma = 0.5)
cbind(weak$efficiency, moderate$efficiency, strict$efficiency)
strict$virtual_output

B <- matrix(c(4, 3, 5), ncol = 1)
G <- Y[, 1, drop = FALSE]
fit <- dir.dea.order(Y = G, B = B,
                     OY = matrix(1, 3, 1), OB = matrix(2, 3, 1))
fit$efficiency
fit$solver_status
```

## Documentation and checking

Each function has its own `man/*.Rd` page. From the parent directory:

```sh
R CMD build PORTO
R CMD check --as-cran PORTO_0.1.0.tar.gz
R CMD Rd2pdf --output=PORTO-reference.pdf PORTO
```

The GitHub Actions workflow runs package checks on Linux, Windows and macOS.
This environment did not have R installed; the tests are supplied for execution
on a real R installation. No passing R CMD check is asserted. The license field
is provisionally GPL-3; confirm it before publishing. See the reference manual
for details of the LPs, orientation naming, tie pooling, zero handling and
limitations. No benchmark dataset is bundled, and no earlier SDG results are
claimed to have been recomputed in this packaging step.
