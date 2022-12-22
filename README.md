
# The profoc Package

## An R package for probabilistic forecast combination
---

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/berrij/profoc/R-CMD-check.yaml?branch=main&style=for-the-badge)](https://github.com/BerriJ/profoc/actions/workflows/R-CMD-check.yaml)
[![GitHub Workflow Status (branch)](https://img.shields.io/github/actions/workflow/status/berrij/profoc/pkgdown.yaml?branch=main&label=Documentation&style=for-the-badge)](https://profoc.berrisch.biz/)
[![Lifecycle: stable](https://img.shields.io/badge/Lifecycle-stable-brightgreen?style=for-the-badge)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The primary function `online` can be used to combine probabilistic forecasts using the CRPS learning algorithm introduced in Berrisch, Ziel (2021): [Pre-Print](https://arxiv.org/pdf/2102.00968.pdf), [Publication](https://doi.org/10.1016/j.jeconom.2021.11.008).
The function `batch` can be used in a similar way for batch optimization. Common methods like `summary`, `print`, `plot`, `update`, and `predict` are available.

Installation
------------

### Install from CRAN

You can install the latest stable release from CRAN using:

``` r
install.packages("profoc")
```

### Install from GitHub

You can install the latest stable release from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("BerriJ/profoc")
```

You can install the latest development release from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("BerriJ/profoc@develop")
```

## Documentation

You can find the documentation at [profoc.berrisch.biz](https://profoc.berrisch.biz/).

## Contributions and Issues

Feel free to [raise an issue](https://github.com/BerriJ/profoc/issues/new) if you find something not working properly.

You are very welcome to contribute to profoc. Please base your pull requests on the develop branch.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
