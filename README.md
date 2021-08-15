
The profoc Package
======================

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/github/workflow/status/berrij/profoc/R-CMD-check?style=for-the-badge)](https://github.com/BerriJ/profoc/actions/workflows/R-CMD-check.yaml)
[![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/berrij/profoc/pkgdown/main?label=Documentation&style=for-the-badge)](https://profoc.berrisch.biz/)
[![Lifecycle: experimental](https://img.shields.io/badge/Lifecycle-experimental-orange?style=for-the-badge)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->


The main function `online` can be used to combine probabilistic forecasts using the CRPS learning algorithm proposed in Berrisch, Ziel (2021).
The function `batch` can be used in a similar way for batch optimization.
Print, plot, update, and predict methods are available.

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

## Contributions and Issues

Feel free to [raise an issue](https://github.com/BerriJ/profoc/issues/new) if you find something not working properly.

You are also very welcome to contribute to profoc. Please base your pull requests on the develop branch.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
