## Test environments

* GitHub Actions - (ubuntu-20.04): release, devel
* GitHub Actions (windows): release
* Github Actions (macOS): release

## R CMD check results

❯ checking installed package size ... **NOTE**

    installed size is 29.0Mb
    sub-directories of 1Mb or more:
      libs  28.4Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Note on the size

The first Note is due to Rcpp, and we cannot solve this because it is caused by the way Rcpp handles header-only libraries.

## Downstream dependencies

There are currently no downstream dependencies for this package.