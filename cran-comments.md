## Test environments

* GitHub Actions - (ubuntu-20.04): release, devel
* GitHub Actions (windows): release
* Github Actions (macOS): release

## R CMD check results

❯ checking installed package size ... **NOTE**

    installed size is 29.1Mb
    sub-directories of 1Mb or more:
      libs  28.5Mb

❯ checking compiled code ... **NOTE**

    File ‘profoc/libs/profoc.so’:
    Found ‘_ZSt4cout’, possibly from ‘std::cout’ (C++)
        Object: ‘oracle.o’
    Found ‘putchar’, possibly from ‘putchar’ (C)
        Object: ‘oracle.o’
    Found ‘puts’, possibly from ‘printf’ (C), ‘puts’ (C)
        Object: ‘oracle.o’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

### Note on the size

The first Note is due to Rcpp, and we cannot solve this because it is caused by the way Rcpp handles header-only libraries.

### Note on C++ and C calls

The second Note is caused by code from the [optim C++ library](https://github.com/kthohr/optim). It's a header-only library that we use for numerical optimization. We are going to contact the maintainer to remove the respective parts, but we hope that you can accept profoc for now despite this Note.

## Downstream dependencies

There are currently no downstream dependencies for this package.