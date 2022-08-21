
# readsas

<!-- badges: start -->

![R-CMD-check](https://github.com/JanMarvin/readspss/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

R package using Rcpp to parse a SAS file into a data.frame(). Currently
`read.sas` is the main function and feature of this package.

The package allows (experimental) reading of sas7bdat files that are

- (un)compressed

As with other releases of the `read` series, focus is again on being as
accurate as possible. Speed is welcome, but a secondary goal.

## Installation

With `drat`:

``` r
drat::addRepo("JanMarvin")
install.packages("readsas")
```

With `devtools`:

``` r
devtools::install_git("https://github.com/JanMarvin/readsas.git")
```

## Basic usage

``` r
fl <- system.file("extdata", "cars.sas7bdat", package="readsas")

dd <- read.sas(fl)

head(dd)
#>   speed dist
#> 1     4    2
#> 2     4   10
#> 3     7    4
#> 4     7   22
#> 5     8   16
#> 6     9   10
```

## Select columns or rows

This is not necessarily faster (only if less than the maximum number of
rows is selected). selected), but it is memory efficient to load only
certain columns or rows.

``` r
fl <- system.file("extdata", "mtcars.sas7bdat", package="readsas")

dd <- read.sas(fl, select.cols = c("VAR1", "mpg", "hp"),
               select.rows = c(2,5), rownames = TRUE)

head(dd)
#>                    mpg  hp
#> Mazda RX4 Wag     21.0 110
#> Datsun 710        22.8  93
#> Hornet 4 Drive    21.4 110
#> Hornet Sportabout 18.7 175
```

## Thanks

The documentation of the sas7bdat package by Matt Shotwell and Clint
Cummins in their R package
[`sas7bdat`](https://github.com/BioStatMatt/sas7bdat), by Jared Hobbs
for the python library
[`sas7bdat`](https://bitbucket.org/jaredhobbs/sas7bdat/src/master/), and
by EPAM in the Java library [`parso`](https://github.com/epam/parso) was
crucial. Without their decryption of the SAS format, this package would
not have been possible.

Further testing was done using the R package
[`haven`](https://github.com/tidyverse/haven) by Hadley Wickam and Evan
Miller.
