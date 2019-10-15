# Readsas [![Build Status](https://travis-ci.org/JanMarvin/readsas.svg?branch=master)](https://travis-ci.org/JanMarvin/readsas) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/JanMarvin/readsas?branch=master&svg=true)](https://ci.appveyor.com/project/JanMarvin/readsas)

R package using Rcpp to parse an SAS file into a data.frame(). Currently 
`read.sas` is the main function and feature of this package.

The package features (experimental) reading of sas7bdat files that are

* uncompressed
* char compressed

As with other releases of the `read` series, focus is again on being as 
exactly as possible. Speed is welcome, but a second seat passenger.

## Installation

With `devtools`:
```R
devtools::install_git("https://github.com/JanMarvin/readsas.git")
```
## Usage
```R
fl <- system.file("extdata", "cars.sas7bdat", package="readsas")

dd <- read.sas(fl)
```


## Thanks

The documentation of the sas7bdat package by Matt Shotwell and Clint Cummins in
their R package [`sas7bdat`](https://github.com/BioStatMatt/sas7bdat) and by 
Jared Hobbs for the python library 
[`sas7bdat`](https://bitbucket.org/jaredhobbs/sas7bdat/src/master/) was crucial!

Additional testing was done using the R package 
[`haven`](https://github.com/tidyverse/haven) by Hadley Wickam and Evan Miller.
