# Reproduction of BasisGraphicalLasso

This package contains code to perform the basis graphical lasso
analysis of:

>Krock, M., Kleiber, W., and Becker, S. (2021), “Nonstationary modeling with sparsity for spatial data via the basis graphical lasso,” J. Comput. Graph. Statist., 30, 375–389, ISSN 1061-8600, URL https://doi.org/10.1080/10618600.2020.1811103

Please install the package with

```r
devtools::install_github("mlkrock/BasisGraphicalLasso")
```

This implementation relies on the `QUIC` function for solving the graphical lasso. The CRAN package is orphaned, so install the latest version of QUIC from the Archive in CRAN.

You need to install the following dependencies also:
```r
install.packages('sp')            #To convert a dataframe to an sptial object
install.packages('LatticeKrig')   #To construct a particular type of basis functions
```
You have to use the latest version of R(4.3) and Rtools(4.3) to install these libraries.

To install the stable/testing version of INLA package in R (will be needed to install FRK package):
```r
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages(FRK)
```
```r
install.packages('psych')   # To determine trace of a matrix
install.packages('sf')
install.packages('huge')    # To create the Hub Graph/ Graphical structure of the precision matrix
install.packages('glasso')  # To perform Graphical lasso algorithm. 
```
The basic function is BGL, and the analysis of the minimum temperature
dataset in the paper can be reproduced as

```r
precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
     distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
     MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)
```
