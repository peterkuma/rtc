Tree Clustering Library (R package)
===================================

RTC performs cluster analysis using trees. A tree approximates
the probability density function from which data has been drawn.
A set of such trees is inferred using a Metropolis-Hastings sampler,
and returned as a partitioning of the parameter space into rectangular segments.

This package uses the [libtc](https://bitbucket.org/pkuma/libtc) library.

* Depends: R (>= 2.15.1)
* License: Proprietary
* Reference manual: [rtc.pdf](/pkuma/rtc/src/master/rtc.pdf)

Installation
------------

Requirements:

* [libtc](https://bitbucket.org/pkuma/libtc)

The package can be installed using the
[devtools](http://www.rstudio.com/products/rpackages/devtools/) package:

    Rscript -e 'install.packages("devtools", repos="http://cran.rstudio.com")'
    Rscript -e 'library(devtools); install_bitbucket("pkuma/rtc")'

or manually:

    git clone https://pkuma@bitbucket.org/pkuma/rtc.git
    R CMD INSTALL rtc

See the reference manual for information about usage.
