
About
=====

A multiscale systems model of bone health and mineral homeostasis.

Documentation
=============

-   Documentation [here](vignettes/modeldoc.Rmd)

Installation
============

``` r
remotes::install_github("metrumresearchgroup/OpenBoneMin")
```

Usage
=====

-   [Simulate teriparatide data](#teri)
-   [Simulate denosumab data](#denos)

``` r
library(OpenBoneMin)
library(ggplot2)
```

<a name="teri"></a>

Simulate teriparatide data
--------------------------

-   `PTHpm` teriparatide concentration (pM)
-   `CaC` calcium concentration (mM)

``` r
out <- sim_teri(dose=c(20,40), dur=9)

plot(out)
```

![](inst/img/README-unnamed-chunk-4-1.png)

<a name="denos"></a>

### Simulate denosumab data

-   `DENCP` denosumab concentration
-   `BMDlsDENchange` lumbar spine change from basline

``` r
out <- sim_denos(dose=c(10,60,210), dur=6)

plot(out, log(DENCP) + BMDlsDENchange ~ time, xlab="Time (months)")
```

![](inst/img/README-unnamed-chunk-5-1.png)

Some helper functions
=====================

Convert `teriparatide` doses
----------------------------

Usually, we think of doses in micrograms. This function turns those doses into `pmol`.

``` r
amt_teri(20)
```

    . [1] 4856.962

Export the model code
---------------------

It's a little hard to see what's happening here. But basically, this grabs the model code and writes it to a file of your choosing. Use this when you want to export the model and start making changes yourself.

``` r
file <- file.path(tempdir(),"my_model.cpp")
file_location <- BoneMin_export(file)
```
