Calcium Bone Model
==================

-   Documentation [here](vignettes/modeldoc.Rmd)

``` r
library(cabone)
```

Simulate teriparatide data
--------------------------

``` r
out <- sim_teri(dose=c(20,40), dur=9)

plot(out)
```

![](inst/img/README-unnamed-chunk-3-1.png)

Simulate denosumab data
-----------------------

``` r
out <- sim_denos(dose=c(10,60,210), dur=6)

plot(out, log(DENCP) + BMDlsDENchange~time)
```

![](inst/img/README-unnamed-chunk-4-1.png)

Simulate secondary hyperparathyroidism
======================================

Have `GFR` decline by a certain amount (`GFRdelta`) over a certail period of time (`GFRtau`).

``` r
sim_2h() %>% plot
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
file_location <- cabone_export(file)
```
