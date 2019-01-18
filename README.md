# soilm_global

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2543324.svg)](https://doi.org/10.5281/zenodo.2543324)

This repository contains all R code to reproduce the analysis used for the paper *Satellite monitoring underestimates the impact of drought on terrestrial primary production*. 

A detailed description of the scope, method, and algorithms is given by the RMarkdown file `si_soilm_global.Rmd`. 

To get this repository, change into a suitable directory and clone by
```bash
cd <suitable-working-direcotry>
git clone https://github.com/stineb/soilm_global.git .
```

To render the RMarkdown file, execute all code, and thereby fully reproduce all analysis, simply enter the following command in RStudio (after making sure the 'rmarkdown' package is installed):
```r
install.packages("rmarkdown", type = "source")
rmarkdown::render_site()
```

Please also note the data use policy, described in `./si_soilm_global.Rmd` and `./LICENSE`. When using this code, please cite the paper Stocker et al. (2019) *Nature Geoscience*. This code is released with a DOI on [Zenodo](http://doi.org/10.5281/zenodo.1286966), where the latest version corresponds to tag `v1.0`. Model outputs from site-scale and global simulations are available on Zenodo with doi:10.5281/zenodo.1423484.

beni stocker, 18.1.2019
