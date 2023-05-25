# "classmapExt"" R package
This R package called "classmapExt" provides supplementary functions designed to extend the functionalities of the ["classmap"](https://cran.r-project.org/web/packages/classmap/index.html) package. ["classmap"](https://cran.r-project.org/web/packages/classmap/index.html), released on CRAN by Raymaekers J. and Rousseeuw P., provides useful and novel graphical displays to visualize results of classification cases (with known class).
The package here: 
- adds support for the visualizations on a classification made by a particular classification algorithm of choice via `vcr.custom.*` .
- adds support for the visualizations for the ["pamr"](https://cran.r-project.org/web/packages/pamr/index.html) (Nearest Shrunken Centroid) ML classifier via `vcr.pamr.*` .
- adds an additional graphical display though `mdscolorscale`. 

This package is intended to be loaded together with the classmap package.

To use this package  in R install using `devtools::install_github("LLazzar/classmapExt")` and load the installed library `library("classmapExt")`

The code available in this package comes from previous work available in [this](https://github.com/LLazzar/classmap-on-pamr) repo and [this](https://github.com/LLazzar/extension-classmap) other one
