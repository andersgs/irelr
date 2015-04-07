# Welcome

Welcome the R package irelr. This package provides functionality formerly available
in the software iRel (Gonçalves da Silva and Russello 2011). It adds to the functionality
by providing the user with additional relatedness indices to those in the original 
publication. 

# How to get it

The package contains C code, so if you are in Windows, you will require the 
RTools package.

I recommend installing the latest version of R from [here](http://www.r-project.org/).
Then installing the latest version of RStudio from [here](http://www.rstudio.com/products/rstudio/download/).

Everyone will require the Hadley's devtools package. From the R command prompt within RStudio, type:

    install.packages('devtools')

Once that is complete, load `devtools`:

    library(devtools)

You can then install `irelr` with the following command:

    install_github("andersgs/irelr", build_vignettes = TRUE)

This might take a few moments as the vignette and source files are built.

# How to use it

The package provides two main functions. A function to estimate relatedness indices
from a set of genotypes, and a function to simulate indices under different 
relatedness categories.

The easiest way to get on your way with `irelr` is to have you data in a 
`genepop` format, with individual identifiers for each sample. `irelr`
ignores population structure. The data file must have a `.gen` extension.
To load your own `genepop` file:

    library(adegenet)
    file_path  <- "<path to file>/mydata.gen"
    mydata <- read.genepop(file_path)

To estimate relatedness values:

    library(adegenet)
    data(nancycats)
    irelr::estimate_rel(nancycats)
    
To simulate relatedness values for the available indices, one needs to define a 
`k`-vector (explained in the documentation and vignette). To simulated 
indices from 10,000 unrelated pairs:

    library(adegenet)
    data(nancycats)
    irelr:sim_rel(nancycats, k_vector = c(1.0, 0.0, 0.0))

An extensive vignette is available with details on how to use the results from
these two functions to reconstruct a pedigree. To access the vignette just type:

    vignette(topic = 'use-irelr', package = 'irelr')

# How to cite

Please cite the original `iRel` publication below, and this code specifically 
with the following doi: https://zenodo.org/badge/doi/10.5281/zenodo.16651.svg

# References

Gonçalves da Silva A, Russello MA (2011) iRel: software for implementing pairwise relatedness 
estimators and evaluating their performance. Conservation Genetics Resources 3: 69-71.
[PDF](http://link.springer.com/article/10.1007/s12686-010-9292-4)

# License

irelr: an R package to reconstruct pedigrees from molecular data

Copyright (C) 2015  Anders Gonçalves da Silva and Michael A. Russello

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.