# Welcome

Welcome the R package irelr. This package provides functionality formerly available
in the software iRel (Gonçalves da Silva and Russello 2011). It adds to the functionality
by providing the user with additional relatedness indices to those in the original 
publication. 

# How to get it on Linux/Mac

I recommend installing the latest version of R from [here](https://www.R-project.org).
Then installing the latest version of RStudio from [here](http://www.rstudio.com/products/rstudio/download/).

You require Hadley Wickham's devtools package. From the R command prompt within RStudio, type:

    install.packages('devtools')

Once that is complete, load `devtools`:

    library(devtools)

You can then install `irelr` with the following command:

    install_github("andersgs/irelr", build_vignettes = TRUE)

This might take a few moments as the vignette and source files are built.

# How to get it on Windows

The latest release (version 0.0.6) has binaries available for `R` version `3.3.2`.

To install: 

 0. Install dependencies: `install.packages(c("adegenet","ape","data.table","ggplot2","gridExtra","hierfstat","moments","reshape2"))`
 1. Download the [file](https://github.com/andersgs/irelr/releases/download/v0.0.6/irelr_0.0.6_win.zip)
 2. Open `R`, click on `Packages` -> `Install package(s) from local files` (If in `RStudio` select `Tools` -> `Install Packages..` --- then select `Install from: Package Archive File`)
 3. Select the file `irelr_0.0.6_win.zip` that you downloaded in step 1.
 
 If all goes well, you should see the following message on the `R` command-prompt:
 
 `package irelr successsfully unpacked and MD5 sums checked`
 
 You should then be able to load `irelr` with the following:
 
 `library("irelr")`
 
 If all works out, you should see the following printed to the console:
 
 `Welcome to irelr (version 0.0.6)`
 

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
with the following: [![DOI](https://zenodo.org/badge/33524849.svg)](https://zenodo.org/badge/latestdoi/33524849)

# References

Gonçalves da Silva A, Russello MA (2011) iRel: software for implementing pairwise relatedness 
estimators and evaluating their performance. Conservation Genetics Resources 3: 69-71.
[PDF](http://link.springer.com/article/10.1007/s12686-010-9292-4)

# Updates

## 2017-01-29
  * Added binary files for Windows for release 0.0.6 (built on version 3.3.2 of `R`) using [win-builder.r-project.org](https://win-builder.r-project.org)
  * Added instructions on how to install the Windows version

## 2015-12-05
  * Updated the vignette to be compliant with the new version of `Genind` object
  taking advantage of the accessor functions.
  * `irelr` now requires at least `adegenet` version 2.0.0
  * Added a function `read.genepop` to replace the `adegenet` version 2.0.0 
  function which does not parse missing data correctly. The included version
  is a copy/paste of the function in the development version of `adegenet`
  found on [here](https://github.com/thibautjombart/adegenet)
  

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