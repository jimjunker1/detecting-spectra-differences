# detecting-spectra-differences

This repository holds the code to reproduce the simulations, results and plots presented in Pomeranz et al. ([Preprint:](https://www.biorxiv.org/content/10.1101/2023.03.14.532592v1), in reviison). 

## Setup 

This repository relies on the `sizeSpectra` package by Andrew Edwards. Instructions on how to download this package from GitHub are available [here](https://github.com/andrew-edwards/sizeSpectra). 

It also relies on a number of bespoke functions written specifically for this analysis, and they are in the `custom_functions.R` script, and are loaded into the appropriate scripts by running the `source("custom_functions.R")` line of code. A full description of the custome functions can be found below. 

Finally, this repository also uses the `tidyverse`, `tidybayes`, and `cowplot` functions, which can all be installed from CRAN using the `install.packages("package-name")` code. 

## Running the code

The scripts are numbered to be run in order. 

`1a_sim_sliding_windows.R`: This is the main analysis presented in the manuscript. It samples observations of body sizes (n=999) from $\lambda$'s ranging from (-0.5, -2.5) and estimates site-specific values of the relationship. This is repeated a number of times to get a distribution of the site-specific $\lambda$ estimates. 
It then samples body size observations from $\lambda$ values varying by $\beta_{env}) = -0.5$ units across a hypothetical gradient, and then tries to recapture the known relationship. This is repeated across three different "windows" (see manuscript) of $\lambda$ values. Likewise, different relationships are also examined ($\beta_{env}) = (-0.25, 0)$

`1b_simulations_sensitivity.R`: This script repeats the main simulation but varies some aspects, such as the number of sites and range of values for the hypothetical gradient, the range and sample size of body size observations. 

`2b_integrated_plots.R`: This script loads the simulated datasets and produces multi-panel, integrated plots for the main text. 

`3_bias.R`: This script estimates the bias in site-specific $\lambda$ values and $\beta_{env}$ estimates across the three methods. 

`4_a_AMD_estimate`: This script loads the [previously published acid mine drainage data](https://onlinelibrary.wiley.com/doi/abs/10.1111/fwb.13196) and estimates the relationship of $\lambda$ using the three methods.  The data used in the present analysis can be found in the `site_dw_pca.csv` file. 

`4_a_NEON_estimate`: This script loads the National Ecological Observatory Network data analyzed in [Pomeranz et al. 2022](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.15862) and estimates the relaiionship of $\lambda$ using the three methods. The data used in the present analysis can be found in the `site_dw_pca.csv` file. 

`4c_empirical_plot.R`: This script takes the estimates from scripts 4a and 4b and makes a single, multi-panel plot for the main text of the manuscript. 

`5a_SI_figs.R`: This script makes various remaining plots for the Supplemental Information section. Not all the plots made it into the final document. 

## Aditional files  

The first draft of this paper was written using RMarkdown. The `Manuscript.Rmd` file is the main document, and the `Manuscript` files with various extensions (i.e., .pdf, .html, .docx, etc.) were the various "knits" of the document made during the writeup process. Likewise, the `author-info-blocks.lua`, `ecology.csl`, `detecting-size-spectra-difference.bib`, etc., are all files necessary for the inclusion of references, citation style, etc. The files are retained here for posterity, and for those interested in seeing an example of integrating RMarkdown in the writing of a manuscript.
* After the first round of reviews, we switched to making changes and edits in Microsoft Word using the "track changes" feature
* The `Manuscript` files here do not reflect changes made to future versions of the submitted manuscript.
* Upon article acceptance, a link to the final version will be included here.
* We have made separate `.rmd` files for the figures and for the supplemental information. 
  
## Custom Functions 

Below is a general description of the custom functions written for this analysis. Please see the `custom_functions.R` script for further details including comments and function annotations.  

`MLE_tidy`: This is a modified function from the `sizeSpectra` package that is compatible with "tidy" data. The data need to include a "grouping" variable.  

`bin_and_center`: This function takes body size data and "bins" it (i.e., the ELBn, and L2n methods in the main text). It takes an argument `breaks` which is a vector of numbers delineating the edges of the bin breaks. `breaks` are calculated for log~2~ bins and equal logarthmic bins (n = 6) in the `compare_slopes()` function (see below)  

`est_ci`: Thgis is a wrapper function to extract the confidence intervals from a lm (linear model) output. It is used internally in the `compare_slopes()` function

`compare_slopes`: This is a wrapper for the main analysis. It estimates $\lambda$'s for a data set using each of the three methods presented in the main text. The "L2n" method was originally reffered to as "NAS", and the internal variables calculated in this function retain the original "NAS" designation. 

`sim_result`: This is a wrapper function which replicates the following process:
1. samples body size data from a bounded power law distribution with a set $\lambda$ value
2. Estimates $\lambda$ values based on the three methods
3. Returns a data.frame with relevant values including the replicate, known and estimated $\lambda$, body size range, hypothetical gradient values, etc.  

`plot_sim`: this function is a wrapper which takes the simulated data, plots the distribution of site-specific $\lambda$ values, estimates relationship across the hypothetical gradient and plots their distribution, and writes out the results in a .csv file.   

`lambda_estimate_density_plot`: This is a function to plot the distributions of site-specific $\lambda$ estimates from the simulated data sets. It is used internally in the `plot_sim()` function  

`calc_relationship_estimate`: This is a function to estimate the relationship of $\lambda$ estimates across the hypothetical gradient from the simulated data sets. It is used internally in the `plot_sim()` function  

`relationship_density_plot`: This is a function to plot the ditribution of the estimated the relationship of $\lambda$ estimates across the hypothetical gradient from the simulated data sets. It is used internally in the `plot_sim()` function  

`main_plot`: This plots the individual $\lambda$ estimates across the hypothetical gradient, and includes an individual regression estimate (`stat_smooth(method = "lm")`) for each replicate in the simulated data. It is used internally in the `plot_sim()` function
