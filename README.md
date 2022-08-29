# csiro22

Course material for the CKMR workshop at CSIRO Battery Point (Hobart), Aug 29th-Sept 2nd 2022.

In the course of the workshop you will learn:

 - How to build and fit basic CKMR models for several situations; 
 - What CKMR can (and sometimes cannot) tell you; 
 - How to avoid modelling mistakes; 
 - A little about modern genetics and using it to find your close-kin; 
 - How to go about Design for a new project.
 

General documents: 

 - [Rcrip22.pdf](https://github.com/markbravington/csiro22/blob/main/Rcrib22.pdf): an overview of the required software and coding framework used in the workshop
 - [turn_off_bc.r](https://github.com/markbravington/csiro22/blob/main/turn_off_bc.r): script to source once before using the offarray package (turns off R's default byte-compiler)
 
Exercises for specific workshop days will be added throughout.

There is a [workshop website](https://markbravington.github.io/csiro22/tutorials/about.html), on which
tips and vignettes for some of the exercises will be posted during the week.


## Day 1 exercise

* Rerun the delfi A example (`data_and_scripts/fit_delfi_A.r`) with a third of the samples

Tip: You'll need to remove two thirds of all the samples, including POPs. You'll also need to adjust
the row numbers in the POPs dataset to match the new row numbers in your `Samps` table. 
You can use the function `subset_samples` to do this for you, e.g. 

`newsamps <- subset_samples(samp_delfi_A, runif(nrow(samp_delfi_A$Samps)) < 1/3)`

Then, move on to `boring_data_prep_delfi_A` using `newsamps` instead of `samp_delfi_A`.
