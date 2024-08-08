Readme
================
2023-04-26

<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesTraitR

<!-- badges: start -->
<!-- badges: end -->

The goal of BayesTraitR is to facilitate phylogenetic comparative
analyses in BayesTraits for R users. It comprises a suite of functions
and tools that will be useful for wrangling datasets and trees into the
appropriate format for analyses within BayesTraits as well as for
interpreting and post-processing the outcomes of such analyses.

**Jump to:**

2)  [Installation](#installation)

3)  [Creating Analyses](#creating-analyses)

4)  [Post-processing](#post-processing)

5)  [Files bundled with BayesTraits](#BTFiles)

## Installation

You can install the development version of BayesTraitR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("joannabaker/BayesTraitR")
```

You can download the latest version of
[BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.1/BayesTraitsV4.0.1.html)
from the University of Reading website. Note that this package has been
developed with the output from BayesTraitsV4.0.1 and different versions
of the program may produce incompatible outputs.

## Creating Analyses

The function is used to generate an input file, data file and tree file
for use in downstream BayesTraits analyses. This can then be run
interactively from within R (using or similar) or externally using the
command-line interface.

For all models that can be run in BayesTraits, we can simply specify a
list of columns found in our dataset that we wish to include using the
argument. However, for regression analyses that may include complex
contrast coding, categorical variables, interactions and transformations
(e.g. polynomials), users may wish to alternatively specify their
analysis variables using the argument.

### Multi-state model

## Post-processing

### The Variable Rates Post-Processor

There is a specific set of commands and functions associated with the
variable rates model and the variable rates regression model. We can use
these functions to import the files output by variable rates analyses
into R and post-process into results ready for interpretation.

There are three types of variable rates models that can be run that each
have very slight differences in the output files.

Here, I will include how to run the variable rates post-processor using
all three examples:

1)  A variable rates model

2)  A variable rates model across a sample of trees

3)  A variable rates model across a tree sample, forcing equal tree
    sampling.

## BayesTraits Files

This package contains the sample files bundled with BayesTraits. This
includes the following trees:

- Artiodactyl_trees (A sample of 500 trees of 17 artiodactyls.)

- Bird_tree (A tree of 5669 birds.)

- Mammal_trees (A sample of 50 trees of 40 mammals.)

- Marsupials_tree (A tree of 246 marsupials.)

- NortheastBantu_tree (A tree of 85 Bantu languages.)

- Primates_trees (A sample of 500 trees of 60 primate species.)

And the following data files:

- Artiodactyl (for multi-state analysis)

- BirdTerritory (for heterogeneous models)

- MammalBody (for continuous models)

- MammalBrainBody (for correlational or regression models)

- MammalBrainBodyGt (for multicorrelational or multivariate regression
  models)

- MammalBrainBodySampleData (an example of how to include multiple
  values per species)

- Marsupials (for continuous models, specifically variable rates)

- MammalModelB (a dataset to depict trend or model B runs - not real
  data)

- NortheastBantu (for geo models)

- Primates (for discrete models)
