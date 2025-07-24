Readme
================
2025-07-15

# BayesTraitR

<!-- badges: start -->
<!-- badges: end -->

The goal of BayesTraitR is to facilitate phylogenetic comparative
analyses in BayesTraits for R users. It comprises a suite of functions
and tools that will be useful for wrangling datasets and trees into the
appropriate format for analyses within BayesTraits as well as for
interpreting and post-processing the outcomes of such analyses.

**This package is currently in development with no formal release.**

**Jump to:**

1)  [Installation](#installation)

2)  [Creating jobs to run in BayesTraits](#jobcreation)

3)  [Models implemented in BayesTraits](#models)

    *Note that in this section, we will include links to resources that
    provide a step-by-step instruction on how to implement the various
    models implemented in BayesTraits using BayesTraitR. These will be
    included as they become available. Ultimately, these help pages will
    become the package vignette upon its official release.*

4)  [Files bundled with BayesTraits](#btfiles)

<a id="installation"></a>

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

<a id="jobcreation"></a>

## Creating Analyses

The function `createBTjob()` is used to generate an input file, data
file and tree file for use in downstream BayesTraits analyses. This can
then be run interactively from within R (using `shell` or similar) or
externally using the command-line interface.

For all models that can be run in BayesTraits, we can simply specify a
list of columns found in our dataset that we wish to include using the
`cols` argument. However, for regression analyses that may include
complex contrast coding, categorical variables, interactions and
transformations (e.g. polynomials), users may wish to alternatively
specify their analysis variables using the `fm` argument.

<a id="models"></a>

**From the manual:**

BayesTraits is a computer package for performing analyses of trait
evolution among groups of species for which a phylogeny or sample of
phylogenies is available. It can be applied to the analysis of traits
that adopt a finite number of discrete states, or to the analysis of
continuously varying traits. The methods can be used to take into
account uncertainty about the model ofevolution and the underlying
phylogeny. Evolutionary hypotheses can be tested including:

- Finding rates of evolution
- Establishing correlations between traits
- Calculating ancestral state values
- Building regression models
- Predicting unknown values
- Testing for modes of evolution
  - Accelerated / decelerated rates of evolution through time
  - Magnitude of phylogenetic signal
  - Variable rate of evolution through time and within the tree
  - Fabric model, directional changes and evolvability change
  - Ornstein-Uhlenbeck processes
  - If trait change is concentrated at speciation events
  - Test for covarion evolution
- Account for uncertainty in trait values
- Interpolate unknown trait values
- Geographic model to estimate ancestral ranges and locations

Each of these are implemented in BayesTraits using a suite of defined
model types. Each of the model types is assigned a number that can be
selected when running BayesTraits to define the type of analyses that is
desired. Here, we give a brief run-down of the available model options.
For some of these, we provide worked examples of how to set up and
interpret the model using BayesTraitR.

### Multi-State models for categorical data with 2 or more states

\<<Brief description to be inserted here>\>

See [here](multistate.md) for an example of how to create and
post-process a multi-state model using BayesTraitR.

### Discrete models for two binary traits

\<<Brief description to be inserted here>\>

- Discrete: Independent

- Discrete: Dependent

- Discrete: Covarion

- Discrete: Heterogeneous

### Simple evolutionary models for continuous data

\<<Brief description to be inserted here>\>

- Random walk

- Directional

### Regression analyses for continuous data (and binary covariates)

\<<Brief description to be inserted here>\>

### Independent contrast methods

\<<Brief description to be inserted here>\>

Note: despite the name, the independent contrasts models do not output
contrasts. The fast-likelihood calculation (contrast) method
(e.g. Felsenstein 1973, Freckleton 2012) is needed for models with
complex calculations like the variable rates, variable rates regression,
and fabric models for continuous data

- Fabric model - tutorial coming soon

- Variable rates - see [here](vr.md) for an example of how to create and
  post-process a variable rates model using BayesTraitR.

- variable rates regression

### Geographic methods (Geo)

\<<Brief description to be inserted here>\>

<a id="btfiles"></a>

## BayesTraits Files

This package contains the sample files bundled with BayesTraits. This
includes the following datasets:

#### Artiodactyls

A example dataset and tree sample (n = 500) for multi-state analysis
with 17 taxa. Trait for example only; not to be used as real data.

#### Birds

An example dataset and tree for heterogeneous or discrete (binary)
analysis with 5669 taxa. Trait for example only; not to be used as real
data.

#### Mammals

Example body size, brain size, and gestation length datasets and
corresponding tree sample (n = 50) for continuous (single-trait),
correlational, and regression analyses with 40 taxa.

#### Marsupials

An example dataset and tree for continuous analysis with 246 taxa. Trait
is body size.

#### Bantu

An example coordinate (latitude/longitude) dataset and tree for 85
analysis with 85 Bantu languages.

#### Primates

An example dataset and tree sample (n = 500) for discrete analysis with
60 taxa. Trait for example only; not to be used as real data.
