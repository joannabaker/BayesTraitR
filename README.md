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

1)  [Installation](#installation)

2)  [Creating jobs to run in BayesTraits](#jobcreation)

3)  [Multi-State Model Example](#multistate)

4)  [Variable-Rates Model Example](#vr)

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

### A worked example of creating and post-processing a multi-state model

We can create a multi-state model using the example files that are
bundled with the BayesTraits download. These example files are also
included as a part of this package.

Please note that for this and all worked examples, the datasets are
small and incomplete and should not be used for any scientific
interpretation. They are included as illustrative examples of how to
interact with BayesTraits and BayesTraitR only.

To **create a maximum-likelihood multi-state model** across the
Artiodactyl dataset and tree sample, we can use the following commands:

``` r
createBTjob(cols = "multistate", model = 1, dataset = Artiodactyl, tree = Artiodactyl_trees, jobname = "MultiStateML", MCMC = F)
```

Note that we specify the following arguments:

- **jobname = “MultiStateML”**: This sets an identifier string with
  which all output files will be tagged.

- **cols = “multistate”**: This specifies the column we wish to run the
  model over in the input dataset.

- **model = 1**: This specifies the multi-state model.

- **MCMC = F**: By selecting to not use MCMC, the model will revert to
  maximum-likelihood.

All other parameters are set to the default values.

Running the above command will create a series of files in your current
working directory:

- MultiStateML-001.txt

- MultiStateML.trees

- MultiStateML.infile

These files can then be used to run BayesTraits, either separately
through your command-line interface or using shell (or similar) options
direction from within R.

For example:

``` r
shell("BayesTraitsV4.exe MultiStateML.trees MultiStateML-001.txt < MultiStateML.infile")
```

Running BayesTraits with these files will create a log file in your
current working directory:

- MultiStateML-001.txt.Log.txt

Now, we can read in the log file into our R workspace using the from
this package.

``` r
log = readBTlog("inst/extdata/MultiStateML-001.txt.Log.txt")
```

We can view this like we would any other R data.frame. Let’s look at
what the output looks like.

``` r
head(log)
#>   Tree.No        Lh      qDG      qGD Root.P.D. Root.P.G.  X
#> 1       1 -7.550301 5.441884 6.066461  0.144555  0.855445 NA
#> 2       2 -7.439013 5.341572 4.576922  0.222662  0.777338 NA
#> 3       3 -7.637118 6.352144 6.476818  0.195341  0.804659 NA
#> 4       4 -7.604123 6.706603 5.772352  0.226221  0.773779 NA
#> 5       5 -7.525569 5.816163 5.410476  0.200482  0.799518 NA
#> 6       6 -7.807101 8.301755 7.480684  0.257565  0.742435 NA
```

**Coming soon…** Plot the transition rates

## A worked example of creating and post-processing a variable rates model

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

### Let’s start with (1).

We can create a single-trait variable rates model over a single tree
using the example files that are bundled with the BayesTraits download.
These example files are also included as a part of this package.

Please note that for this and all worked examples, the datasets are
small and incomplete and should not be used for any scientific
interpretation. They are included as illustrative examples of how to
interact with BayesTraits and BayesTraitR only.

To **create a single-trait variable rates model** across the Mammal
dataset and tree sample, we can use the following commands:

``` r
createBTjob(cols = "Body", model = 7, dataset = MammalBody, tree = Mammal_trees, jobname = "MammalBody_VR", MCMC = T, optarg = c("varrates", "stones 500 10000"))
```

Note that we specify the following arguments, along with the dataset and
tree objects:

- **jobname = “MammalBody_VR”**: This sets an identifier string with
  which all output files will be tagged.

- **cols = “Body”**: This specifies the column we wish to run the model
  over in the input dataset.

- **model = 7**: This specifies the independent-contrast model for a
  single trait. Note that this does not output contrasts but instead
  uses the fast likelihood calculation in order to estimate complex
  model parameterization such as the variable rates models.

- **MCMC = T**: The variable rates model is only available for MCMC
  analyses.

We also specify the following optional arguments using the **optarg =
c(…)** argument:

- **varrates** Specifies the model should estimate variable rates.

- **stones 500 10000** Specifies that the model should be run with
  stepping-stone sampling. The model will run 500 stones for 10,000
  iterations per stone.

All other parameters are set to the default values.

This will create a series of files in your current working directory:

- MammalBody_VR-001.txt

- MammalBody_VR.trees

- MammalBody_VR.infile

These files can then be used to run BayesTraits, either separately
through your command-line interface or using shell (or similar) options
direction from within R.

``` r
shell("BayesTraitsV4.exe MammalBody.trees MammalBody-001.txt < MammalBody.infile")
```

This will create a number of output files in your current working
directory:

- MammalBody-001.txt.Log.txt

- MammalBody-001.txt.Output.trees

- MammalBody-001.txt.Schedule.txt

- MammalBody-001.txt.Stones.txt

- MammalBody-001.txt.VarRates.txt

Now, we can read in the log file into our R workspace using .

``` r
log = readBTlog("inst/extdata/MultiStateML-001.txt.Log.txt")
```

We can view this like we would any other R data.frame. Let’s look at
what the output looks like.

``` r
head(log)
#>   Tree.No        Lh      qDG      qGD Root.P.D. Root.P.G.  X
#> 1       1 -7.550301 5.441884 6.066461  0.144555  0.855445 NA
#> 2       2 -7.439013 5.341572 4.576922  0.222662  0.777338 NA
#> 3       3 -7.637118 6.352144 6.476818  0.195341  0.804659 NA
#> 4       4 -7.604123 6.706603 5.772352  0.226221  0.773779 NA
#> 5       5 -7.525569 5.816163 5.410476  0.200482  0.799518 NA
#> 6       6 -7.807101 8.301755 7.480684  0.257565  0.742435 NA
```

Let’s summarize the variable rates output, and create a stretched tree :

``` r
res = summarizeVR(vrfile = "inst/extdata/MammalBody_VR-001.txt.VarRates.txt", treefile = "inst/extdata/Mammal_Consensus.trees")
#> Extracting VR information from file...
#> Extracting VR information from file...
#>  Linking to tree...
#>      Summarizing Rates...

stretchedtree = read.nexus("inst/extdata/Mammal_Consensus.trees")
stretchedtree$edge.length = stretchedtree$edge.length * res$ratesummary$medianscalar[-1]

plot(stretchedtree)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

**NOTE**: because we input a sample of trees, we must specify only a
single tree to summarize on. Here we calculate and use a consensus tree
but this can be any tree of your choice - as long as it contains the
same taxa that were included in the analysis.

Any branches that were stretched that do not exist in your single
summarizing tree will be assigned to their most recent common ancestor.

## BayesTraits Files

This package contains the sample files bundled with BayesTraits. This
includes the following datasets:

\*[Artiodactyls](#artiodactyls): A example dataset and tree sample (n =
500) for multi-state analysis with 17 taxa. I don’t know what the data
is and I think I will replace this dataset with activity pattern at some
point in the near future.

\*[Birds](#birds): An example dataset and tree for heterogeneous or
discrete (binary) analysis with 5669 taxa. I don’t know what the data
is; this will be updated in due course.

\*\[Mammals\]: Example body size, brain size, and gestation length
datasets and corresponding tree sample (n = 50) for continuous
(single-trait), correlational, and regression analyses with 40 taxa.
Includes

\*\[Marsupials\]: An example body size dataset and tree for continuous
analysis with 246 taxa.

\*\[Bantu\]: An example coordinate (latitude/longitude) dataset and tree
for 85 analysis with 85 Bantu languages.

\*\[Primates\]: An example dataset and tree sample (n = 500) for
discrete analysis with 60 taxa. I don’t know what the data is; this will
be updated in due course.

### Artiodactyls

This is what the consensus artiodactyl tree looks like:
<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

And here is the data

``` r
class(Artiodactyl)
#> [1] "data.frame"
head(Artiodactyl)
#>    tip_label multistate
#> 1 Chevrotain          D
#> 2    Giraffe          D
#> 3       Goat          D
#> 4      Sheep          D
#> 5  Pronghorn          G
#> 6    Buffalo          D
```

### Birds

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" /> \###
Mammals

This is what the consensus mammal tree looks like:
<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

And here is the data

``` r
class(MammalBody)
#> [1] "data.frame"
head(MammalBody)
#>          tip_label     Body
#> 1          Opossum 3.454326
#> 2   Diprotodontian 3.680879
#> 3         Elephant 6.564903
#> 4            Hyrax 3.390582
#> 5         Tenrecid 2.065953
#> 6 Lo_Ear_Ele_shrew 1.667453
```

In addition to the data files, we also have the following objects:

- MammalBrainBodySampleData (an example of how to include multiple
  values per species)

- MammalModelB (a dataset to depict trend or model B runs - not real
  data)
