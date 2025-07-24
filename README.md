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

For now, tutorials are included here on the homepage but they will soon
be moved to separate tutorial files and upon package release,
incorporated into the vignette.

**Jump to:**

1)  [Installation](#installation)

2)  [Creating jobs to run in BayesTraits](#jobcreation)

3)  [Multi-State Model Example](#multistate)

4)  [Variable-Rates Model Example](vr.Rmd)

5)  [Geo Model](#geo)

6)  [Fabric model](#fabric)

7)  [Files bundled with BayesTraits](#btfiles)

<a id="installation"></a>

## 1 Installation

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

### 2 Creating Analyses

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

More detailed descriptions of how to create all model types implemented
in BayesTraits can be found in the tutorial pages (linked above).

<a id="multistate"></a>

### 3 A worked example of creating and post-processing a multi-state model

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
createBTjob(cols = "multistate", dataset = Artiodactyl, tree = Artiodactyl_trees, jobname = "MultiStateML", model = 1, MCMC = F)
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

Running BayesTraits with these files will create a log file in your
current working directory:

- MultiStateML-001.txt.Log.txt

Now, we can read in the log file into our R workspace using the function
`readBTlog` from this package.

``` r
log = readBTlog("inst/extdata/MultiStateML-001.txt.Log.txt")
```

We can view this like we would any other R data.frame. Let’s look at
what the output looks like.

``` r
head(log)
```

**Coming soon…** Plot the transition rates

<a id="vr"></a>

## 4 Creating and post-processing variable rates models

There is a specific set of commands and functions associated with the
variable rates model and the variable rates regression model. We can use
these functions to import the files output by variable rates analyses
into R and post-process into results ready for interpretation.

There are three types of variable rates models that can be run that each
have very slight differences in the output files.

1)  A variable rates model

2)  A variable rates model across a sample of trees

3)  A variable rates model across a tree sample, forcing equal tree
    sampling. Note that this option is not limited to variable rates
    models and is equivalent to running analyses sequentially across
    each tree.

### Example: A single-trait variable rates model

We can create a single-trait variable rates model over a single tree
using the example files that are bundled with the BayesTraits download.
These example files are also included as a part of this package.

Please note that for this and all worked examples, the datasets are
small and incomplete and should not be used for any scientific
interpretation. They are included as illustrative examples of how to
interact with BayesTraits and BayesTraitR only.

#### Step 1: Generate the input files

To **create a single-trait variable rates model** across the Mammal
dataset and tree sample, we can use the following commands:

``` r
createBTjob(cols = "Body", dataset = MammalBody, tree = Mammal_trees, jobname = "MammalBody_VR", model = 7, MCMC = T, optarg = c("varrates", "stones 500 10000"))
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

#### Step 2: Run the analysis

These files can then be used to run BayesTraits, either separately
through your command-line interface or using shell (or similar) options
direction from within R.

This will create a number of output files in your current working
directory:

- MammalBody-001.txt.Log.txt

- MammalBody-001.txt.Output.trees

- MammalBody-001.txt.Schedule.txt

- MammalBody-001.txt.Stones.txt

- MammalBody-001.txt.VarRates.txt

#### Step 3: Post-process the analysis

Now, we can read in the log file into our R workspace using `readBTlog`.

``` r
log = readBTlog("inst/extdata/MultiStateML-001.txt.Log.txt")
```

We can view this like we would any other R data.frame. Let’s look at
what the output looks like.

``` r
head(log)
```

Let’s summarize the variable rates output, and create a stretched tree :

``` r
res = summarizeVR(vrfile = "inst/extdata/MammalBody_VR-001.txt.VarRates.txt", tree = read.nexus("inst/extdata/Mammal_Consensus.trees"))

stretchedtree = scaleTree(res, frequency = 95, magnitude = 1, type = "mean")


plot(stretchedtree)
```

**NOTE**: because we input a sample of trees, we must specify only a
single tree to summarize on. Here we calculate and use a consensus tree
but this can be any tree of your choice - as long as it contains the
same taxa that were included in the analysis.

Any branches that were stretched that do not exist in your single
summarizing tree will be assigned to their most recent common ancestor.

<a id="geo"></a>

## 5 Geo Model

coming soon…

<a id="fabric"></a>

## 6 Fabric Model

coming soon..

<a id="btfiles"></a>

## 7 BayesTraits Files

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
