Creating and post-processing multi-state models
================
2025-07-15

## Step 1: Generate the input files

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

## Step 2: Run the analysis

These files can then be used to run BayesTraits, either separately
through your command-line interface or using shell (or similar) options
direction from within R.

``` r
shell("BayesTraitsV4.exe MultiStateML.trees MultiStateML-001.txt < MultiStateML.infile")
```

Running the above command will create a series of files in your current
working directory:

- MultiStateML-001.txt

- MultiStateML.trees

- MultiStateML.infile

## Step 3: Post-process the analysis

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

## Step 4: Plot the transition rates

**Coming soon…**
