#' Example Tree: Artiodactyl_trees
#'
#' A sample of 500 trees spanning 17 artiodactyl taxa. Used in some of the BayesTraits example analyses.
#'
#' @format ## `Artiodactyl_trees`
#' A \code{multiPhylo} format sample of 500 phylogenetic trees with 17 tips:
#' @source BayesTraits v4.1 example files.
"Artiodactyl_trees"

#' Example Tree: Bird_tree
#'
#' A phylogeny of 5669 bird species.
#'
#' @format ## `Bird_tree`
#' A \code{Phylo} format phylogenetic tree of birds with 5669 species:
#' @source BayesTraits v4.1 example files.
"Bird_tree"

#' Example Tree: Mammal_trees
#'
#' A sample of 50 trees spanning 40 mammal taxa. Used in some of the BayesTraits example analyses.
#'
#' @format ## `Mammal_trees`
#' A \code{multiPhylo} format sample of 50 phylogenetic trees with 40 tips:
#' @source BayesTraits v4.1 example files.
"Mammal_trees"

#' Example Tree: Marsupial_tree
#'
#' A phylogeny of 246 marsupial species.
#'
#' @format ## `Marsupial_tree`
#' A \code{Phylo} format phylogenetic tree of marsupials with 246 species:
#' @source BayesTraits v4.1 example files.
"Marsupial_tree"

#' Example Tree: NortheastBantu_tree
#'
#' A phylogeny of 85 Northeastern Bantu languages.
#'
#' @format ## `NortheastBantu_tree`
#' A \code{Phylo} format phylogenetic tree of 85 Northeastern bantu languages:
#' @source BayesTraits v4.1 example files.
"NortheastBantu_tree"

#' Example Tree: Primates_trees
#'
#' A sample of 500 phylogenetic trees sampling 60 primate species.
#'
#' @format ## `Primates_trees`
#' A \code{multiPhylo} format phylogenetic tree of 60 primate species:
#' @source BayesTraits v4.1 example files.
"Primates_trees"

#' Example Dataset: Artiodactyl
#'
#' A dataset used for multi-state analysis.
#'
#' @format ## `Artiodactyl`
#' A \code{data.frame} of categorical data for 17 species.
#' \describe{
#'   \item{tip_label}{Common names of 17 artiodactyl species linked to the tree.}
#'   \item{multistate}{Categorical data column.}
#' }
#' @source BayesTraits v4.1 example files.
"Artiodactyl"

#' Example Dataset: ArtiodactylMLIn
#' An input file used for multi-state analysis.
#' @format ## `ArtiodactylMLIn`
#' An input file used to run multistate analyses.
#' @source BayesTraits v4.1 example files.
"ArtiodactylMLIn"

#' Example Dataset: BirdHetCom
#' An input file used to run a heterogeneous model.
#' @format ## `BirdHetCom`
#' An input file used to run a heterogeneous model, where a different evolutionary model is fitted to the Passeriformes.
#' @source BayesTraits v4.1 example files.
"BirdHetCom"

#' Example Dataset: BirdTerritory
#' An input data file used to run a heterogeneous model.
#' @format ## `BirdTerritory`
#' A \code{data.frame} of binary data for birds used to fit a heterogeneous model, usinga lard bird phylogeny.
#' \describe{
#'   \item{tip_label}{Species names of 5669 artiodactyl species linked to the tree.}
#'   \item{Territory}{Binary data column.}
#' }
#' @source BayesTraits v4.1 example files.
"BirdTerritory"

#' Example Dataset: HarmonicMeanLh
#' An example of 10,000 log likelihoods.
#' @format ## `HarmonicMeanLh`
#' A file containing a list of log-likelihoods that can be used with the web-calculator to estimate the log harmonic mean. However, please note that this functionality has been removed from current versions of BayesTraits in favour of stepping stone sampling.
#' @source BayesTraits v4.1 example files.
"HarmonicMeanLh"

#' Example Dataset: MammalBody
#' An input file used for continuous models.
#' @format ## `MammalBody`
#' A \code{data.frame} of body size data for mammals.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{Body}{logged body size data.}
#' }
#' @source BayesTraits v4.1 example files.
"MammalBody"

#' Example Dataset: MammalBrainBody
#' An input file used for corr or reg models.
#' @format ## `MammalBrainBody`
#' A \code{data.frame} of brain and body size data for mammals.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{Brain}{logged brain size data.}
#'   \item{Body}{logged body size data.}
#' }
#' @source BayesTraits v4.1 example files.
"MammalBrainBody"

#' Example Dataset: MammalBrainBodyGt
#' An input file used for corr or reg models.
#' @format ## `MammalBrainBodyGt`
#' A \code{data.frame} of body, brain, and gestation time for mammals.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{Brain}{logged brain size data.}
#'   \item{Body}{logged body size data.}
#' }
#' @source BayesTraits v4.1 example files.
"MammalBrainBodyGt"

#' Example Dataset: MammalBrainBodySampleData
#' An example DistData file
#' @format ## `MammalBrainBodySampleData`
#' An example distData object
#' @source BayesTraits v4.1 example files.
"MammalBrainBodySampleData"

#' Example Dataset: MammalModelB
#' An example trend data file.
#' @format ## `MammalModelB`
#' A \code{data.frame} of fictional data with a trend.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{trend}{simulated trend data}
#' }
#' @source BayesTraits v4.1 example files.
"MammalModelB"

#' Example Dataset: Marsupials
#' An example continuous data file.
#' @format ## `Marsupials`
#' A \code{data.frame} of marsupial body size.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{Body}{body size data}
#' }
#' @source BayesTraits v4.1 example files.
"Marsupials"

#' Example Dataset: NortheastBantu
#' An example geo data file.
#' @format ## `NortheastBantu`
#' A \code{data.frame} of Bantu co-ordinates.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{long}{longitude}
#'   \item{lat}{latitude}
#' }
#' @source BayesTraits v4.1 example files.
"NortheastBantu"

#' Example Dataset: Primates
#' An example discrete data file.
#' @format ## `Primates`
#' A \code{data.frame} of binary data for discrete models.
#' \describe{
#'   \item{tip_label}{Species name.}
#'   \item{trait1}{first binary trait}
#'   \item{trait2}{second binary trait}
#' }
#' @source BayesTraits v4.1 example files.
"Primates"

