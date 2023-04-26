#' Example Dataset: Marsupials
#'
#' A subset of the mammalian brain size dataset presented in Venditti et al 2023
#'
#' @format ## `marsupials`
#' A data frame with 173 rows and 6 columns:
#' \describe{
#'   \item{Species}{Binomial as found in the Time Tree of Life, downloaded from timetree.org Jan 2020.}
#'   \item{Order, Family}{Taxonomic information}
#'   \item{Br,Bd}{Brain and body size data in log10 grams}
#'   \item{Ref}{Source for the data.}
#'   ...
#' }
#' @source Venditti, Baker & Barton (2023) forthcoming publication
"marsupials"

#' Example Tree: Marsupials
#'
#' A subset of the timetree.org mammalian phylogeny representing 173 marsupials
#'
#' @format ## `marsupial_tree`
#' A "phylo" format Phylogenetic tree with 173 tips and 160 internal nodes:
#' \describe{
#'   \item{Tip Labels}{Binomial as found in the Time Tree of Life, downloaded from timetree.org Jan 2020.}
#'   ...
#' }
#' @source Downloaded from timetree.org Jan 2020
"marsupials"
