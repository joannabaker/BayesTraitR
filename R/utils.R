#' @title Make colour transparent.
#' @description Convert named colour (or hexadecimal code) to transparent version of the same colour.
#' @author Originally taken from [stackoverflow](https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color)
#' @param colour Either a named R colour or a hexidecimal code of format "#000000".
#' @param alpha The level of desired transparency.
#' @return Hexidecimal code for converted transparent colour.
#' @importFrom grDevices col2rgb rgb
#' @export
#' @keywords internal
#' @examples
#' makeTransparent("red", 0.5)
#' makeTransparent("#000000", 0.1)
makeTransparent = function(colour, alpha=0.5) {

  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")

  alpha = floor(255*alpha)
  newColor = col2rgb(col=colour, alpha=FALSE)

  .makeTransparent = function(col, alpha) {rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)}

  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)

}



#' @title Pairwise comparison between two sets of values.
#' @author Jo Baker
#' @description Calculates the difference between 2 columns of data and the proportion of the new values that are greater than the original reference set. Note that these values have to be pairwise otherwise this doesn't make much sense! This can be useful for comparing differences amongst parameter values in MCMC regressions.
#' @param x A numeric vector with the reference values.
#' @param y A numeric vector with the comparison values.
#' @param pmcmc If true, instead of returning the proportion of y that is greater than x it will return the proportion of y that is DIFFERENT to x i.e. a pmcmc value (proportion of differences overlapping zero).
#' @return The proportion of the new values that are greater than (or different to, if pmcmc = t) the reference values.
#' @export
#
overlap = function(x, y, pmcmc = F){

  change = y-x
  p = length(change[change > 0]) / length(change)
  if(pmcmc) p = ifelse(p > 0.5, 1-p, p)
  return(p)}


#' @title Pairwise comparison among columns.
#' @description For every pair of values in a matrix, calculates the overlap.
#' Overlap is the proportional difference between each 2 columns of data. This can be used to identify pairwise difference matrices amongst group-level parameter estimates.
#' Note that these values have to be paired otherwise this doesn't make much sense!
#' @param mat A matrix of dimension n>1,n>2 with the values of interest. This MUST have column names.
#' @param pmcmc If true, instead of returning the proportion of y that is greater than x it will return the proportion of y that is DIFFERENT to x i.e. a pmcmc value (proportion of differences overlapping zero).
#' @return A matrix identifying the proportion of the values that overlap between each pair (pmcmc).
#' @export
#'
matrixoverlap = function(mat, pmcmc = F){

  res = matrix(data = NA, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(names(mat), names(mat)))
  for(i in 1:nrow(res)){
    for(j in 1: ncol(res)){
      if(i == j) next
      # mat[colnames(x)[i], colnames(x[j])]
      res[i,j] = 	overlap(x = mat[,names(mat)[i]], y = mat[,colnames(mat[j])], pmcmc)
    }
  }
  return(res)
}

#' @title Plot rotated histograms.
#' @description Allows for plotting of histograms in any orientation: up, down, left, or right.
#' @param A An R object obtained from the hist(plot = F) function.
#' @param direction The direction for the histogram. Can be any of "up", "down", "left", or "right".
#' @return Plots a rotated histogram.
#' @importFrom graphics arrows axis hist polygon rect text
#' @export
# functions for plotting histograms
histRotate = function(A, direction = "up"){

  if(!direction %in% c("up", "down", "left", "right")) stop("Direction must be one of \"up\", \"down\", \"left\", or \"right\"")


  if(direction == "up"){
    plot(NULL, type = "n", ylim = c(0,max(A$counts)), xlim = c(range(A$breaks)), axes = F, xlab = NA, ylab = NA)
    rect(A$breaks[1:(length(A$breaks) - 1)], 0, A$breaks[2:length(A$breaks)], A$counts)
    axis(side = 1)
  }


  if(direction == "down"){
    plot(NULL, type = "n", ylim = c(max(A$counts), 0), xlim = c(range(A$breaks)), axes = F, xlab = NA, ylab = NA)
    rect(A$breaks[1:(length(A$breaks) - 1)], 0, A$breaks[2:length(A$breaks)], A$counts)
    axis(side = 3)
  }

  if(direction == "left"){
    plot(NULL, type = "n", xlim = c(max(A$counts), 0), ylim = c(range(A$breaks)), axes = F, xlab = NA, ylab = NA)
    rect(0, A$breaks[1:(length(A$breaks) - 1)], A$counts, A$breaks[2:length(A$breaks)])
    axis(side = 4, las = 1)
  }

  if(direction == "right"){
    plot(NULL, type = "n", xlim = c(0, max(A$counts)), ylim = c(range(A$breaks)), axes = F, xlab = NA, ylab = NA)
    rect(0, A$breaks[1:(length(A$breaks) - 1)], A$counts, A$breaks[2:length(A$breaks)])
    axis(side = 2, las = 1)
  }
}



#' @title Create teletubby plots for data distributions.
#' @description This function will plot a bunch of DENSITY distributions given a series of ROWS of data.
#'     For each row in the provided matrix / data frame, this function will calculate a density distribution.
#'     Optionally, it will generate polygons to add to a new or pre-existing plot.
#' @param x A matrix or data.frame of dimension i>1,j>2 with the values of interest.
#' @param col A named R colour or hexadecimal value that defines the colour of the polygon to be plotted (or saved).
#' @param alpha The transparency value of the density distribution - for many columns, it might be useful to set this to a low value. Defaults to 0.1
#' @param log Boolean operator defining whether or not to log the data provided before calculating the density distribution. Defaults to true.
#' @param plot Boolean operator defining whether or not to generate a new plot for the current distributions. Defaults to true. If this is set to false, the function can be used to overlay new polygons onto a pre-existing plot.
#' @return A list of densities (one for each row in the original provided data.frame or matrix).
#' @export
#
# This function will plot a bunch of DENSITY distributions given a series of ROWS of data.
teletubby = function(x, col = "indianred", alpha = 0.1, log = T, plot = T){

  if(length(col) == 1) col = rep(col, nrow(x))

  # Generate list of densities (we will generate plot margins from this)
  .dlist = list()
  for(i in 1:nrow(x))
    if(log == T) .dlist[[i]] <- density(log10(x[i,])) else .dlist[[i]] = density(x[i,])

  # Calculate axis limits + create starting plot if we are not adding teletubbies to an existing plot
  if(plot){
    xlim = c(min(unlist(lapply(1:length(.dlist), function(x)min(.dlist[[x]]$x)))), max(unlist(lapply(1:length(.dlist), function(x)max(.dlist[[x]]$x)))))
    ylim = c(min(unlist(lapply(1:length(.dlist), function(x)min(.dlist[[x]]$y)))), max(unlist(lapply(1:length(.dlist), function(x)max(.dlist[[x]]$y)))))

    plot(xlim, xlim = xlim, ylim = ylim, las = 1, bty = "l", ylab = "Density", xlab = "Rate")
  }

  # Add the teletubbies!
  for(i in 1:length(.dlist))
    polygon(.dlist[[i]], col = makeTransparent(col[i], alpha = alpha), border = NA)

  return(.dlist)


}
