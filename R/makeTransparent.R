
#' @title Make colour transparent.
#' @description Convert named colour (or hexadecimal code) to transparent version of the same colour.
#' @author Originally taken from [stackoverflow](https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color)
#' @param colour Either a named R colour or a hexidecimal code of format "#000000".
#' @param alpha The level of desired transparency.
#' @return Hexidecimal code for converted transparent colour.
#' @importFrom grDevices col2rgb rgb
#' @export
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
