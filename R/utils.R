#' @keywords internal
last <- function(x) { tail(x, n = 1) }

#' @keywords internal
h <- function(w) if( any( grepl( "Recycling array of length 1 in vector-array arithmetic is deprecated", w) ) ) invokeRestart( "muffleWarning" )
