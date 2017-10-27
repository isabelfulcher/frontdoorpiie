#' piieffect
#'
#' Long description for function goes here
#' and here...
#'
#' @param value A character that will be returned
#' @import stats
#' @import numDeriv
#' @export
#' @author Isabel Fulcher
#' @examples
#'
#' output <- piieffect("hello")
#'
setGeneric("piieffect",
function(value) standardGeneric("piieffect"))

#' @describeIn piieffect Generic/Function
#' @export
setMethod("piieffect", c(value = "character"),
          function(value){ return(value) })
