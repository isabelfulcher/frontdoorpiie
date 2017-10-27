#' piieffect
#'
#' Long description for function goes here
#' and here...
#'
#' @param data A dataframe
#' @param outcome The variable name for the outcome variable in data
#' @param intermediate The variable name for the intermediate variable in data
#' @param exposure The variable name for the exposure variable in data
#' @param covariates.outcome A vector of variable names for covariates to be included in the outcome model
#' @param covariates.intermediate A vector of variable names for covariates to be included in the outcome model
#' @param covariates.exposure A vector of variable names for covariates to be included in the exposure model
#' @param interaction A binary variable indicating if an interaction term between intermediate and exposure is needed
#' @param astar A numeric value for the level of the exposure wanted for comparison
#' @import stats
#' @import numDeriv
#' @export
#' @author Isabel Fulcher
#' @examples
#'
#' simdata <- readRDS(system.file("rds","simdata1.rds",package="frontdoorpiie"))
#' simdata$c1c2 <- simdata$c1*simdata$c2
#' output <- piieffect(data=simdata,outcome="y",intermediate="m",exposure="a",
#' covariates.outcome=c("c1","c2"),covariates.intermediate=c("c1"),covariates.exposure=c("c1","c2","c1c2"),interaction=1,astar=0)
#'
setGeneric("piieffect",
function(data,outcome,intermediate,exposure,covariates.outcome,covariates.intermediate,covariates.exposure,interaction,astar) standardGeneric("piieffect"))

#' @describeIn piieffect Generic/Function
#' @export
setMethod("piieffect", c(data = "data.frame",outcome = "character",intermediate="character",exposure="character",covariates.outcome="vector",covariates.intermediate="vector",covariates.exposure="vector",interaction="numeric",astar="numeric"),
          function(data,outcome,intermediate,exposure,covariates.outcome,covariates.intermediate,covariates.exposure,interaction,astar){
            return(mean(data[,outcome]))
            })
