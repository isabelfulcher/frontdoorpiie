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
#' output <- piieffect(data=simdata,outcome="y",intermediate="m",exposure="a",covariates.outcome=c("c1","c2"),covariates.intermediate=c("c1"),covariates.exposure=c("c1","c2","c1c2"),interaction=1,astar=0)
#'
setGeneric("piieffect",
function(data,outcome,intermediate,exposure,covariates.outcome,covariates.intermediate,covariates.exposure,interaction,astar) standardGeneric("piieffect"))

#' @describeIn piieffect Generic/Function
#' @export
setMethod("piieffect", c(data = "data.frame",outcome = "character",intermediate="character",exposure="character",covariates.outcome="vector",covariates.intermediate="vector",covariates.exposure="vector",interaction="numeric",astar="numeric"),
          function(data,outcome,intermediate,exposure,covariates.outcome,covariates.intermediate,covariates.exposure,interaction,astar){

            ################
            ##   SETUP    ##
            ################

            n <- ncol(data)

            ##### functions #####
            expit <- function(x){
              output <- exp(x) / (1 + exp(x))
              return(output)
            }

            ##### model fits ####
            # formulas
            if (interaction==1){
              data$interaction <- data[,intermediate]*data[,exposure]
              outcome_model <- as.formula(paste(outcome,"~",intermediate,"+",exposure,"+","interaction","+",paste(covariates.outcome,collapse="+")))
            }else{
              outcome_model <- as.formula(paste(outcome,"~",intermediate,"+",exposure,"+",paste(covariates.outcome,collapse="+")))
            }

            intermediate_model <- as.formula(paste(intermediate,"~",exposure,"+",paste(covariates.intermediate,collapse="+")))

            exposure_model <- as.formula(paste(exposure,"~",paste(covariates.exposure,collapse="+")))

            #fits
            fit_y <- lm(outcome_model,data=data)
            fit_z <- lm(intermediate_model, data=data)
            fit_a <- glm(exposure_model,data=data,family="binomial")

            #design matrices
            matrix_y <- data.frame(model.matrix(fit_y))
            matrix_z <- data.frame(model.matrix(fit_z))
            matrix_a <- data.frame(model.matrix(fit_a))

            #### estimates ####
            sigma <- summary(fit_z)$sigma

            theta_hat <- summary(fit_y)$coefficients[,1]
            beta_hat <- summary(fit_z)$coefficients[,1]
            alpha_hat <- summary(fit_a)$coefficients[,1]


            ################
            ## ESTIMATION ##
            ################

            ##  PARAMETRIC  ##

            # 1) MLE #




            ## SEMIPARAMETRIC ##

            matrix_z_astar <- matrix_z
            matrix_z_astar[,2] <- 0

            z_mean_astar <- as.matrix(matrix_z_astar)%*%beta_hat
            z_mean <- as.matrix(matrix_z)%*%beta_hat
            a_mean <- expit(as.matrix(matrix_a)%*%alpha_hat)
            y_mean <- as.matrix(matrix_y)%*%theta_hat

            matrix_y_a_mean <- matrix_y
            if (interaction == 1){
              matrix_y_a_mean[,exposure] <- a_mean
              matrix_y_a_mean[,interaction] <- a_mean*matrix_y[,intermediate]
            }else{matrix_y_a_mean[,exposure] <- a_mean}

            matrix_y_z_mean <- matrix_y
            if(interaction == 1){
              matrix_y_z_mean[,intermediate] <- z_mean_astar
              matrix_y_z_mean[,interaction] <- z_mean_astar*matrix_y[,exposure]
            } else{matrix_y_z_mean[,intermediate] <- z_mean_astar}

            matrix_y_az_mean <- matrix_y
            if(interaction == 1){
              matrix_y_az_mean[,exposure] <- a_mean
              matrix_y_az_mean[,intermediate] <- z_mean
              matrix_y_az_mean[,interaction] <- a_mean*z_mean
            } else{matrix_y_az_mean[,exposure] <- a_mean
            matrix_y_az_mean[,intermediate] <- z_mean}

            sum_a <- as.matrix(matrix_y_a_mean)%*%theta_hat
            sum_z <- as.matrix(matrix_y_z_mean)%*%theta_hat
            sum_a_z <- as.matrix(matrix_y_az_mean)%*%theta_hat

            # 2) SP 1  #

            psi_sp1 <- data[,outcome]* (dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))

            # 3) SP 2  #

            psi_sp2 <- ((1-data[,intermediate])/(1-a_mean))*sum_a

            # 4) SP DR #

            psi_dr <- (data[,outcome] - y_mean)*
                           (dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))
                            + ((1-data[,exposure])/(1-a_mean))*(sum_a - sum_a_z)
                            + sum_z


            ################
            ##  INFERENCE ##
            ################

            return(c(mean(psi_dr),mean(psi_sp1)))
            })
