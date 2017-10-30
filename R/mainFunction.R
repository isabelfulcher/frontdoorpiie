#' piieffect
#'
#' Long description for function goes here
#' and h
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
#' #Load example dataset
#' simdata <- readRDS(system.file("rds","simdata1.rds",package="frontdoorpiie"))
#' #Create an interaction term among covariates
#' simdata$c1c2 <- simdata$c1*simdata$c2
#' #Apply the function to estimate PIIE
#' example <- piieffect(data=simdata,outcome="y",intermediate="m",exposure="a",
#' covariates.outcome=1,covariates.intermediate=c("c1"),covariates.exposure=c("c1","c2","c1c2"),
#' interaction=1,astar=0)
#'
setGeneric("piieffect",
function(data,outcome,intermediate,exposure,covariates.outcome,covariates.intermediate,covariates.exposure,interaction,astar) standardGeneric("piieffect"))

#' @describeIn piieffect Generic/Function
#' @export
setMethod("piieffect", c(data = "data.frame",outcome = "character",intermediate="character",exposure="character",covariates.outcome="vector",covariates.intermediate="vector",covariates.exposure="vector",interaction="numeric",astar="numeric"),
          function(data,outcome,intermediate,exposure,covariates.outcome=c(1),covariates.intermediate=c(1),covariates.exposure=c(1),interaction=1,astar=0){

            ################
            ##   ERRORS   ##
            ################

            if (length(covariates.outcome) == 0 | covariates.outcome == 0){covariates.outcome <- 1}
            if (length(covariates.intermediate) == 0 | covariates.intermediate == 0){covariates.intermediate <- 1}
            if (length(covariates.exposure) == 0 | covariates.exposure == 0){covariates.exposure <- 1}

            if ( sum(is.na(data[,c(outcome,intermediate,exposure)])) > 0)
              stop("Variables cannot have missing values")

            if ( interaction != 1 & interaction != 0)
              stop("interaction must be a binary value")

            if ( astar != 1 & astar != 0)
              stop("astar must be a binary value")

            if ( is.numeric(data[,outcome])=="FALSE")
              stop("Outcome variable must be numeric")

            if ( is.numeric(data[,intermediate])=="FALSE")
              stop("Intermediate variable must be numeric")

            if ( is.numeric(data[,exposure])=="FALSE")
              stop("Exposure variable must be numeric")


            ################
            ##   SETUP    ##
            ################

            n <- nrow(data)

            ##### model fits ####
            # model form for outcome
            if (interaction==1){
              data$interaction <- data[,intermediate]*data[,exposure]
              outcome_model <- as.formula(paste(outcome,"~",exposure,"+",intermediate,"+","interaction","+",paste(covariates.outcome,collapse="+")))
            }else{
              outcome_model <- as.formula(paste(outcome,"~",exposure,"+",intermediate,"+",paste(covariates.outcome,collapse="+")))
            }

            # model form for intermediate
            intermediate_model <- as.formula(paste(intermediate,"~",exposure,"+",paste(covariates.intermediate,collapse="+")))

            # model form for exposure
            exposure_model <- as.formula(paste(exposure,"~",paste(covariates.exposure,collapse="+")))

            #model fits
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

            #using new variable names to account for 0s in the covariate models
            #only necessary for MLE due to the estimator form
            theta_hat_m <- theta_hat
            beta_hat_m <- beta_hat

            #estimate of E(A)
            mean_exposure <- mean(data[,exposure])

            #estimate of E(C)
            if (length(covariates.outcome)>1){
              mean_covariates_outcome <- apply(data[,covariates.outcome],2,mean)
            } else if(covariates.outcome==1) {
              mean_covariates_outcome <- 0
              theta_hat_m <- c(theta_hat,0)
            } else {mean_covariates_outcome <- mean(data[,covariates.outcome])}

            #estimate of E(AC)
            if (length(covariates.intermediate)>1){
              mean_covariates_exposure <- apply(data[,exposure]*data[,covariates.intermediate],2,mean)
            } else if (covariates.intermediate==1){
              mean_covariates_exposure <- 0
              beta_hat_m <- c(beta_hat,0)
            } else {mean_covariates_exposure <- mean(data[,exposure]*data[,covariates.intermediate]) }

            #put everything together to estimate psi
            if (interaction == 1){psi_mle <- (theta_hat_m[1] + theta_hat_m[3]*beta_hat_m[1] + theta_hat_m[3]*beta_hat_m[2]*astar
                                          + (theta_hat_m[2] + theta_hat_m[4]*beta_hat_m[1] + theta_hat_m[4]*beta_hat_m[2]*astar)*mean_exposure
                                          + (theta_hat_m[3]*beta_hat_m[3:length(beta_hat_m)] + theta_hat_m[5:length(theta_hat_m)])%*%t(t(mean_covariates_outcome))
                                          + theta_hat_m[4]*beta_hat_m[3:length(beta_hat_m)]%*%t(t(mean_covariates_exposure)))
            } else {psi_mle <- (theta_hat_m[1] + theta_hat_m[3]*beta_hat_m[1] + theta_hat_m[3]*beta_hat_m[2]*astar
                            + theta_hat_m[2]*mean_exposure
                            + (theta_hat_m[3]*beta_hat_m[3:length(beta_hat_m)] + theta_hat_m[4:length(theta_hat_m)])%*%t(t(mean_covariates_outcome)))}


            ## SEMIPARAMETRIC ##

            # design matrix for Z setting the exposure to astar for everyone
            matrix_z_astar <- matrix_z
            matrix_z_astar[,exposure] <- astar

            # predicted values for each individual
            z_mean_astar <- as.matrix(matrix_z_astar)%*%beta_hat
            z_mean <- as.matrix(matrix_z)%*%beta_hat
            a_mean <- expit(as.matrix(matrix_a)%*%alpha_hat)
            y_mean <- as.matrix(matrix_y)%*%theta_hat

            # design matrix for Y setting the exposure A to the predicted value
            matrix_y_a_mean <- matrix_y
            if (interaction == 1){
              matrix_y_a_mean[,exposure] <- a_mean
              matrix_y_a_mean[,interaction] <- a_mean*matrix_y[,intermediate]
            }else{matrix_y_a_mean[,exposure] <- a_mean}

            # design matrix for Y setting the intermediate Z to the predicted value under regime astar
            matrix_y_z_mean <- matrix_y
            if(interaction == 1){
              matrix_y_z_mean[,intermediate] <- z_mean_astar
              matrix_y_z_mean[,interaction] <- z_mean_astar*matrix_y[,exposure]
            } else{matrix_y_z_mean[,intermediate] <- z_mean_astar}

            # design matrix for Y setting the intermediate A and Z to their predicted values
            matrix_y_az_mean <- matrix_y
            if(interaction == 1){
              matrix_y_az_mean[,exposure] <- a_mean
              matrix_y_az_mean[,intermediate] <- z_mean
              matrix_y_az_mean[,interaction] <- a_mean*z_mean
            } else{matrix_y_az_mean[,exposure] <- a_mean
            matrix_y_az_mean[,intermediate] <- z_mean}

            # three sum terms from DR estimator equation
            sum_a <- as.matrix(matrix_y_a_mean)%*%theta_hat
            sum_z <- as.matrix(matrix_y_z_mean)%*%theta_hat
            sum_a_z <- as.matrix(matrix_y_az_mean)%*%theta_hat

            # 2) SP 1  #

            psi_sp1_i <- data[,outcome]* (dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))

            piie_sp1_i <- data[,outcome] - psi_sp1_i

            # 3) SP 2  #

            psi_sp2_i <- ((1-data[,intermediate])/(1-a_mean))*sum_a

            piie_sp2_i <- data[,outcome] - psi_sp2_i

            # 4) SP DR #

            psi_dr_i <- (data[,outcome] - y_mean)*
                           (dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))
                            + ((1-data[,exposure])/(1-a_mean))*(sum_a - sum_a_z)
                            + sum_z

            piie_dr_i <- data[,outcome] - psi_dr_i

            # Population level estimates for psi #
            psi_sp1 <- mean(psi_sp1_i)
            psi_sp2 <- mean(psi_sp2_i)
            psi_dr  <- mean(psi_dr_i)

            # Population level estimates for PIIE #
            piie_mle <- mean(data[,outcome]) - psi_mle
            piie_sp1 <- mean(data[,outcome]) - psi_sp1
            piie_sp2 <- mean(data[,outcome]) - psi_sp2
            piie_dr  <- mean(data[,outcome]) - psi_dr

            ################
            ##  INFERENCE ##
            ################

            # 1) MLE   #

            #closed form expression
            if(interaction ==1){piie_var_mle <- ( (var(data[,exposure])/n)*(beta_hat[2]^2)*(theta_hat[3] + theta_hat[4])^2
                                              + (mean_exposure^2)*((theta_hat[3] + theta_hat[4])^2)*vcov(fit_z)[2,2]
                                              + mean_exposure*beta_hat[2]*(mean_exposure*beta_hat[2]*vcov(fit_y)[3,3] + mean_exposure*beta_hat[2]*vcov(fit_y)[3,4])
                                              + mean_exposure*beta_hat[2]*(mean_exposure*beta_hat[2]*vcov(fit_y)[3,4] + mean_exposure*beta_hat[2]*vcov(fit_y)[4,4]) )
            } else { piie_var_mle <- ( ((mean_exposure*theta_hat[3])^2)*vcov(fit_z)[2,2]
                                   + ((mean_exposure*beta_hat[2])^2)*vcov(fit_y)[3,3]
                                   + ((beta_hat[2]*theta_hat[3])^2)*var(data[,exposure])/n )  }


            # 2) SP 1  #

            # function will take derivative of the score
            deriv_sp1 <- numDeriv::jacobian(U.sp1,c(alpha_hat,beta_hat,theta_hat,piie_sp1),data=data,outcome=outcome,intermediate=intermediate,exposure=exposure,interaction=interaction,astar=astar,matrix_y=matrix_y,matrix_z=matrix_z,matrix_a=matrix_a,sigma=sigma)

            # score function
            score_sp1 <- as.matrix(cbind( matrix_a*c(data[,exposure] - a_mean),
                                matrix_z*c(data[,intermediate] - z_mean),
                                matrix_y*c(data[,outcome] - y_mean),
                                (piie_sp1_i - piie_sp1)))

            # calculate variance matrix
            var_sp1 <- (solve(deriv_sp1)%*%t(score_sp1)%*%score_sp1%*%t(solve(deriv_sp1)))

            # variance for PIIE
            piie_var_sp1 <- var_sp1[length(c(alpha_hat,beta_hat,theta_hat,piie_sp1)),length(c(alpha_hat,beta_hat,theta_hat,piie_sp1))]

            # 3) SP 2  #

            # function will take derivative of the score
            deriv_sp2 <- numDeriv::jacobian(U.sp2,c(alpha_hat,beta_hat,theta_hat,piie_sp2),data=data,outcome=outcome,intermediate=intermediate,exposure=exposure,interaction=interaction,astar=astar,matrix_y=matrix_y,matrix_z=matrix_z,matrix_a=matrix_a,sigma=sigma)

            # score function
            score_sp2 <- as.matrix(cbind( matrix_a*c(data[,exposure] - a_mean),
                                          matrix_z*c(data[,intermediate] - z_mean),
                                          matrix_y*c(data[,outcome] - y_mean),
                                          (piie_sp2_i - piie_sp2)))

            # calculate variance matrix
            var_sp2 <- (solve(deriv_sp2)%*%t(score_sp2)%*%score_sp2%*%t(solve(deriv_sp2)))

            # variance for PIIE
            piie_var_sp2 <- var_sp2[length(c(alpha_hat,beta_hat,theta_hat,piie_sp2)),length(c(alpha_hat,beta_hat,theta_hat,piie_sp2))]


            # 4) SP DR #

            # function will take derivative of the score
            deriv_dr <- numDeriv::jacobian(U.dr,c(alpha_hat,beta_hat,theta_hat,piie_sp1),data=data,outcome=outcome,intermediate=intermediate,exposure=exposure,interaction=interaction,astar=astar,matrix_y=matrix_y,matrix_z=matrix_z,matrix_a=matrix_a,sigma=sigma)

            # score function
            score_dr <- as.matrix(cbind( matrix_a*c(data[,exposure] - a_mean),
                                          matrix_z*c(data[,intermediate] - z_mean),
                                          matrix_y*c(data[,outcome] - y_mean),
                                          (piie_dr_i - piie_dr)))

            # calculate variance matrix
            var_dr <- (solve(deriv_sp1)%*%t(score_dr)%*%score_sp1%*%t(solve(deriv_dr)))

            # variance for PIIE
            piie_var_dr <- var_dr[length(c(alpha_hat,beta_hat,theta_hat,piie_dr)),length(c(alpha_hat,beta_hat,theta_hat,piie_dr))]


            ################
            ##   OUTPUT   ##
            ################

            #95% Wald-type confidence intervals
            ci_piie_mle <- c(piie_mle-qnorm(.975)*sqrt(piie_var_mle),piie_mle+qnorm(.975)*sqrt(piie_var_mle))
            ci_piie_sp1 <- c(piie_sp1-qnorm(.975)*sqrt(piie_var_sp1),piie_sp1+qnorm(.975)*sqrt(piie_var_sp1))
            ci_piie_sp2 <- c(piie_sp2-qnorm(.975)*sqrt(piie_var_sp2),piie_sp2+qnorm(.975)*sqrt(piie_var_sp2))
            ci_piie_dr  <- c(piie_dr-qnorm(.975)*sqrt(piie_var_dr),piie_dr+qnorm(.975)*sqrt(piie_var_dr))

            output <- matrix( c(c(psi_mle,piie_mle,sqrt(piie_var_mle),ci_piie_mle),
                              c(psi_sp1,piie_sp1,sqrt(piie_var_sp1),ci_piie_sp1),
                              c(psi_sp2,piie_sp2,sqrt(piie_var_sp2),ci_piie_sp2),
                              c(psi_dr,piie_dr,sqrt(piie_var_dr),ci_piie_dr)),4,4,byrow=TRUE)

            colnames(output) <- c("Psi","PIIE","Standard Error","95% CI")
            rownames(output) <- c("MLE","SP1","SP2","DR")

            return(output)
            })


# expit function
expit <- function(x){
  output <- exp(x) / (1 + exp(x))
  return(output)
}

# Score functions for semiparametric estimators
U.sp1 <- function(estimates,data,outcome,intermediate,exposure,interaction,astar,matrix_y,matrix_z,matrix_a,sigma){

  n <- nrow(data)

  len_y <- ncol(matrix_y)
  len_z <- ncol(matrix_z)
  len_a <- ncol(matrix_a)

  alpha_hat <- estimates[1:len_a]
  beta_hat <- estimates[(len_a+1):(len_a+len_z)]
  theta_hat <- estimates[(len_a+len_z+1):(length(estimates) - 1)]
  piie_est <- estimates[length(estimates)]

  matrix_z_astar <- matrix_z
  matrix_z_astar[,exposure] <- astar

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

  psi_sp1_i <- data[,outcome]*(dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))

  piie_sp1_i<- data[,outcome] - psi_sp1_i

  piie_sp1 <- mean(data[,outcome]) - mean(psi_sp1_i)

  score_sp1 <- cbind(matrix_a*c(data[,exposure] - a_mean),
                     matrix_z*c(data[,intermediate] - z_mean),
                     matrix_y*c(data[,outcome] - y_mean),
                     (piie_sp1_i - piie_est))

  deriv <- matrix(1,1,n)%*%as.matrix(score_sp1)

  return(deriv)

}

U.sp2 <- function(estimates,data,outcome,intermediate,exposure,interaction,astar,matrix_y,matrix_z,matrix_a,sigma){

  n <- nrow(data)

  len_y <- ncol(matrix_y)
  len_z <- ncol(matrix_z)
  len_a <- ncol(matrix_a)

  alpha_hat <- estimates[1:len_a]
  beta_hat <- estimates[(len_a+1):(len_a+len_z)]
  theta_hat <- estimates[(len_a+len_z+1):(length(estimates) - 1)]
  piie_est <- estimates[length(estimates)]

  matrix_z_astar <- matrix_z
  matrix_z_astar[,exposure] <- astar

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

  psi_sp2_i <- ((1-data[,intermediate])/(1-a_mean))*sum_a

  piie_sp2_i<- data[,outcome] - psi_sp2_i

  piie_sp2 <- mean(data[,outcome]) - mean(psi_sp2_i)

  score_sp2 <- cbind(matrix_a*c(data[,exposure] - a_mean),
                     matrix_z*c(data[,intermediate] - z_mean),
                     matrix_y*c(data[,outcome] - y_mean),
                     (piie_sp2_i - piie_est))

  deriv <- matrix(1,1,n)%*%as.matrix(score_sp2)

  return(deriv)

}

U.dr <- function(estimates,data,outcome,intermediate,exposure,interaction,astar,matrix_y,matrix_z,matrix_a,sigma){

  n <- nrow(data)

  len_y <- ncol(matrix_y)
  len_z <- ncol(matrix_z)
  len_a <- ncol(matrix_a)

  alpha_hat <- estimates[1:len_a]
  beta_hat <- estimates[(len_a+1):(len_a+len_z)]
  theta_hat <- estimates[(len_a+len_z+1):(length(estimates) - 1)]
  piie_est <- estimates[length(estimates)]

  matrix_z_astar <- matrix_z
  matrix_z_astar[,exposure] <- astar

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

  psi_dr_i <- (data[,outcome] - y_mean)*
               (dnorm(data[,intermediate],z_mean_astar,sigma)/dnorm(data[,intermediate],z_mean,sigma))
                + ((1-data[,exposure])/(1-a_mean))*(sum_a - sum_a_z)
                + sum_z

  piie_dr_i<- data[,outcome] - psi_dr_i

  piie_dr <- mean(data[,outcome]) - mean(psi_dr_i)

  score_dr <- cbind(matrix_a*c(data[,exposure] - a_mean),
                     matrix_z*c(data[,intermediate] - z_mean),
                     matrix_y*c(data[,outcome] - y_mean),
                     (piie_dr_i - piie_est))

  deriv <- matrix(1,1,n)%*%as.matrix(score_dr)

  return(deriv)

}


