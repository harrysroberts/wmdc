
# wmdc.R
# Main function for calibrating an impedance function using the WMDC method

#' Calibrate an impedance function using a marginal impedance distribution
#'
#' This function calibrates an impedance function using the marginal
#' distribution of an impedance variable, expressed as a frequency table or
#' vector of observations.
#'
#' @param data An observed marginal impedance distribution in one of the 
#' following two forms:
#' \itemize{
#'    \item A data frame containing a frequency table for the impedance
#' variable. The data frame should have two columns: the first column contains
#' the impedance values (as individual numeric values, not as ranges), and the
#' second column contains the corresponding frequencies.
#' \item A vector containing numeric impedance values
#' }
#'
#' @param shape A numeric value representing the fixed or starting value for the
#' shape parameter r in the power exponential form. A shape parameter of 1
#' results in an exponential form, while one of 2 results in a Gaussian form
#' The default value is 1, i.e. exponential.
#'
#' @param dimension A numeric value representing the fixed effective spatial
#' dimension of the study environment, or the starting value for the dimension
#' parameter when `free_dimension = TRUE`. Expected range is [1,2], including
#' non-integer values, as networks have a fractal structure spanning a subset
#' the 2-D plane. Default is 2, reflecting free movement across 2-D space.
#'
#' @param free_shape A Boolean value indicating whether to allow the shape to be
#' a free parameter in the calibration process. If TRUE, the shape parameter
#' will be estimated from the data rather than being fixed at the value
#' indicated by `shape`. Default is FALSE.
#'
#' @param free_dimension A Boolean value indicating whether to allow the spatial
#' dimension to be a free parameter in the calibration process. If TRUE, the
#' spatial dimension will be estimated from the data rather than being fixed at
#' the value indicated by `dimension`. Default is FALSE.
#'
#' @param maxLik_output A Boolean value indicating whether the full output from
#' the maxLik operation should be returned in the output list
#'
#' @return A list containing the following elements:
#' \itemize{
#'    \item \code{impedance_function}: A function that computes the calibrated
#'    impedance decay for a given impedance value, based on the specified form
#'    and calibrated parameters.
#'    \item \code{weight_function}: A weight function capturing the 
#'    distortion induced by systematic variation in the spatial structure of 
#'    origins and destinations across impedance.
#'    \item \code{density_function}: The probability density function for the
#'    marginal distribution of impedance, based on the calibrated parameters and
#'    specified form.
#'    \item \code{survival_function}: The survival function (P(X>=x)) for the
#'    marginal distribution of impedance, based on the calibrated parameters and
#'    specified form.
#'    \item \code{parameters}: A list of the calibrated parameters of the
#'    impedance function.
#'    \item \code{parameters_se}: A list of standard errors of the calibrated
#'    parameters of the impedance function.
#'    \item \code{dimension}: The effective spatial dimension used in the
#'    calibration, either fixed or estimated.
#'    \item \code{dimension_se}: The standard error of the estimated spatial
#'    dimension, if it was estimated; otherwise NA.
#'    \item \code{log_likelihood}: The log-likelihood value of the calibrated 
#'    model, indicating the goodness of fit of the model to the data.
#'    \item \code{AIC}: The Akaike Information Criterion for the fitted model
#'    indicating the relative quality of the model fit, with lower values
#'    indicating a better fit relative to other models.
#'    \item \code{survival_R_sq}: The R-squared value for the fit of the
#'    estimated survival function to the observed survival function, indicating
#'    how well the calibrated model captures the distribution of impedance in
#'    the data. This is not recommended for model selection but can be used as a
#'    descriptive measure of fit, and aids comparison with the commonly used
#'    method of fitting impedance functions to the observed survival function of
#'    impedance values.
#'    \item \code{maxLik_output}: (optional) The full output from the maxLik
#'    operation, including parameter estimates, standard errors, convergence
#'    information, and more. This is included in the output list if
#'    `maxLik_output`is set to TRUE, and is useful for users who want to inspect
#'    the details of the optimization process or perform additional analyses
#'    based on the maxLik results.
#' }
#'
#' @importFrom maxLik maxLik
#' @importFrom stats AIC coef
#' @importFrom miscTools stdEr
#'
#' @export
wmdc <- function(
    data,
    shape = 1,
    dimension = 2,
    free_shape = FALSE,
    free_dimension = FALSE,
    maxLik_output = FALSE
    ) {

  # ---------------------------------------------------------------------------
  # Data validation and wrangling
  # ---------------------------------------------------------------------------
  
  # Coerce data to data frame
  data <- tryCatch(
    as.data.frame(data),
    error = function(e)
      stop("Data must be coercible to a data frame.")
  )
  
  # Check that data has at least one column
  if (ncol(data) == 0) {
    stop("Data contains no columns.")
  }
  
  # Remove invalid values from first column
  bad_impedance <- !is.finite(data[[1]]) | data[[1]] < 0
  
  n_bad_impedance <- sum(bad_impedance)
  
  if (n_bad_impedance > 0) {
    warning(sprintf(
      "Removed %d rows due to NA, Inf, or negative values in the first column.",
      n_bad_impedance
    ))
  }
  
  data <- data[!bad_impedance, , drop = FALSE]
  
  # Check whether frequency column exists
  if (ncol(data) >= 2 && is.numeric(data[[2]])) {
    
    # Remove rows where frequencies are not positive numbers
    bad_freq <- !is.finite(data[[2]]) | data[[2]] < 0
    
    n_bad_freq <- sum(bad_freq)
    
    if (n_bad_freq > 0) {
      warning(sprintf(
        "Removed %d rows due to NA, Inf, or negative values in the frequency column.",
        n_bad_freq
      ))
    }
    
    data <- data[!bad_freq, , drop = FALSE]
    
  } else {
    
    # Extract frequency table
    counts <- table(data[[1]])
    
    # Convert data to frequency table data frame
    data <- data.frame(impedance = as.numeric(names(counts)), count = as.numeric(counts))
    
  }
  
  # ---------------------------------------------------------------------------
  # Other input validation
  # ---------------------------------------------------------------------------
  
  # Check shape is a positive numeric value
  if (!is.numeric(shape) || shape <= 0) {
    stop("`shape` must be a positive numeric value.")
  }

  # Check that dimension is a positive numeric value
  if (!is.numeric(dimension) || dimension <= 0) {
    stop("`dimension` must be a positive numeric value.")
  }

  # Check that free_shape is a logical value
  if (!is.logical(free_shape)) {
    stop("`free_shape` must be a logical value (TRUE or FALSE).")
  }

  # Check that free_dimension is a logical value
  if (!is.logical(free_dimension)) {
    stop("`free_dimension` must be a logical value (TRUE or FALSE).")
  }

  # Check that maxLik_output is a logical value
  if (!is.logical(maxLik_output)) {
    stop("`maxLik_output` must be a logical value (TRUE or FALSE).")
  }

  # ---------------------------------------------------------------------------
  # Extract required log-likelihood, gradient and Hessian functions
  # ---------------------------------------------------------------------------



  if(!free_shape & !free_dimension) {

    likelihood <- build_likelihood(ggamma_model, fixed = c("r","alpha"))

    logLik <- function(free_params,data){

      params <- list(
        alpha = dimension,
        r = shape,
        lambda = free_params[1]
        )

      likelihood$logLik(params, data)

      }

    grad <- function(free_params,data){

      params <- list(
        alpha = dimension,
        r = shape,
        lambda = free_params[1]
        )

      likelihood$grad(params, data)

      }

    hess <- function(free_params,data){

      params <- list(
        alpha = dimension,
        r = shape,
        lambda = free_params[1]
        )

        likelihood$hess(params, data)
    }

    start = c(1)

    parameters = list(lambda = NA)

    parameters_se = list(lambda = NA)

  } else if (!free_shape & free_dimension) {

    likelihood <- build_likelihood(ggamma_model, fixed = c("r"))

    logLik <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = shape,
        lambda = free_params[2]
      )

      likelihood$logLik(params, data)

    }

    grad <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = shape,
        lambda = free_params[2]
      )

      likelihood$grad(params, data)

    }

    hess <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = shape,
        lambda = free_params[2]
      )

      likelihood$hess(params, data)

    }

    start = c(dimension,1)

    parameters = list(alpha = NA, lambda = NA)

    parameters_se = list(alpha = NA, lambda = NA)

  } else if (free_shape & !free_dimension) {

    likelihood <- build_likelihood(ggamma_model, fixed = c("alpha"))

    logLik <- function(free_params,data){

      params <- list(
        alpha = dimension,
         r = free_params[1],
        lambda = free_params[2]
      )

      likelihood$logLik(params, data)

    }

    grad <- function(free_params,data){

      params <- list(
        alpha = dimension,
        r = free_params[1],
        lambda = free_params[2]
      )

      likelihood$grad(params, data)

    }

    hess <- function(free_params,data){

      params <- list(
        alpha = dimension,
        r = free_params[1],
        lambda = free_params[2]
      )

      likelihood$hess(params, data)

    }

    start = c(shape,1)

    parameters = list(r = NA, lambda = NA)

    parameters_se = list(r = NA, lambda = NA)

  } else {

    likelihood <- build_likelihood(ggamma_model, fixed = c(""))

    logLik <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = free_params[2],
        lambda = free_params[3]
        )

      likelihood$logLik(params, data)

    }

    grad <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = free_params[2],
        lambda = free_params[3]
        )

      likelihood$grad(params, data)

    }

    hess <- function(free_params,data){

      params <- list(
        alpha = free_params[1],
        r = free_params[2],
        lambda = free_params[3]
        )

      likelihood$hess(params, data)

    }

    start = c(dimension,shape,1)

    parameters = list(alpha = NA, r = NA, lambda = NA)

    parameters_se = list(alpha = NA, r = NA, lambda = NA)

  }

  # ---------------------------------------------------------------------------
  # Perform maximum likelihood estimation (with NaN warnings suppressed)
  # ---------------------------------------------------------------------------

  result <- suppressWarnings(
    maxLik(
      logLik = logLik,
      grad = grad,
      hess = hess,
      start = start,
      data = data
    )
  )

  # ---------------------------------------------------------------------------
  # Extract numerical outputs
  # ---------------------------------------------------------------------------

  #present dimension (alpha) separately from the impedance function parameters
  #note that when it is estimated, alpha is always defined to be the first parameter
  if("alpha" %in% names(parameters)){
    dimension = coef(result)[1]
    dimension_se = stdEr(result)[1]
  } else {
    dimension_se = NA
  }

  #extract the parameters estimates in list format
  for(i in seq_along(parameters)){
    parameters[[i]] = coef(result)[i]
    parameters_se[[i]] = stdEr(result)[i]
  }

  #remove alpha from the impedance parameter list
  parameters$alpha = NULL
  parameters_se$alpha = NULL

  # ---------------------------------------------------------------------------
  # Define impedance, spatial weight, density and survival functions
  # ---------------------------------------------------------------------------

  weight_body <- substitute(x^alpha_minus_1, list(alpha_minus_1 = dimension-1))

  if (!free_shape & shape == 1) {

    impedance_body <- substitute(
      exp(-lambda * x),
      list(lambda = parameters$lambda)
      )

    constant <- parameters$lambda^dimension/gamma(dimension)

  } else if (!free_shape & shape != 1) {

    impedance_body <- substitute(
      exp(-(lambda * x)^r),
      list(lambda = parameters$lambda, r = shape)
      )

    constant <- shape*parameters$lambda^dimension/gamma(dimension/shape)

  } else {

    impedance_body <- substitute(
      exp(-(lambda * x)^r),
      list(lambda = parameters$lambda,r = parameters$r)
      )

    constant <- parameters$r*parameters$lambda^dimension/gamma(dimension/parameters$r)

  }

  density_body <- Reduce(
    function(a, b) call("*", a, b),
    list(constant, weight_body, impedance_body)
    )

  if (!free_shape & shape == 1 & !free_dimension & dimension == 1){

    survival_body = impedance_body

  } else if (!free_shape & shape == 1 & !free_dimension & dimension == 2){

    survival_body = substitute(
      (1 + lambda * x) * exp(-(lambda * x)),
      list(lambda = parameters$lambda)
      )

  } else if (!free_shape & shape == 2 & !free_dimension & dimension == 1) {

    survival_body = substitute(
      2 * pnorm(- lambda * sqrt(2) * x , lower.tail = FALSE),
      list(lambda = parameters$lambda)
      )

  } else if (!free_shape & shape == 2 & !free_dimension & dimension == 2) {

    survival_body = impedance_body

  } else if (!free_shape & !(shape %in% c(1,2))) {

    survival_body = substitute(
      pgamma((lambda * x)^r,c,lower.tail = FALSE),
      list(
        lambda = parameters$lambda,
        r = shape,
        c = dimension/shape
        )
      )

  } else {

    survival_body = substitute(
      pgamma((lambda * x)^r,c,lower.tail = FALSE),
      list(
        lambda = parameters$lambda,
        r = parameters$r,
        c = dimension/parameters$r
        )
      )

  }

  impedance_function <- function(x) {}
  body(impedance_function) <- impedance_body

  weight_function <- function(x) {}
  body(weight_function) <- weight_body

  density_function <- function(x) {}
  body(density_function) <- density_body

  survival_function <- function(x) {}
  body(survival_function) <- survival_body

  # ---------------------------------------------------------------------------
  # Compute survival R^2
  # ---------------------------------------------------------------------------

  x = c(0:max(data[[1]]))

  S_estimated = survival_function(x)

  S_observed = sapply(x, function(y) sum(data[[2]][data[[1]] >= y]))/sum(data[[2]])

  survival_R_sq = 1 - sum((S_estimated - S_observed)^2)/sum((S_observed - mean(S_observed))^2)


  # ---------------------------------------------------------------------------
  # Create list output
  # ---------------------------------------------------------------------------

  output <- list(

    impedance_function = impedance_function,

    weight_function = weight_function,

    density_function = density_function,

    survival_function = survival_function,

    parameters = parameters,

    parameters_se = parameters_se,

    dimension = dimension,

    dimension_se = dimension_se,

    log_likelihood = result$maximum,

    AIC = AIC(result)[1],

    survival_R_sq = survival_R_sq

  )

  if(maxLik_output){

    output$maxLIk_output <- result

  }

  output

}







