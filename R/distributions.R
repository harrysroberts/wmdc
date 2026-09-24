
# distributions.R
# Specifications of the distribution models used for fitting marginal impedance.

# ==============================================================================
# Generalised gamma impedance distribution model
# ==============================================================================

#' Generalised gamma impedance distribution model
#'
#' A flexible impedance model based on the generalised gamma distribution.
#' This model is used to represent marginal impedance distributions of travel
#' times and decomposes naturally into:
#'
#' * an impedance function: \eqn{f(x) = \exp\{-(\lambda x)^r\}}
#' * a spatial weight function: \eqn{\g(x) = x^{\alpha - 1}}
#'
#' As such it assumes that the impedance function is in power exponential form,
#' and that the spatial weight grows as a power function of impedance.
#'
#' The full density is:
#' \deqn{
#' \p(x) =
#'   r \frac{\lambda^\alpha}{\Gamma(\alpha / r)}
#'   x^{\alpha - 1} \exp\{-(\lambda x)^r\}
#' }
#'
#' @format
#' A named list defining the model form with components:
#' \describe{
#'   \item{name}{Character string identifying the model.}
#'   \item{parameters}{Character vector of parameter names.}
#'   \item{log_density}{Function computing the log-density.}
#'   \item{gradient}{Named list of score functions.}
#'   \item{hessian}{Named list of second derivatives.}
#' }
#'
#' @details
#' Parameter interpretation:
#' \itemize{
#'   \item \code{r}: functional form of impedance (shape).
#'   \item \code{alpha}: effective spatial dimension
#'   \item \code{lambda}: rate parameter controlling travel-time decay.
#' }
#'
#' @keywords internal
ggamma_model <- list(

  name = "generalised_gamma",

  parameters = c("r", "alpha", "lambda"),

  # ---------------------------------------------------------------------------
  # Log-density
  # ---------------------------------------------------------------------------

  log_density = function(x, params) {

    if(!("alpha" %in% names(params))){
      stop("Missing parameter alpha for computing log density", call. = FALSE)
    }

    if(!("r" %in% names(params))){
      stop("Missing parameter r for computing log density", call. = FALSE)
    }

    if(!("lambda" %in% names(params))){
      stop("Missing parameter lambda for computing log density", call. = FALSE)
    }

    log(params$r) +
      params$alpha * log(params$lambda) -
      lgamma(params$alpha / params$r) +
      (params$alpha - 1) * log(x) -
      (params$lambda * x)^params$r

  },

  # ---------------------------------------------------------------------------
  # Gradient (first derivative)
  # ---------------------------------------------------------------------------

  gradient = list(

    alpha = function(x, params) {

      if(!("alpha" %in% names(params))){
        stop("Missing parameter alpha for computing alpha first derivative", call. = FALSE)
      }

      if(!("r" %in% names(params))){
        stop("Missing parameter r for computing alpha first derivative", call. = FALSE)
      }

      if(!("lambda" %in% names(params))){
        stop("Missing parameter lambda for computing alpha first derivative", call. = FALSE)
      }

      log(params$lambda) -
        (1 / params$r) * digamma(params$alpha / params$r) +
        log(x)

    },

    r = function(x, params) {

      if(!("alpha" %in% names(params))){
        stop("Missing parameter alpha for computing r first derivative", call. = FALSE)
      }

      if(!("r" %in% names(params))){
        stop("Missing parameter r for computing r first derivative", call. = FALSE)
      }

      if(!("lambda" %in% names(params))){
        stop("Missing parameter lambda for computing r first derivative", call. = FALSE)
      }

      1 / params$r +
        (params$alpha / params$r^2) * digamma(params$alpha / params$r) -
        (params$lambda * x)^params$r * log(params$lambda * x)

    },

    lambda = function(x, params) {

      if(!("alpha" %in% names(params))){
        stop("Missing parameter alpha for computing lambda first derivative", call. = FALSE)
      }

      if(!("r" %in% names(params))){
        stop("Missing parameter r for computing lambda first derivative", call. = FALSE)
      }

      if(!("lambda" %in% names(params))){
        stop("Missing parameter lambda for computing lambda first derivative", call. = FALSE)
      }

      params$alpha / params$lambda -
        params$r / params$lambda * (params$lambda * x)^params$r

    }
  ),

  # ---------------------------------------------------------------------------
  # Hessian (second derivatives)
  # ---------------------------------------------------------------------------

  hessian = list(

    alpha = list(

      alpha = function(x, params) {

        if(!("alpha" %in% names(params))){
          stop("Missing parameter alpha for computing alpha-alpha second derivative", call. = FALSE)
        }

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing alpha-alpha second derivative", call. = FALSE)
        }

        -1 / params$r^2 * trigamma(params$alpha / params$r)

      },

      r = function(x, params) {

        if(!("alpha" %in% names(params))){
          stop("Missing parameter alpha for computing alpha-r second derivative", call. = FALSE)
        }

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing alpha-r second derivative", call. = FALSE)
        }

        1 / params$r^2 * digamma(params$alpha / params$r) +
          params$alpha / params$r^3 * trigamma(params$alpha / params$r)

      },

      lambda = function(x, params) {

        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing alpha-lambda second derivative", call. = FALSE)
        }

        1 / params$lambda

      }

    ),

    r = list(

      alpha = function(x, params) {

        if(!("alpha" %in% names(params))){
          stop("Missing parameter alpha for computing r-alpha second derivative", call. = FALSE)
        }

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing r-alpha second derivative", call. = FALSE)
        }

        1 / params$r^2 * digamma(params$alpha / params$r) +
          params$alpha / params$r^3 * trigamma(params$alpha / params$r)

      },

      r = function(x, params) {

        if(!("alpha" %in% names(params))){
          stop("Missing parameter alpha for computing r-r second derivative", call. = FALSE)
        }

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing r-r second derivative", call. = FALSE)
        }

        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing r-r second derivative", call. = FALSE)
        }

        -1 / params$r^2 -
          2 * params$alpha / params$r^3 * digamma(params$alpha / params$r) -
          params$alpha^2 / params$r^4 * trigamma(params$alpha / params$r) -
          (params$lambda * x)^params$r * (log(params$lambda * x))^2

      },

      lambda = function(x, params) {

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing r-lambda second derivative", call. = FALSE)
        }

        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing r-lambda second derivative", call. = FALSE)
        }

        -(params$lambda * x)^params$r / params$lambda *
          (1 + params$r * log(params$lambda * x))

      }

    ),

    lambda = list(

      alpha = function(x, params) {

        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing lambda-alpha second derivative", call. = FALSE)
        }

        1 / params$lambda

      },

      r = function(x, params) {

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing lambda-r second derivative", call. = FALSE)
        }
        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing lambda-r second derivative", call. = FALSE)
        }

        -(params$lambda * x)^params$r / params$lambda *
          (1 + params$r * log(params$lambda * x))

      },

      lambda = function(x, params) {

        if(!("alpha" %in% names(params))){
          stop("Missing parameter alpha for computing lambda-lambda second derivative", call. = FALSE)
        }

        if(!("r" %in% names(params))){
          stop("Missing parameter r for computing lambda-lambda second derivative", call. = FALSE)
        }

        if(!("lambda" %in% names(params))){
          stop("Missing parameter lambda for computing lambda-lambda second derivative", call. = FALSE)
        }

        -params$alpha / params$lambda^2 -
          params$r * (params$r - 1) / params$lambda^2 * (params$lambda * x)^params$r

      }
    )
  )
)