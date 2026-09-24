
# likelihood.R
# Generic likelihood function constructors for impedance distributions


# ==============================================================================
# Likelihood constructor
# ==============================================================================

#' Construct likelihood functions for an impedance distribution model
#'
#' Builds log-likelihood, gradient, and Hessian functions suitable for use with
#' \code{\link[maxLik]{maxLik}} for given modelling assumptions.
#'
#' @param model An impedance model definition, as returned by
#'   \code{\link{get_impedance_model}}.
#' @param fixed Named list of fixed parameters. Any parameter not listed is
#'   treated as free and estimated.
#'
#' @return
#' A list with elements:
#' \describe{
#'   \item{logLik}{Log-likelihood function.}
#'   \item{grad}{Gradient (score) function, in vector form.}
#'   \item{hess}{Hessian function, in matrix form}
#' }
#'
#' @keywords internal
build_likelihood <- function(model, fixed) {

  # ---------------------------------------------------------------------------
  # Parameter checks
  # ---------------------------------------------------------------------------

  all_params  <- model$parameters
  free_params <- setdiff(all_params, fixed)

  if (length(free_params) == 0) {
    stop("At least one parameter must be free.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Log-likelihood
  # ---------------------------------------------------------------------------

  logLik <- function(params, data) {

    if(!all(free_params %in% names(params))) {
      stop("Missing parameters", call. = FALSE)
    }

    ll <- model$log_density(
      x      = data[[1]],
      params = params
    )

    sum(data[[2]] * ll)
  }

  # ---------------------------------------------------------------------------
  # Gradient
  # ---------------------------------------------------------------------------

  grad <- function(params, data) {

    if(!all(free_params %in% names(params))) {
      stop("Missing parameters", call. = FALSE)
    }

    grad_vector <- vapply(
      free_params,
      function(p) {
        grad_function <- model$gradient[[p]]
        g <- grad_function(
          x      = data[[1]],
          params = params
        )
        sum(data[[2]] * g)
      },
      numeric(1)
    )

    grad_vector

  }

  # ---------------------------------------------------------------------------
  # Hessian
  # ---------------------------------------------------------------------------

  hess <- function(params, data) {

    if(!all(free_params %in% names(params))) {
      stop("Missing parameters", call. = FALSE)
    }

    hess_matrix <- matrix(NA, nrow = length(free_params), ncol = length(free_params))

    for (i in seq_along(free_params)) {
      for (j in seq_along(free_params)) {
        if(j < i) {
          hess_matrix[i, j] <- hess_matrix[j, i]
          next
        }
        p_i <- free_params[i]
        p_j <- free_params[j]
        hess_function <- model$hessian[[p_i]][[p_j]]
        hess_matrix[i, j] <- sum(data[[2]] * hess_function(
          x      = data[[1]],
          params
        ))
      }
    }

    colnames(hess_matrix) <- free_params

    rownames(hess_matrix) <- free_params

    hess_matrix

  }


  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  list(
    logLik = logLik,
    grad   = grad,
    hess   = hess
  )

}
