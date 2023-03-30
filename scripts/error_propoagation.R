#' Propagate error by formula using formula: sigma = sqrt( (df/dx)^2 sx^2 + (df/dy)^2 sy^2 + ... ). 
#' This only works for independent variables as there are no covariance terms!
#' All errors must be supplied as named \code{...} parameters where the names all start with \code{d} and match the variable names in the \code{formula}.
#' @param eq the equation for calculating a new quantity that should be used to propagate error
#' @param ... the error terms for the variables in the formula. Must all start with \code{d} and match all variables in the code. Omit a term by setting its error to 0/
#' @examples
#' tibble::tibble(
#'   x  = runif(5), 
#'   y  = runif(5), 
#'   x_err = runif(5)/10, 
#'   y_err = runif(5)/10,
#'   z = (x-y)/(x+y)^2,
#'   z_err = propagate_error((x-y)/(x+y)^2, dx = x_err, dy = y_err)
#' )
propagate_error <- function(eq, ..., quiet = FALSE) {
  # expression, variables and error terms
  eq_quo <- enquo(eq)
  eq_expr <- quo_get_expr(eq_quo)
  eq_env <- quo_get_env(eq_quo)
  variables <- all.vars(eq_expr)
  errors <- paste0("d", variables)
  error_terms <- rlang::enexprs(...)
  
  # check error terms are provided
  if (length(missing_terms <- setdiff(errors, names(error_terms))) > 0) {
    sprintf(
      "missing error term(s): '%s'. All terms must be defined but can be set to 0 to omit individual errors.",
      paste(missing_terms, collapse = "', '")
    ) %>%
      stop(call. = FALSE)
  } else if (length(extra_terms <- setdiff(names(error_terms), errors)) > 0) {
    sprintf(
      "extra error term(s): '%s'. These terms do not match the error terms expected for the formula.",
      paste(extra_terms, collapse = "', '")
    ) %>%
      stop(call. = FALSE)
  }
  
  # construct propagated error
  derivatives <- map(variables, ~D(eq_expr, .x))
  
  # terms for each variables: (df/dx * err)^2
  prop_error_terms <- map2(derivatives, error_terms[errors], ~rlang::expr(((!!.x) * (!!.y))^2))
  
  # sum
  sum_expr <- prop_error_terms[[1]]
  for (prop_error_term in prop_error_terms[-1]) sum_expr <- rlang::expr(!!sum_expr + !!prop_error_term)
  
  # sqrt
  final_expr <- rlang::expr(sqrt(!!sum_expr))
  
  # info
  info <- sprintf("Info: calculating propagated error for formula '%s' using equation '%s'",
                  paste0(deparse(eq_expr, width = 500), collapse=""), 
                  paste0(deparse(final_expr, width = 500), collapse="")) 
  if (!quiet) message(info)
  rlang::eval_tidy(final_expr, eq_env)
}