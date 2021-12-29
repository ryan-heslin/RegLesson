#' @title Center a Vector or Matrix
#' @param X Numeric vector or matrix
#' @details This function subtracts the mean from each element of a vector or matrix.
#'
#' @return Centered vector or matrix
#' @export
#'
#' @examples
#' center(1:10)
center <- function(X) {
  X - mean(X)
}


#' @title Compute the Total Sum of Squares
#' @param X Vector or matrix.
#' @details Computes the sum of the squared distances of a vector or matrix's elements from its mean.
#'
#' @return Scalar representing the \eqn{SSTO}.
#' @export
#'
#' @examples
#' with(iris, ssto(cbind(Sepal.Width, Sepal.Length)))
ssto <- function(X) {
  (t(X) %*% X) - (t(X) %*% matrix(1, nrow = NROW(X), ncol = NROW(X)) %*% X)
}

#' @title Compute Two-sided \eqn{p}-Values under the Normal Distribution
#' @param x A numeric vector
#' @param deg_free Degrees of freedom under the null hypothesis distribution
#' @details This convenience function computes the probabilities of a vector of values under the distribution assumed by the null hypothesis. I am legally obliged to note this is not remotely the same thing as the probability of a genuine effect.
#'
#' @return Vector of computed \eqn{p}-values
#' @export
#'
#'
p_value <- function(x, deg_free) {
  2 * stats::pt(-abs(x), deg_free)
}

#' @title Compute the Sum of Squares
#' @param x Numeric vector
#' @param y Numeric scalar vector of the same length as \code{x}. Defaults to
#' the mean of \code{x}.
#' @details This function computes the sum of squared differences between \code{x} and \code{y}
#'
#' @return Sum of squares of the deviation vector.
#' @export
#'
#'
sum_square <- function(x, y = mean(x)) {
  validate_same_length(x, y, .recycle_ok = TRUE)
  sum((x - y)^2)
}

#' @title Compute the Bounds of a Confidence Interval
#' @param center Vector of interval centers
#' @param bound  \eqn{\pm} boundary to define the interval
#' @details This convenience function creates intervals from vectors of centers and boundary distances.
#'
#' @return  \eqn{n \times 3} matrix of interval boundaries with the center, where \eqn{n} is the length of \code{center}.
#' @export
#'
#' @examples
#' ci(c(1.2, 3.4, 5), 1.7)
ci <-
  function(center, bound) {
    cbind(
      lower = center - bound,
      estimate = center,
      upper = center + bound
    )
  }

#' @inheritParams new_LinRegProblem
#'
#' @title Conduct Simple Linear Regression
#' @details This function finds the minimum-norm least squares solution of \eqn{X^T X \beta = X^T Y}.
#'
#' @return Vector of estimated beta coefficients.
#' @export
#'
#'
simple_ols <- function(X, Y, add_intercept = FALSE) {
  validate_same_length(X, Y)
  if (add_intercept) {
    X <- cbind(1, X)
  }
  # Find minimal solution if need be
  betas <- as.vector(pseudoinverse(X) %*% Y)
  # Crude hack to check if intercept column present
  stats::setNames(betas, paste0("beta", seq_along(betas) - has_intercept(X)))
}

#' @inheritParams new_LinRegProblem
#'
#' @title Fit a Linear Model
#' @param betas Vector of coefficients for the model.
#' @details This function finds the least-squares estimates given a design matrix and bet coefficients.
#'
#' @return Vector of least-squares estimates.
#' @export
#'
#'
simple_fit <- function(X, betas) {
  stopifnot(NCOL(X) == NROW(betas))
  X %*% as.matrix(betas)
}


#' @title Compute the Mean Squared Error
#' @param Y Response vector.
#' @param Y_hat Vector of predictions.
#' @param n Number of rows in design matrix.
#' @param p Number of columns in design matrix.
#' @details This function computes the mean squared error, given coefficient and prediction vectors and the parameters for computing degrees of freedom.
#'
#' @return Estimated mean squared error.
#' @export
#'
#'
mse <- function(Y, Y_hat, n = NROW(Y), p) {
  sum_square(Y, Y_hat) / (n - p)
}

#' @title Compute Confidence Intervals for Estimated Coefficients
#' @param betas Vector of estimated beta coefficients. \eqn{beta+0} is valid.
#' @param MSE Estimated mean squared error for the model
#' @param n Number of rows of the design matrix
#' @param p Number of parameters in the model, counting \eqn{beta_0}.
#' @param covariance Covariance matrix of the estimated coefficients
#' @param alpha Alpha level for the test
#' @param simultaneous Logical. Apply Bonferroni correction for simultaneous intervals? Default \code{FALSE}.
#' @details This function computes confidence intervals for the estimated beta coefficients of a linear model, optionally simultaneously.
#'
#' @return A \eqn{q \times 3} matrix of the estimated coefficients, with lower and upper bounds for the interval.
#' @export
#'
#'
confint_betas <-
  function(betas,
           MSE,
           n,
           p = NROW(covariance),
           covariance,
           alpha = .05,
           simultaneous = FALSE) {
    stat <-
      stats::qt(1 - alpha / (2 * ifelse(simultaneous, p, 1)), df = n - p) * sqrt(diag(covariance))
    ci(betas, stat)
  }
# Compute the confidence interval for the mean of
# Y_h (i.e., an existing observation)
# Set m to Inf to get CI for E of a value for new
# X - same as for existing value
# X_h is m x p, betas p x 1, covariance p x p

#'
#'
#' @title Compute Confidence Intervals for Regression Estimates
#' @param X_h \eqn{h \times p} matrix, or single \eqn{h} length vector, containing data for observations. Be sure to include the
#' intercept column if present.
#' @param betas Vector of estimated coefficients.
#' @param MSE Estimated mean squared error.
#' @param covariance Covariance matrix of the estimated coefficients.
#' @param n Number of observations in the design matrix. Computed if omitted.
#' @param p Number of parameters in the model. Computed if omitted.
#' @param m Number of new observations to predict.
#' @param alpha Alpha level for constructing confidence intervals.
#' @param new Logical. Are the values new observations, or were they used to fit the model?
#' Confidence intervals are narrower in the latter case. Default \code{FALSE}.
#' @param simultaneous Logical. Should simultaneous Bonferroni-corrected intervals be constructed instead of separate intervals? Default \code{FALSE}.
#' @details This function finds the predicted values for given observations predicted
#' by an estimated linear model. It also computes confidence intervals, the width
#' varying by whether the observations are new and/or whether simultaneous intervals are desired.
#'
#' @return \eqn{h \times 3} matrix of predictions with interval bounds.
#' @export
#'
#'
confint_EY_h <-
  function(X_h,
           betas,
           MSE,
           covariance,
           n,
           p = NCOL(X_h),
           m = 1,
           alpha = .05,
           new = FALSE,
           simultaneous = FALSE) {
    Y_h_hat <- c(X_h %*% betas)

    # Add MSE/m if new levels of X
    s_Y_h_hat <-
      diag(X_h %*% covariance %*% t(X_h)) + (MSE * ((1 / m) * new))
    bound <-
      stats::qt(1 - alpha / (2 * ifelse(simultaneous, NROW(Y_h_hat), 1)), df = n - p) * sqrt(s_Y_h_hat)
    ci(Y_h_hat, bound)
  }

#' @title Conduct the \eqn{F} Test for Existence of a Regression Relation
#' @param MSR Estimated mean square of regression for the model
#' @param MSE Estimated mean squared error for the model
#' @param n Number of rows of the data matrix
#' @param p Number of parameters in the model
#' @param alpha Alpha level to use for the test.
#' @details This function tests the null hypothesis that no
#' regression relation exists, given the relevant mean squares and degrees of freedom.
#'
#' @return Value of the computed \eqn{F} statistic
#' @export
#'
#'
regression_relation <- function(MSR, MSE, n, p, alpha = .05) {
  stat <- MSR / MSE
  announce_result(
    stat,
    stats::qf(1 - alpha, p - 1, n - p, lower.tail = FALSE),
    "F test of regression relation",
    alpha
  )
  stat
}

# Given positional subscripts of beta coefficients for full and reduced
# models on a data matrix and response, computes SSE

#' @title Decompose Sums of Squares within a Linear Model
#' @param X  A design matrix for a linear model
#' @param Y  A response vector for a linear model
#' @param reduced_betas Indices for beta coefficients, indexed by 0 (e.g., 0 for \eqn{\beta_0}, etc. ).
#' @param full_betas  Indices for coefficients of the full model. Defaults to all coefficients.
#' @param intercept_fitted Logical. Does the model include an intercept? Default \code{TRUE}.
#' @details This function decomposes the sums of squares within a reduced and full model and computes the full model's coefficient of partial determination: the portion of \eqn{SSE} of the reduced model it explains.
#'
#' @return Object of \code{SSRDecomp} class. It contains the \eqn{SSR} and \eqn{SSE} of the full and reduced models, along with the coefficient of partial correlation.
#' @export
#'
#'
SSRDecomp <-
  function(X,
           Y,
           reduced_betas,
           full_betas,
           intercept_fitted = TRUE) {
    # Convert betas to valid indices accounting for intercept
    # at index 1 - 2 -> beta 1, etc.
    if (!intercept_fitted) {
      X <- cbind(1, X)
    }
    model_reduced_betas <- reduced_betas + 1

    if (missing(full_betas)) {
      full_betas <- seq(0, NCOL(X) - 1)
    }
    model_full_betas <- full_betas + 1

    validate_same_length(X, Y)
    invisible(sapply(
      list(model_reduced_betas, model_full_betas),
      validate_subscript,
      x = X,
      negative_ok = FALSE
    ))
    validate_subset(model_full_betas, model_reduced_betas, improper_ok = FALSE)

    Y_bar <- mean(Y)
    reduced <-
      simple_fit(X[, model_reduced_betas], simple_ols(X[, model_reduced_betas], Y, add_intercept = FALSE))
    full <-
      simple_fit(X[, model_full_betas], simple_ols(X[, model_full_betas], Y, add_intercept = FALSE))
    # Convert beta subscript vectors back to zero-index for notation
    partition <-
      paste0(
        "R_squared_Y",
        paste(sort(setdiff(
          full_betas[full_betas != 0], reduced_betas[reduced_betas != 0]
        )), collapse = ""),
        "|",
        paste(sort(reduced_betas[reduced_betas != 0]), collapse = "")
      )

    # Create instance of SSRDecomp, with relevant sums
    # of squares and attributes representing model
    # specification
    out <- within(list(), {
      SSR_full <- sum_square(full, Y_bar)
      SSR_reduced <- sum_square(reduced, Y_bar)
      SSE_full <- sum_square(Y, full)
      SSE_reduced <- sum_square(Y, reduced)
      SSR_marginal <- SSE_reduced - SSE_full
      R_squared_partial <- SSR_marginal / SSE_reduced
    })
    structure(
      out,
      class = "SSRDecomp",
      reduced_betas = reduced_betas,
      full_betas = full_betas,
      reduced_model = paste("Y =", paste(
        sprintf("b%dX%d", reduced_betas, reduced_betas),
        collapse = " + "
      )),
      full_model = paste("Y =", paste(
        sprintf("b%dX%d", full_betas, full_betas),
        collapse = " + "
      )),
      partition = partition
    )
  }

#' @inheritParams new_LinRegProblem
#'
#' @title Conduct the Correlation Transformation
#' @details Applies the correlation transformation to \code{X} ,
#' standardizing each variable.
#'
#' @return The transformed data matrix.
#' @export
#'
#'
correlate_transform <- function(X) {
  (center(X) / stats::sd(X)) / sqrt(length(X) - 1)
}


#' @inheritParams new_LinRegProblem
#'
#' @title Conduct the Generalized \eqn{F} Test
#' @param reduced_betas Vector of indices of beta coefficients for the reduced model, indexed from zero for \eqn{\beta_0}, 1 for \eqn{\beta_1}, etc.
#' @param full_betas Vector of indices of beta coefficients for the full model, indexed the same way.
#' @param n Number of rows in the design matrix.
#' @param p Number of columns in the design matrix.
#' @param intercept_fitted Logical. Was an intercept included in the model? Default \code{TRUE}.
#' @param alpha Alpha level of the significance test.
#' @details This function conducts the generalized \eqn{F} test comparing a model with
#' \eqn{q} parameters to one with \eqn{p} parameters. The extra degrees of freedom are therefore \eqn{q-p}, where \eqn{q -p \geq 1}.
#'
#' @return \eqn{F} statistic of the test.
#' @export
#'
#'
extra_SS <-
  function(X,
           Y,
           reduced_betas,
           full_betas,
           n = NROW(Y),
           p = NCOL(X),
           intercept_fitted = TRUE,
           alpha = .05) {
    stopifnot(NCOL(X) > 1)
    q <- NROW(reduced_betas) + intercept_fitted
    decomp <-
      SSRDecomp(X, Y, reduced_betas, full_betas, intercept_fitted)
    # Remember reduced betas are RETAINED, not dropped
    cat("Betas tested:\n")
    pp_list(paste0("b", setdiff(full_betas, reduced_betas)))
    stat <-
      ((decomp$SSE_reduced - decomp$SSE_full) / (p - q)) / (decomp$SSE_full / (n - p))
    announce_result(stat, stats::qf(1 - alpha, p - q, n - p, lower.tail = FALSE), "F-test of beta signficance", alpha)
    stat
  }


#' @inheritParams new_LinRegProblem
#'
#' @title Conduct the \eqn{F} Test for Lack of Fit
#' @param Y_hat Vector of fitted values
#' @param resid Vector of residuals
#' @param n Number of rows in the design matrix. Computed if not provided.
#' @param p Number of columns in the design matrix. Computed if not provided.
#' @param alpha Alpha level for the significance test.
#' @details This function carries out the \eqn{F} test for lack of fit, which
#' essentially compares the model estimates to the local means of the response
#' at each distinct observation. The test requires duplicate observations.
#'
#' @return Computed \eqn{F} statistic for the test.
#' @export
#'
#'
lack_of_fit <-
  function(X,
           Y,
           Y_hat,
           resid = sum_square(Y, Y_hat),
           n = NROW(Y),
           p = NCOL(X),
           alpha = .05) {
    validate_same_length(X, Y, Y_hat)
    X_j <-
      do.call(interaction, c(asplit(X, MARGIN = 2), drop = TRUE))
    unique_xes <- length(unique(X_j))
    stopifnot("No repeat X values in data" = unique_xes < n, unique_xes > 2)
    Y_bar_j <- tapply(Y, X_j, mean)
    SSPE <- sum(sapply(split(Y, X_j), sum_square))
    stat <-
      (resid - SSPE) / (unique_xes - 2) / (SSPE / (n - unique_xes))
    announce_result(
      stat,
      stats::qf(1 - alpha, unique_xes - p, n - unique_xes, lower.tail = FALSE),
      "F test of lack of fit",
      alpha
    )
  }
# Be sure to add intercept column if present!

#' @title Compute Working-Hoteling Confidence Bands
#' @param X Design matrix of the linear model.
#' @param X_h Matrix of observations to fit the interval for.
#' @param betas Estimated regression coefficients.
#' @param Y_h \eqn{1 \times h}
#' @param covariance \eqn{p \times p} covariance matrix of beta coefficients.
#' @param alpha Significance level for the interval.
#' @param n Number of rows in the data matrix.
#' @param p Number of parameters in the model.
#' @param MSE Estimated mean squared error.
#' @details Given fitted values or enough data to compute them, computes the boundaries of the Working-Hoteling interval estimates.
#'
#' @return \eqn{h \times 3} matrix of interval bounds
#' @export
#'
#'
Working_Hoteling <-
  function(X,
           X_h,
           betas,
           Y_h = simple_fit(X_h, betas),
           covariance,
           alpha = .05,
           n,
           p,
           MSE) {
    validate_same_length(X_h, Y_h)
    s_Y_h_hat <-
      diag(X_h %*% covariance %*% t(X_h))
    bound <- sqrt(p * stats::qf(1 - alpha, p, n - p, lower.tail = FALSE) * s_Y_h_hat)
    t(mapply(ci, Y_h, bound))
  }

# Theoretical normal distribution of ranked values

#' @title Find the Expected Values of Residuals Assuming Normality
#' @param ranks Ranks of model residuals.
#' @param mse Estimated mean squared error of the model.
#' @details This function computes the expected values of model residuals under a normal distribution, as the ordinary model assumptions say they should. The results can be used in normality tests. See _ALSM_, p. 110-111.
#'
#' @return Vector of the expected value approximation described above.
#' @export
#'
#'
expected_normal <-
  function(ranks, mse) {
    sqrt(mse) * stats::qnorm((ranks - .375) / (length(ranks) + .25))
  }

#' @title Conduct the Brown-Forsythe Test of Equal Variances on Residuals
#' @param resid vector of model residuals.
#' @param groups Factor or logical of the same length as \code{resid} denoting groups into which to split the residuals. Defaults to a logical vector indicating whether each residual is below the median residual.
#' @param n Number of rows in the design matrix.
#' @param p Number of columns in the design matrix..
#' @param alpha Alpha level for the significance test.
#' @details This function performs the Brown-Forsythe test of equal variance on specified groups of model residuals. A significant result provides evidence of unequal variances.
#'
#' @return Value of the Brown-Forsythe statistic.
#' @export
#'
#'
brown_forsythe <-
  function(resid,
           groups = resid <= stats::median(resid),
           n = NROW(resid),
           p,
           alpha = .05) {
    validate_same_length(resid, groups)
    resid <-
      lapply(split(resid, groups), function(x) {
        abs(x - stats::median(x))
      })
    ns <- lengths(resid)
    # pooled SE
    s_squared <-
      sum(sapply(resid, function(x) {
        sum_square(x, mean(x))
      })) / (sum(ns) - p)
    stat <-
      abs(diff(sapply(resid, mean)) / sqrt(s_squared * sum(1 / ns)))
    announce_result(
      stat,
      stats::qt(1 - alpha / 2, df = sum(ns) - p),
      "Brown-Forsythe test of constant variance",
      alpha
    )
    unname(stat)
  }

#' @title Conduct the Breusch-Pagan Test of Heteroscedascity
#' @param X Design matrix of a linear model.
#' @param resid Vector of residuals from the fitted model.
#' @param SSE Sum of squared errors of the fitted model.
#' @param n Number of rows of the design matrix. Computed if not provided.
#' @param q Number of parameters in the model.
#' @param alpha Alpha level of the significance test.
#' @param has_intercept Logical. Was the model fitted with an intercept?
#' @details This function carries out the Breusch-Pagan test, which conducts a chi-squared test to assess whether the residuals have constant variance.
#'
#' @return Value of the Breusch-Pagan statistic.
#' @export
#'
#'
breusch_pagan <-
  function(X,
           resid,
           SSE,
           n = NROW(X),
           q = seq_len(NCOL(X)),
           alpha = .05,
           has_intercept = TRUE) {
    validate_same_length(X, resid)
    if (!has_intercept) X <- cbind(1, X)
    validate_subscript(X, q, 2)
    resid_fit <- LinRegProblem(X[, q], resid^2, add_intercept = !has_intercept(X))
    stat <- (resid_fit$SSR / (length(q) - 1)) / ((SSE / n)^2)
    announce_result(
      stat,
      stats::qchisq(1 - alpha, length(q) - 1),
      "Breusch-Pagan test for constancy of error variance",
      alpha
    )
    stat
  }


#' @inheritParams new_LinRegProblem
#'
#' @title Compare the Box-Cox Transformation at Different Levels of \eqn{\lambda}
#' @param lambdas Vector of values of \eqn{\lambda}. The transformation entails raising \eqn{Y} to the power of \eqn{\lambda}, or logging it in the case of \eqn{\lambda = 0}.
#' @param n Number of rows in the design matrix.
#' @details This function assess the Box-Cox transformation at different values of
#' \eqn{\lambda}, the tuning parameter. The best corresponds to the lowest model \eqn{SSE}.
#'
#' @return Two-column matrix with the vector of \eqn{\lambda} in one column and the estimated \eqn{SSE} in the other.
#' @export
#'
#'
box_cox <- function(X, Y, lambdas, n = NROW(Y)) {
  validate_same_length(X, Y)
  K2 <- prod(Y)^(1 / n)

  bc_standardize <- function(lambda, Y, K2) {
    if (lambda == 0) {
      K2 * log(Y)
    } else {
      (1 / (lambda * K2^(lambda - 1))) * (Y^(lambda) - 1)
    }
  }
  standardizations <-
    lapply(lambdas, bc_standardize, Y = Y, K2 = K2)
  fits <-
    lapply(standardizations, function(x) {
      simple_fit(X, betas = simple_ols(X, x))
    })
  SSEs <- mapply(sum_square, fits, standardizations)
  matrix(c(lambdas, SSEs),
    nrow = length(lambdas),
    dimnames = list(NULL, c("lambda", "SSE"))
  )
}

#' @title Solve Problems Given a Linear Model
#' @param instance A \code{\link{LinRegProblem}} instance, which is passed as the \code{envir} argument to \code{\link{substitute_call}}.
#' @param ... Expressions or symbols to evaluate after substitution in the context of \code{instance}.
#' @param more_subs Passed directly to \code{\link{substitute_call}}.
#' @param .print Logical. Should the results be printed (useful if the output goes in an Rmarkdown file), or returned as a list?
#' @details This function takes arbitrary expressions via \code{...} and evaluates them in the context of a \code{\link{LinRegProblem}} instance. The results are simply printed, if \code{print} is \code{TRUE}, and otherwise returned as a list.
#' Because \code{\link{substitute_call}} will complete partial calls, this allows for compact specification of solutions.
#'
#' @return \code{NULL} if \code{.do_print} is \code{TRUE}, otherwise a list containing results of the expressions.
#' @export
#'
#' @examples
#' X <- LinRegProblem(mtcars[, -1], mtcars[, 1], Y_transform = scale)
#' do_problems(X, regression_relation())
# Common_args must be quoted, others unquoted
do_problems <- function(instance, ..., more_subs = list(), .print = TRUE) {
  stopifnot("LinRegProblem" %in% class(instance))
  expressions <- as.list(substitute(list(...)))[-1L]
  substituted <- lapply(expressions, substitute_call,
    envir = instance,
    quote = FALSE, more_subs = more_subs
  )
  results <- lapply(substituted, eval, envir = instance)
  if (.print) {
    pp_list(results)
  } else {
    results
  }
}

#' @title Substitute Values from an Environment into a Call
#' @param sub_call R expression, quoted if \code{quote = FALSE}, unquoted otherwise. May also be a symbol or constant.
#' @param envir An environment from which to substitute values.
#' @param more_subs Optional named list of additional values to substitute. If any names are shared with \code{envir}, those in \code{more_subs} take precedence.
#' @param depth Used internally to track recursion depth. Do not specify this argument.
#' @param quote Logical. Should the call be quoted (i.e., it is passed as an unevaluated expression?).
#' @details This function recursively traverses an R call object (a list-like representation of the abstract syntax tree) and replaces names with corresponding values in \code{envir}. Mandatory non-default arguments to any function in the call are substituted the same way, even if omitted; thus, this function can
#' convert semantically invalid calls into valid ones. The modified call may still be invalid, however, if some mandatory arguments are unspecified and neither \code{envir} or \eqn{more_subs} contains substitutions for them.
#'
#' This function is meant to be used with calls to functions from this package, using a \code{\link{LinRegProblem}} instance as \code{envir}. Because most function arguments are also attributes of \code{LinRegProblem} objects, doing so will automatically yield a valid call.
#' @return Call object with all matched names substituted with corresponding names.
#' @export
#'
#'
substitute_call <-
  function(sub_call,
           envir,
           more_subs = list(),
           depth = 1,
           quote = TRUE) {
    # Quote call only on entering
    if (depth == 1) {
      if (quote) sub_call <- substitute(sub_call)
      more_subs <- as.environment(more_subs)
      parent.env(more_subs) <- envir
      # Eval and return if single symbol, just return if constant
      if (is.symbol(sub_call)) {
        return(eval(sub_call, more_subs))
      } else if (is.atomic(sub_call)) {
        return(sub_call)
      }
    }
    sub_call <- as.list(sub_call)
    current_fun <- eval(sub_call[[1]])
    # Take union of missing args lacking defaults and those
    # equal to missing symbol and substitute values in envir
    subscript <- drop_if(length(sub_call) > 1)
    missing <- names(sub_call[subscript][sapply(sub_call[subscript], is_missing_symbol)])
    no_defaults <-
      setdiff(get_default_args(current_fun, invert = TRUE), names(sub_call[subscript]))

    # more_subs preempts default args as well
    to_substitute <- union(
      missing,
      no_defaults
    ) |>
      union(intersect(names(more_subs), get_default_args(current_fun)))
    available <- union(names(envir), names(more_subs))
    sub_call[intersect(to_substitute, available)] <- lapply(intersect(to_substitute, available), get, pos = more_subs)
    calls <- sapply(sub_call[subscript], is.call)
    # Recurse over calls in arguments of current call
    sub_call[subscript][calls] <-
      lapply(sub_call[subscript][calls],
        substitute_call,
        more_subs = more_subs,
        envir = envir,
        depth = depth + 1
      )
    as.call(sub_call)
  }

#' @title Estimate Observations of \eqn{X} from a Given \code{Y}
#' @param Y Vector of \eqn{Y} observations used to fit the model
#' @param Y_0 Hypothetical \eqn{Y} observations on which to predict \code{X}
#' @param betas Estimated regression coefficients.
#' @param Y_bar Mean of the response vector \eqn{Y}. Will be computed if not provided.
#' @param X Design matrix used to fit a model.
#' @param m Number of new observations at each level of \code{Y_0} to predict.
#' @param n Number of rows in the design matrix.
#' @param p Number of parameters in the model.
#' @param MSE Estimated mean squared error.
#' @param alpha Alpha level for confidence intervals.
#' @param simultaneous Logical. Should Bonferroni-corrected simultaneous intervals be constructed for each prediction? Default \code{FALSE}.
#' @param intercept_fitted Logical. Does the model contain an intercept parameter? Default \code{FALSE}.
#' @details this function does inverse prediction, finding observations of \eqn{X} given observations of \eqn{Y}.
#'
#' @return \eqn{3 \times h} matrix of predictions with interval bounds, where \eqn{h} is n the number of observations predicted.
#' @export
#'
#'
calibrate <- function(Y, Y_0, betas,
                      Y_bar = mean(Y),
                      X,
                      m = 1,
                      n = NROW(Y),
                      p = NROW(betas),
                      MSE,
                      alpha = .05,
                      simultaneous = FALSE,
                      intercept_fitted = TRUE) {
  stopifnot("Intercept term cannot exist in model with p = 1" = intercept_fitted && p > 1)
  # Irritating special case
  if (intercept_fitted) {
    p <- length(betas) - 1
    beta_0 <- betas[1]
    betas <- betas[-1]
  } else {
    beta_0 <- 0
  }
  Ys <- diag(x = Y_0 - unname(beta_0), nrow = NROW(Y_0)) %*% matrix(1, ncol = p, nrow = NROW(Y_0))
  preds <- Ys %*% ifelse(betas == 0, 0, 1 / betas)
  # See https://ntrs.nasa.gov/api/citations/20110016499/downloads/20110016499.pdf
  s_X_0 <- (sqrt(MSE) / sum(betas^2)) * ((1 / m) + (1 / n) + (((preds - mean(X))^2) / sum_square(X)))
  stat <-
    stats::qt(1 - alpha / (2 * ifelse(simultaneous, p, 1)), df = n - p) * s_X_0
  ci(c(preds), stat)
}

#' @title Find the Pseudoinverse of a Matrix
#' @param X A numeric matrix.
#' @details The pseudoinverse of a matrix is \eqn{V \Sigma^+ U^T} - essentially, its inverse over its row space. It always exists even if no left or right inverse does.
#' This function is partly based on the \code{\link[MASS]{ginv}}, which does the same thing.
#'
#' @return The pseudoinverse of \eqn{X}.
#' @export
#'
#' @examples
#' pseudoinverse(matrix(rnorm(12), ncol = 3))
pseudoinverse <- function(X) {
  stopifnot(is.numeric(X))
  if (!is.matrix(X)) X <- as.matrix(X)
  eigens <- eigen(t(X) %*% X, symmetric = TRUE)
  sigmas <- sqrt(round(eigens$values, digits = 6))
  signif <- which(!is_zero(sigmas))
  if (!length(signif)) {
    return(matrix(0, nrow = ncol(X), ncol = nrow(X)))
  }
  # The "economy" SVD - (n x r) (r x r) (r x m)
  svd(X)$v[, signif, drop = FALSE] %*% ((1 / sigmas[signif]) * t(svd(X)$u[, signif, drop = FALSE]))
}

is_zero <- function(x) {
  sapply(x, function(x) isTRUE(all.equal(x, 0)))
}
