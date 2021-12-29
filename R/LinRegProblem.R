#' @title Construct a LinRegProblem Instance
#'
#' @param X  An \eqn{n \times m} design matrix, optionally with an intercept column present.
#' @param Y  An \eqn{n}-length vector or \eqn{n \times 1} matrix representing the response.
#' @param add_intercept Logical. Should an intercept column (a column of 1s) be added to \code{X}?
#' @param X_transform One-argument transformation function to apply to \code{X} before fitting the model.
#' Defaults to \code{\link{identity}}.
#' @param Y_transform One-argument transformation function to apply to \code{Y} before fitting the model. Defaults to \code{\link{identity}}.
#'
#' @details This is the constructor function for the \code{LinRegProblem} class. Users
#' should not use this function and should instead call \code{\link{LinRegProblem}}.
#'
#' @return A \code{LinRegProblem} object.
#' @export
#'
new_LinRegProblem <- function(X, Y, add_intercept = TRUE, X_transform = identity, Y_transform = identity) {
  validate_LinRegProblem(X, Y)
  X_transform <- match.fun(X_transform)
  Y_transform <- match.fun(Y_transform)
  intercept <- NULL

  instance <- within(list(), {
    X <- apply(X, MARGIN = 2, FUN = X_transform)
    Y <- apply(Y, MARGIN = 2, FUN = Y_transform)
    if (add_intercept) {
      X <- cbind(1, X)
    }
    hat <- unname(X %*% pseudoinverse(X))
    n <- nrow(Y)
    p <- ncol(X)
    # Ignored if false
    colnames(X) <- rep("", p)
    intercept <- has_intercept(X)
    if (intercept) colnames(X)[1] <- "(Intercept)"
    colnames(X)[(1 + intercept):p] <- default_colnames(X[, (1 + intercept):p, drop = FALSE])
    betas <- simple_ols(X, Y)
    Y_hat <- simple_fit(X, betas)
    resid <- unname(Y - Y_hat)
    colnames(Y) <- default_colnames(Y, "Y")
    colnames(Y_hat) <- default_colnames(Y_hat, "Y_hat")
    colnames(resid) <- default_colnames(resid, "resid")
    X_bar <- colMeans(X[, drop_if(intercept), drop = FALSE])
    Y_bar <- mean(Y)
    SSR <- sum_square(Y_hat, Y_bar)
    SSE <- sum((resid)^2)
    SSTO <- SSE + SSR
    MSE <- SSE / (n - p)
    MSR <- SSR / (p - 1)
    # Getting (X^TX)^-1 from the pseudoinverse
    covariance <- local({
      SVD <- svd(X, nv = ncol(X))
      MSE * SVD$v %*% ((ifelse(abs(SVD$d) < 1e-6, 0, 1 / (SVD$d^2))) * t(SVD$v))
    })
    r <- stats::cor(X[, drop_if(intercept)], Y)
    R_squared <- SSR / SSTO
    R_squared_adj <- 1 - (SSE / (n - p)) / (SSTO / (n - 1))
    leverage <- diag(hat)
    resid_var <- MSE * (1 - leverage)
    semistudentized <- resid / sqrt(resid_var)
    studentized_deleted <- resid * sqrt((n - p - 1) / (SSE * (1 - leverage) - resid^2))
    DFFITS <- studentized_deleted * sqrt(leverage / (1 - leverage))
    cook_sd <- (resid^2 / (p * MSE)) * (leverage / ((1 - leverage)^2))
    VIF <- tryCatch(diag(solve(stats::cor(X[, drop_if(intercept)]))), error = function(e) rep(NA_real_, p))
    plot_data <- data.frame(X[, drop_if(intercept), drop = FALSE], Y, Y_hat, resid)
  })

  instance <- as.environment(instance)
  parent.env(instance) <- globalenv()
  structure(
    instance,
    class = "LinRegProblem",
    predictors  = colnames(instance$X),
    response = colnames(instance$Y),
    intercept = intercept,
    multiple = ncol(instance$X) > 1,
    X_transform = deparse(substitute(X_transform)),
    Y_transform = deparse(substitute(Y_transform))
  )
}

validate_LinRegProblem <- function(X, Y) {
  allowed_types <- c("integer", "double")
  validate_same_length(X, Y)
  stopifnot(
    "X and Y must both be either integer or double" = typeof(X) %in% allowed_types && typeof(Y) %in% allowed_types,
    "Response must have only one column" = NCOL(Y) == 1
  )
}

#' @title Create a New  \code{LinRegProblem} Object
#' @inheritParams new_LinRegProblem
#' @details This function constructs a new \code{LinRegProblem} object from its arguments. Instances of the class contain the design matrix, response vector, fitted values, and estimated coefficients for a linear model. They also contain many relevant model statistics, such as \eqn{DFFITS} or the covariance matrix of the coefficients.
#' Members of the class are environments and therefore have reference semantics. Their main use is to illustrate the process of linear modeling and gather useful data about given models. By passing a \code{LinRegProblem} object as the \code{envir} argument to \code{sub_call}, you can evaluate arbitrary expressions in the context of the object, making it easy to do textbook problems.
#'
#' @return \code{LinRegProblem} object
#' @export
#'
#' @examples
#' LinRegProblem(mtcars[, -1], mtcars[, 1], Y_transform = scale)
LinRegProblem <- function(X, Y, add_intercept = TRUE, X_transform = identity, Y_transform = identity) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  new_LinRegProblem(X, Y, add_intercept, X_transform, Y_transform)
}
