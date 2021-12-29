print.LinRegProblem <- function(x, ...) {
  cat("A linear regression problem.\n")

  # Index of last x column in dataframe
  last_x <- ncol(x$plot_data) - 3
  print(stats::setNames(x$plot_data, c(
    paste0(colnames(x$plot_data)[1:last_x], " (X", 1:last_x, ")"),
    colnames(x$plot_data)[(last_x + 1):ncol(x$plot_data)]
  )))
  cat("\n")
  cat("Betas \n")
  print(x$betas)
  cat("\nStatistics \n")
  atts <-
    setdiff(
      names(x),
      c("X", "Y", "Y_hat", "fitted", "resid", "betas", "plot_data")
    )
  pp_list(lapply(atts, get, envir = x), atts)
}

print.SSRDecomp <- function(x, ...) {
  cat("Sum of squared residuals decomposition", sep = "\n")
  cat("Reduced model:", attr(x, "reduced_model"), "\n")
  cat("Full model:", attr(x, "full_model"), "\n")
  pp_list(x, replace(
    names(x),
    which(names(x) == "R_squared_partial"),
    attr(x, "partition")
  ))
}

#' @title Update a \code{LinRegProblem} Object with a New Formula
#' @param object A \code{LinRegProblem} object.
#' @param new_formula A new model formula. \code{\link{model.matrix}} will attempt to evaluate it in the context of the existing object's column names.
#' @param ... Additional arguments to \code{\link{LinRegProblem}} to use when creating the new object, or to pass on to methods.
#' @details This function refits an existing \code{LinRegProblem} object to a new formula.The formula must be valid in the context of existing column names.
#'
#' @return If successful, the updated \code{LinRegProblem} object.
#' @export
#'
#' @examples
#' X <- LinRegProblem(mtcars[, -1], mtcars[, 1], Y_transform = scale)
#' update(X, Y ~ X1 + X2)
update.LinRegProblem <- function(object, new_formula, ...) {
  stopifnot("LinRegProblem" %in% class(object))
  new_X <-
    object$plot_data[, setdiff(colnames(object$plot_data), c("Y_hat", "resid")), drop = FALSE]
  new_Y <- object$plot_data[, "Y", drop = FALSE]
  new_X <- stats::model.matrix(new_formula, new_X)
  new_LinRegProblem(as.matrix(new_X), as.matrix(new_Y), ...)
}

anova.LinRegProblem <- function(object, ...) {
  anova_table <- with(object, matrix(c(
    SSR, p - 1, MSR,
    SSE, n - p, MSE,
    SSTO, n - 1, NA_real_,
    n * Y_bar^2, 1, NA_real_,
    sum(Y^2), n, NA_real_
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(
    c("Regression", "Error", "Total", "Mean-corrected", "Total uncorrected"),
    c("SS", "df", "MS")
  )
  ))

  cat("ANOVA Decomposition", sep = "\n")
  print_table(anova_table, upper_left = "Source of Variation")
}
registerS3method("print", "LinRegProblem", print.LinRegProblem)

registerS3method("print", "SSRDecomp", print.SSRDecomp)

registerS3method("update", "LinRegProblem", update.LinRegProblem)
registerS3method("anova", "LinRegProblem", anova.LinRegProblem)
