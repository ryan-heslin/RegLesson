current_envir <- function() {
  parent.frame()
}

has_intercept <- function(X) {
  NCOL(X) > 1 && all(X[, 1] == 1)
}

#' @title Name a List's Elements After Itself
#' @param ... Elements of the list. If unnamed, they are named after themselves
#' @details This function gathers a collection of objects in a list. Names are quoted directly from the expressions passed via \code{....}.
#'
#' @return A list containing the elements of \code{...}, named as described above.
#' @export
#'
#' @examples selfname_list(mtcars, iris, y = 1:5)
selfname_list <- function(...) {
  dots <- list(...)
  dots_exprs <- sapply(as.list(substitute(list(...))[-1]), deparse)
  if (is.null(names(dots))) {
    bad_names <- rep(TRUE, length(dots))
  } else {
    bad_names <- names(dots) %in% c("", NA_character_)
  }
  stats::setNames(dots, ifelse(bad_names, dots_exprs, names(dots)))
}

# Create vector of spaces to pad elements of a vector when printed so all
# have the same indentation.
pad_print <- function(x) {
  x <- as.character(x)
  # literal string NA, not the value
  x[is.na(x)] <- "NA"
  spaces <- max(nchar(x)) + 1 - nchar(x)
  sapply(spaces, function(y) {
    paste(rep(" ", y), collapse = "")
  })
}

announce_result <- function(stat, threshold, name, alpha, reject_above_threshold = TRUE) {
  validate_same_length(stat, threshold, name, alpha, reject_above_threshold, .length = 1)
  cat(name, paste("alpha =", alpha), paste("Test statistic:", stat),
    paste("Threshold:", threshold), paste("Decision:", ifelse((stat > threshold & reject_above_threshold) |
      (stat < threshold & !reject_above_threshold),
    "Reject null", "Fail to reject null"
    )),
    sep = "\n"
  )
}
# Print members of list, prepending by name if it exists
pp_list <- function(lst, print_names = names(lst)) {
  if (is.null(print_names)) {
    print_names <- rep("", length(lst))
  }
  validate_same_length(lst, print_names)
  spaces <- pad_print(print_names)
  nonempty <- nzchar(print_names)
  print_names[nonempty] <- paste0(print_names[nonempty], ":")
  print_names <- paste0(print_names, spaces)

  invisible(mapply(function(x, y) {
    cat(y)
    if (!is.null(attributes(x))) {
      cat("\n")
      print(x)
    } else {
      cat(x, "\n")
    }
  }, lst, print_names))
}
# If recycle_ok, 1-length vectors don't throw error
validate_same_length <-
  function(..., .length, .recycle_ok = FALSE) {
    dots_lengths <- unique(sapply(list(...), NROW))
    stopifnot((if (missing(.length)) length(dots_lengths) == 1L else dots_lengths == .length) ||
      (.recycle_ok &&
        length(dots_lengths) == 2L &&
        1L %in% dots_lengths))
  }

# Nicely print a matrix or array, justifying left and optionally
# including dimnames.
print_table <- function(x, upper_left = "", use_names = TRUE, blank_NA = TRUE) {
  if (use_names) {
    if (!is.null(rownames(x))) {
      x <- cbind(rownames(x), x)
      rownames(x) <- NULL
      rownames_added <- TRUE
    }
    if (!is.null(colnames(x))) {
      if (rownames_added) {
        colnames(x)[1] <- upper_left
      }
      x <- rbind(colnames(x), x)
      colnames(x) <- NULL
    }
  }
  if (is.character(x) && blank_NA) {
    x[is.na(x)] <- ""
  }

  spaces <- apply(x, MARGIN = 2, pad_print)
  x <- apply(x, MARGIN = 2, \(x) paste0(x, pad_print(x)))
  invisible(apply(x, MARGIN = 1, \(x) cat(x, "\n")))
}

validate_determined <- function(X) {
  tryCatch(
    solve(t(X) %*% X),
    error = function(e) {
      stop("Linear dependency in data matrix")
    }
  )
}

# Confirm vector only natural numbers
validate_natural <- function(x) {
  all(x > 0 && x %% 1 == 0)
}

# Confirm all elements of subscript vector correspond to valid indices
validate_subscript <- function(x, subscript, dimension = 1, negative_ok = TRUE, duplicate_ok = FALSE) {
  if (!negative_ok) stopifnot(min(subscript) > 0)
  if (!duplicate_ok) stopifnot(anyDuplicated(subscript) == 0)
  subscript <- abs(subscript)
  stopifnot(validate_natural(subscript))
  if (is.null(x_dim <- dim(x))) {
    stopifnot(dimension == 1, max(subscript) <= max(seq_along(x)))
  } else {
    stopifnot(
      validate_subscript(x_dim, dimension, negative_ok = FALSE),
      max(subscript) <= x_dim[dimension]
    )
  }
  invisible(TRUE)
}

# Confirm one set is (optionally proper) subset of another
validate_subset <- function(set, sub_set, improper_ok = TRUE) {
  stopifnot(all(intersect(sub_set, set) == sub_set) && (improper_ok || length(set) > length(sub_set)))
}

# Sets have no more than threshold items in intersect
validate_intersect_threshold <- function(A, B, threshold = 1L) {
  stopifnot(length(intersect(A, B)) < threshold)
}

#' @title Assign Generic Column Names to an Object
#' @param X An object that can receive a `colnames` attribute
#' @param name_base Stem for each default column name. Default `"X"`.
#' @details Old column names are retained unless they are the empty string. Any such column names are replaced by new default names numbered from 1. All columns are treated this way if `colnames(X)` is `NULL`.
#'
#' @return \code{X}, with its colnames attribute modified as described.
#' @export
#'
#' @examples
#' X <- matrix(1:16, nrow = 2)
#' default_colnames(X)
default_colnames <- function(X, name_base = "X") {
  stopifnot(is.matrix(X))
  old_colnames <- colnames(X)
  new_colnames <- paste0(name_base, 1:ncol(X))
  if (is.null(old_colnames) || all(old_colnames == "")) {
    if (ncol(X) == 1) {
      name_base
    } else {
      new_colnames
    }
    # Replace blank or NA colnames with defaults
  } else {
    ifelse(
      sapply(old_colnames, is.na, USE.NAMES = FALSE) |
        old_colnames == "",
      new_colnames,
      old_colnames
    )
  }
}

# Return negative index only if predicate is satisfied
drop_if <-
  function(predicate, drop_index = 1) {
    if (predicate) {
      -drop_index
    } else {
      substitute()
    }
  }

# Convert scraped problem data to LinRegProblem instances
data2problem <- function(data, add_intercept = TRUE) {
  tryCatch(
    LinRegProblem(
      X = data[, -1, drop = FALSE],
      Y = data[, 1, drop = FALSE],
      add_intercept = add_intercept
    ),
    error = function(e) {
      cat("Error", e$message, "\n")
    }
  )
}

# Detach environment in caller environment, update with new expression, reattach
update_env_mask <- function(envir, update_expr) {
  envir <- substitute(envir)
  update_expr <- substitute(update_expr)
  as.function(list(bquote({
    evalq(tryCatch(detach(.(envir)), error = function(e) {
      invisible(NULL)
    }), envir = parent.frame())

    evalq(.(update_expr), envir = parent.frame())
    evalq(attach(.(envir)), envir = parent.frame())
  })))
}

is_missing_symbol <- function(x) {
  identical(x, substitute())
}


#' @title Get a Function's Default Arguments
#' @param fun Function object, or the name of a function.
#' @param invert Logical. Obtain only non-default arguments?
#' @param ignore_dots Logical. Consider \code{...} (if present) a non-default argument?
#' @details This convenience function extracts the default arguments and values from a function, or the non-default arguments if \code{invert = TRUE}.
#'
#' @return Character vector of the names of the specified arguments, \code{NULL} if none meet the conditions.
#' @export
#'
#' @examples
#' get_default_args(lm)
#'
# Extract a function's arguments with default values, or without
get_default_args <- function(fun, invert = FALSE, ignore_dots = TRUE) {
  args <- formals(fun)
  if (ignore_dots) args <- args[names(args) != "..."]
  names(args[invert == sapply(args, is_missing_symbol)])
}
