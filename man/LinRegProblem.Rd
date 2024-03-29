% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinRegProblem.R
\name{LinRegProblem}
\alias{LinRegProblem}
\title{Create a New  \code{LinRegProblem} Object}
\usage{
LinRegProblem(
  X,
  Y,
  add_intercept = TRUE,
  X_transform = identity,
  Y_transform = identity
)
}
\arguments{
\item{X}{An \eqn{n \times m} design matrix, optionally with an intercept column present.}

\item{Y}{An \eqn{n}-length vector or \eqn{n \times 1} matrix representing the response.}

\item{add_intercept}{Logical. Should an intercept column (a column of 1s) be added to \code{X}?}

\item{X_transform}{One-argument transformation function to apply to \code{X} before fitting the model.
Defaults to \code{\link{identity}}.}

\item{Y_transform}{One-argument transformation function to apply to \code{Y} before fitting the model. Defaults to \code{\link{identity}}.}
}
\value{
\code{LinRegProblem} object
}
\description{
Create a New  \code{LinRegProblem} Object
}
\details{
This function constructs a new \code{LinRegProblem} object from its arguments. Instances of the class contain the design matrix, response vector, fitted values, and estimated coefficients for a linear model. They also contain many relevant model statistics, such as \eqn{DFFITS} or the covariance matrix of the coefficients.
Members of the class are environments and therefore have reference semantics. Their main use is to illustrate the process of linear modeling and gather useful data about given models. By passing a \code{LinRegProblem} object as the \code{envir} argument to \code{sub_call}, you can evaluate arbitrary expressions in the context of the object, making it easy to do textbook problems.
}
\examples{
LinRegProblem(mtcars[, -1], mtcars[, 1], Y_transform = scale)
}
