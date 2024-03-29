% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operations.R
\name{ci}
\alias{ci}
\title{Compute the Bounds of a Confidence Interval}
\usage{
ci(center, bound)
}
\arguments{
\item{center}{Vector of interval centers}

\item{bound}{\eqn{\pm} boundary to define the interval}
}
\value{
\eqn{n \times 3} matrix of interval boundaries with the center, where \eqn{n} is the length of \code{center}.
}
\description{
Compute the Bounds of a Confidence Interval
}
\details{
This convenience function creates intervals from vectors of centers and boundary distances.
}
\examples{
ci(c(1.2, 3.4, 5), 1.7)
}
