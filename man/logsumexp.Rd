% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_helper.R
\name{logsumexp}
\alias{logsumexp}
\title{Logsumexp trick}
\usage{
logsumexp(log_vect)
}
\arguments{
\item{log_vect}{Input numerical vector.}
}
\value{
A probability vector, each element in [0,1] and all sum to 1
}
\description{
`logsumexp()` computes the softmax of a vector, using the logsumexp trick
}
