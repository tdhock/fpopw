% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compression.R
\name{uncompress.smt}
\alias{uncompress.smt}
\title{decompress.smt}
\usage{
uncompress.smt(smt.CP, vec.rep)
}
\arguments{
\item{smt.CP}{smoothed and compressed profile}

\item{vec.rep}{weights to use for decompression}
}
\value{
a vector to replicate duplicated datapoints
}
\description{
vector to decompress a compressed smoothed profile (a call to rep)
}
