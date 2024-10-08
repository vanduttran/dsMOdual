% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client_func.R
\name{federateComDim}
\alias{federateComDim}
\title{Federated ComDim}
\usage{
federateComDim(loginFD,
                      logins,
                      func,
                      symbol,
                      ncomp = 2,
                      scale = FALSE,
                      scale.block = TRUE,
                      threshold = 1e-8,
                      chunk = 500L,
                      mc.cores = 1,
                      width.cutoff = 500L
                      )
}
\arguments{
\item{loginFD}{Login information of the FD server}

\item{logins}{Login information of data repositories}

\item{func}{Encoded definition of a function for preparation of raw data
matrices. Two arguments are required: conns (list of Opal connections), 
symbol (name of the R symbol) (see datashield.assign).}

\item{symbol}{Encoded vector of names of the R symbols to assign in the
DataSHIELD R session on each server in \code{logins}.
The assigned R variables will be used as the input raw data.
Other assigned R variables in \code{func} are ignored.}

\item{ncomp}{Number of common dimensions}

\item{scale}{A logical value indicating if variables are scaled to unit
variance. Default, FALSE. See \code{MBAnalysis::ComDim}.}

\item{scale.block}{A logical value indicating if each block of variables is
divided by the square root of its inertia. Default, TRUE.
See \code{MBAnalysis::ComDim}.}

\item{threshold}{if the difference of fit<threshold then break the iterative
loop (default 1e-8)}

\item{chunk}{Size of chunks into what the resulting matrix is partitioned.
Default, 500L.}

\item{mc.cores}{Number of cores for parallel computing. Default, 1.}

\item{TOL}{Tolerance of 0. Default, \code{.Machine$double.eps}.}

\item{width.cutoff}{Default, 500L. See \code{deparse1}.}
}
\value{
A \code{ComDim} object. See \code{MBAnalysis::ComDim}.
}
\description{
Function for ComDim federated analysis on the virtual cohort
combining multiple cohorts: finding common dimensions in multiblock data.
}
