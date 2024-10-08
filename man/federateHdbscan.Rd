% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client_func.R
\name{federateHdbscan}
\alias{federateHdbscan}
\title{Federated hdbscan}
\usage{
federateHdbscan(loginFD,
                       logins,
                       func,
                       symbol,
                       metric = 'euclidean',
                       minPts = 10,
                       chunk = 500L,
                       mc.cores = 1,
                       TOL = .Machine$double.eps,
                       width.cutoff = 500L,
                       ...)
}
\arguments{
\item{loginFD}{Login information of the FD server}

\item{logins}{Login information of data repositories}

\item{func}{Encoded definition of a function for preparation of raw data
matrices. Two arguments are required: conns (list of DSConnection-classes), 
symbol (name of the R symbol) (see datashield.assign).}

\item{symbol}{Encoded vector of names of the R symbols to assign in the
DataSHIELD R session on each server in \code{logins}.
The assigned R variables will be used as the input raw data.
Other assigned R variables in \code{func} are ignored.}

\item{metric}{Either \code{euclidean} or \code{correlation} for distance
metric between samples. 
For Euclidean distance, the data from each cohort will be centered (not
scaled) for each variable.
For correlation-based distance, the data from each cohort will be centered
scaled for each sample.}

\item{minPts}{Minimum size of clusters, see \code{dbscan::hdbscan}.
Default, 10.}

\item{chunk}{Size of chunks into what the resulting matrix is partitioned.
Default, 500L.}

\item{mc.cores}{Number of cores for parallel computing. Default, 1.}

\item{TOL}{Tolerance of 0. Default, \code{.Machine$double.eps}.}

\item{width.cutoff}{Default, 500L. See \code{deparse1}.}

\item{...}{see \code{dbscan::hdbscan}}
}
\value{
An object of class \code{hdbscan}.
}
\description{
Function for hdbscan federated analysis on the virtual cohort
combining multiple cohorts.
}
