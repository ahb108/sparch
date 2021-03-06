% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.R
\name{spnetBuild}
\alias{spnetBuild}
\title{Build a geographical network from a network edge-list.}
\usage{
spnetBuild(sites, siteid, edges, edgew, cutoff = 0, verbose = TRUE)
}
\arguments{
\item{sites}{To Add.}

\item{siteid}{To Add.}

\item{edges}{To Add.}

\item{edgew}{To Add.}

\item{cutoff}{To Add.}

\item{verbose}{A logical variable indicating whether extra information on progress should be reported. Default is TRUE.}
}
\value{
An object of class SpatialLines.
}
\description{
Build a geographical network from a network edge-list.
}
\details{
This function takes a point dataset that can be mapped in geogrpahical space (e.g. a SpatialPointsDataFrame) and a basic set of from-and-to edges (optionally with weights) that might used in a standard network model, and it builds line data that can be mapped corrctly as a geographical layout.
}
