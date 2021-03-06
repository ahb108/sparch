% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertformats.R
\name{ngrResolve}
\alias{ngrResolve}
\title{Resolve a Britsh or Irish National Grid Reference.}
\usage{
ngrResolve(ngr, centre = FALSE, grid = "British")
}
\arguments{
\item{ngr}{A character string with a British National Grid Reference}

\item{centre}{Return the coordinate of the origin of the local square (default) or the centre.}

\item{grid}{Input grid system (currently British or Irish)}
}
\value{
A numeric coordinate pair
}
\description{
Convert a Britsh or Irish National Grid Reference to a coordinate pair.
}
\examples{
ngrResolve("SJ456782")
}
