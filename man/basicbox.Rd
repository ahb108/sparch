% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygonprocessing.R
\name{basicbox}
\alias{basicbox}
\title{Produce a basic square or rectangular mapping polygon}
\usage{
basicbox(ll, ur, ul = c(ll[1], ur[2]), lr = c(ur[1], ll[2]), proj4string = NA)
}
\arguments{
\item{ll}{coordinates of the lower left corner}

\item{ur}{coordinates of the upper right corner}

\item{ul}{coordinates of the upper left corner (by default set automatically)}

\item{lr}{coordinates of the lower right corner (by default set automatically)}

\item{proj4string}{projection string of class CRS.}
}
\value{
An object of class SpatialPolygonsDataFrame
}
\description{
Function to produce a basic square or rectangular box (e.g. for polygon clipping or plot borders)
}
\examples{
utm34n <- CRS("+init=epsg:32634") # Greek UTM (west)
b <- basicbox(c(681900,4010900), c(682500,4011400), proj4string=utm34n)
plot(b, axes=TRUE)
}
