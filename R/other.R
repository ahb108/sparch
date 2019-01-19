
#' @export
isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  #Checks whether a value is a whole number (used in nnhistMC reporting)
  abs(x - round(x)) < tol
}

#' @export
# Check for odd number
isOdd <- function(x){ x %% 2 != 0 }

## 

#' @export
# Find-Replace function
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
##

#' @export
spJitter <- function(pts, xamount, yamount=xamount){
    proj <- NA
    if (!is.na(proj4string(pts)) | proj4string(pts)!="NA"){
        proj <- proj4string(pts)
    }
    if (class(pts) == "SpatialPointsDataFrame"){
        df <- cbind(coordinates(pts),pts@data)
        df[,1] <- jitter(df[,1],amount=xamount)
        df[,2] <- jitter(df[,2],amount=yamount)
        coordinates(df) <- df[,1:2]
        proj4string(df) <- proj
    } else if (class(pts) == "SpatialPoints"){
        df <- coordinates(pts)
        df[,1] <- jitter(df[,1],amount=xamount)
        df[,2] <- jitter(df[,2],amount=yamount)
        df <- SpatialPoints(df, proj4string=CRS(proj))
    } else {
        stop("Only works for SpatialPoints* at present.")
    }
    return(df)
}

#' @export
rybcolourmap <- function(range, ...) {
  col <- rybcolours(range, ...)
  z <- colourmap(col, range=range)
  return(z)
}

#' @export
rybcolours <- function(range, sealevel=0, ncolours=100, nbeach=0){
    ## modified from a routine by A. Baddeley
    stopifnot(is.numeric(range) && length(range)==2)
    stopifnot(all(is.finite(range)))
    yr <- colorRampPalette(c("yellow","orangered","darkred"), space="rgb")
    cb <- colorRampPalette(c("blue","cyan","yellow"), space="rgb")
    depths <- range[1]
    peaks <- range[2]
    dv <- diff(range)/(ncolours - 1)
    epsilon <- nbeach * dv/2
    lowtide <- max(sealevel - epsilon, depths)
    hightide <-  min(sealevel + epsilon, peaks)
    countbetween <- function(a, b, delta) { max(0, round((b-a)/delta)) }
    nsea <- countbetween(depths, lowtide, dv)
    nbeach <- countbetween(lowtide,  hightide, dv)
    nland <- countbetween(hightide,  peaks, dv)
    colours <- character(0)
    if(nsea > 0)  colours <- cb(nsea) # cyan/blue
    if(nbeach > 0)  colours <- c(colours,rep("yellow",nbeach)) # yellow
    if(nland > 0)  colours <- c(colours, yr(nland)) # darkred/yellow
    return(colours)
}
