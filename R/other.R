#' @title Handle the positioning and labelling of a (spatstat) colour ramp more easily 
#' @description Handle the positioning and labelling of a (spatstat) colour ramp more easily.  
#' @param x x-axis location of the colourramp lower-left corner
#' @param x y-axis location of the colourramp lower-left corner.
#' @param xprop proportion of plot width to be allocated to the colour ramp width
#' @param yprop proportion of plot height to be allocated to the colour ramp height
#' @param zlim minimum and maximum values for the colour ramp
#' @param cm object of class 'colourmap' to plot (or default spatstat colour ramp if left null)
#' @param steps how many numerical steps to label in the ramp legend
#' @param sigdigits how to round the digits in the labels
#' @param cex.axis size of the labels
#' @param las orientation of the labels
#' @param vertical whether the ramp is places as a vertical or horizontal ribbon
#' @return a image ribbon added to an existing plot.
#' @examples
#' exampledensitymap <- density(cells, 0.05)
#' plot(exampledensitymap, ribbon=FALSE)
#' ribbonplot(x=0.1, y=0.1, xprop=0.03, yprop=0.2, zlim=c(min(exampledensitymap),max(exampledensitymap)), cex.axis=0.6, col.axis="white", col.ticks="white", las=2)
#' @import stats
#' @export
ribbonplot <- function(x, y, xprop, yprop, zlim, cm=NULL, steps=4, rounding=3, vertical=TRUE,...){
    plotlim <- par("usr")
    if (is.null(cm)){
        cm <- colourmap(Kovesi$values[[29]],range=c(zlim[1],zlim[2]))
    }
    xlimramp <- c(x, x+((plotlim[2]-plotlim[1])*xprop))
    ylimramp <- c(y, y+((plotlim[4]-plotlim[3])*yprop))
    xticks <- seq(zlim[1], zlim[2], zlim[2]/steps)
    xticks.text <- round(xticks, rounding)
    plot(cm, vertical=vertical, main="", ylim=ylimramp, xlim=xlimramp, add=TRUE, at=reScale(xticks,to=c(ylimramp[1],ylimramp[2])), labels=xticks.text, ...)
}

#' @title Rescale a numeric vector to a specified minimum and maximum 
#' @description Rescale a numeric vector to a specified minimum and maximum.  
#' @param x numeric vector to smooth.
#' @param type what kind of rescaling to perform. Current options are 'simple' (default) and 'normal' which produces a z-score and 'custom' for which the 'to' argument must be specified.
#' @param to numeric vector of length 2 specifying the minimum and maximum value to perform a linear rescale between (default is 0 and 1)
#' @param na.rm Set to TRUE,this removes NAs before rescaling.
#' @return A numeric vector of rescaled values.
#' @examples
#' reScale(15:200)
#' @import stats
#' @export
reScale <- function(x, type="simple", to=c(0,1), na.rm=TRUE){

    types <- c("simple","normal")
    if (!type %in% types){
        stop("The rescale type you have chosen is not currently an option.")
    }
    if (max(x)-min(x)==0){
        warning("All the values in x are the same, and will just be recentred on 0 if type='normal' or max(to) if type='simple'.")
        if (type=="normal"){ res <- rep(0,length(x)) } else { res <- rep(max(to), length(x)) }
        return(res)
    }
    if (na.rm){ x <- na.omit(x) }
    if (type=="normal"){
        res <- (x-mean(x))/sd(x)
    } else {
        xrange <- range(x)
        mfac <- (to[2] - to[1])/(xrange[2] - xrange[1])
        res <- to[1] + (x - xrange[1]) * mfac
    }
    return(res)
}

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

#' @export
greyscales <- function(x, n, start=0, end=1, gamma=1, alpha=1, setrange=NULL){
    tmpcols <- gray.colors(n=n, start=start, end=end, gamma=gamma, alpha=alpha)
    if (!is.null(setrange)){
        ticks <- c(setrange[2],0,setrange[1])
        if (length(setrange)!=2 | any(!is.numeric(setrange))){
            stop("setrange must be a numeric vector of length 2.")
        } else {
            if (setrange[2] <= setrange[1]){
                stop("setrange must be a numeric vector of length 2 in ascending order.")
            }
            if (minValue(x) < setrange[1] & maxValue(x) > setrange[2]){
                tmpn <- n-2
                mybrks <- c(minValue(x), seq(setrange[1],setrange[2],(setrange[2]-setrange[1])/tmpn), maxValue(x))
                mycolrs <- c(tmpcols[1], gray.colors(tmpn, start=start, end=end, gamma=gamma, alpha=alpha), tmpcols[length(tmpcols)])
            } else if (minValue(x) < setrange[1]){
                tmpn <- n-1
                mybrks <- c(minValue(x), seq(setrange[1],setrange[2],(setrange[2]-setrange[1])/tmpn))
                mycolrs <- c(tmpcols[1], gray.colors(tmpn, start=start, end=end, gamma=gamma, alpha=alpha))
            } else if (maxValue(x) > setrange[2]){
                tmpn <- n-1
                mybrks <- c(seq(setrange[1],setrange[2],(setrange[2]-setrange[1])/tmpn), maxValue(x))
                mycolrs <- c(gray.colors(tmpn, start=start, end=end, gamma=gamma, alpha=alpha), tmpcols[length(tmpcols)])
            } else {
                tmpn <- n
                mybrks <- seq(setrange[1],setrange[2],(setrange[2]-setrange[1])/tmpn)
                mycolrs <- gray.colors(n, start=start, end=end, gamma=gamma, alpha=alpha)
            }
        }
    } else {
        ticks <- c(maxValue(x),0,minValue(x))
        mybrks <- seq(minValue(x),maxValue(x),(maxValue(x)-minValue(x))/n)
        mycolrs <- gray.colors(n, start=start, end=end, gamma=gamma, alpha=alpha)
    }
    return(list(breaks=mybrks, cols=mycolrs, ticks=ticks))
}
