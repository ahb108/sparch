#' Produce a very simple map scale-bar
#'
#' This function draw a simple scale-bar in the form of a box with a certain width and height (in map units) and hung on its top right corner.
#'
#' @param x Top left x coordinate
#' @param y Top left y coordinate
#' @return a plotted polygon
#' @export
scaleBar <- function(x, y, width, height,...){
  polygon(x=c(x, (x+width), (x+width), x), y=c(y, y, (y-height), (y-height)), ...)
}

#' Produce an circle or ellipse as a SpatialPolygons object
#'
#' This function draw a simple circle or ellipse, centred on a certain coordinate pair and of a certain size (the circle's radius or height, width and azimuth for an ellipse). The smooothness of the circle can be controlled by the chosen number of vertices placed along the perimeter, and the circle can be assigned a projection system if desired. 
#'
#' @param xc x coordinate of circle or ellipse centre
#' @param yc y coordinate of circle or ellipse centre
#' @param a the size of the circle radius in map units or the length of the longer axis of the ellipse
#' @param b the size of the shorter axis of the ellipse, if necessary (default is circular with b=a).
#' @param angle the map orientation of the shorter axis of the ellipse, if necessary (default is longer axis is oriented north).
#' @param nsteps the number of vertices used to describe the circle of the ellipse, if necessary (i.e. the smoothness of the plotted circle, default is 360)
#' 
#' @return a plotted polygon
#' @export
spEllipse <- function(xc, yc, a, b=a, angle=0, nsteps=360, proj4string=CRS(as.character(NA)), noSp=FALSE, ...){

  theta <- seq(0,(2*pi),len=nsteps)
  angle <- (90 - angle) * (pi /180)

  x <- xc + ((a * cos(theta) * cos(angle)) - (b * sin(theta) * sin(angle)))
  y <- yc + ((a * cos(theta) * sin(angle)) + (b * sin(theta) * cos(angle)))

  cpol<-cbind(x,y)
  cpol <- rbind(cpol, cpol[1,])
  if (noSp){
    return(cpol)
  } else {
    cpolp <- Polygons(list(Polygon(cpol)), ID="1")
    cpolsp <- SpatialPolygons(list(cpolp), proj4string=proj4string)
    return(cpolsp)
  }
}


#' Produce a simple north arrow
#'
#' This function draw a simple north arrow of a certain size and at a certain map location.
#'
#' @param x x coordinate for the centre of the bounding circle of the north arrow
#' @param y y coordinate for the centre of the bounding circle of the north arrow
#' @param r Radius of the bounding circle of the north arrow
#' 
#' @return a plotted polygon
#' @export

northArrow <- function(x, y, r, type="simple", fill="black", bkgrd=NA, ...){
  
  if (type=="simple"){
    polygon(spEllipse(x,y,r, noSp=TRUE), col=bkgrd, ...)
    bl <- c((x + r * sin(210*pi/180)),(y + r * cos(210*pi/180)))
    br <- c((x + r * sin(150*pi/180)),(y + r * cos(150*pi/180)))
    polygon(x=c(x, br[1], x,  bl[1],x),y=c(y+r, br[2], (y-(0.5 * r)), bl[2],y+r), col=fill, ...)
  } else {stop("This north arrow type does not exist. Try leaving `type' blank.")}
}

##

#' @export
# Create a graduated symbol function
gradSymbols <- function(x, mincex, maxcex, valrange=NULL, sqrttrans=FALSE,...){
  if (sqrttrans) { x <- sqrt(x) }
  if (is.null(valrange)){
    minx <- min(x)
    maxx <- max(x)
  } else if (is.vector(valrange) & length(valrange)==2){
    minx <- valrange[1]
    maxx <- valrange[2]
  } else { stop("valrange must be NULL or a vector of length 2.\n") }
  y <- mincex + (((x-minx)*maxcex - (x-minx)*mincex) / (maxx-minx))
  return(y)
}


#' Addition of both a north arrow and scalebar to map
#'
#' This is a helper function to add both north arrow and scalebar to a map.
#'
#' @param x x coordinate for the centre of the added elements (the centre of the text label)
#' @param y y coordinate for the centre of the added elements (the centre of the text label)
#' @param scalesize width of the scalebar in map units (this function assumes projcted units but will not throw a warning if unprojected units are used.
#' @param scalelabel text label for the state width of the scalebar (it requires the user to enter this and is not automatically generated from scalesize)
#' 
#' @return a plotted north arrow and scale-bar
#' @export
#' 
mapElements <- function(x, y, scalesize, scalelabel, cex.label=0.75, font.label=1, col.label="black", type.north="simple", lwd.north=0.5, offset.north=height.scale*1.5, r.north=scalesize/5, bkgrd.north="white", fill.north="black", offset.scale=height.scale, height.scale=scalesize/8, border.scale="black", col.scale="white", lwd.scale=0.5){
    
    northArrow(x, y+(r.north+offset.north), r=r.north, type=type.north, lwd=lwd.north, fill=fill.north, bkgrd=bkgrd.north)
    scaleBar(x-(scalesize/2), y-(offset.scale), width=scalesize, height=height.scale, lwd=lwd.scale, col=col.scale, border=border.scale)
    text(x, y, labels=scalelabel, cex=cex.label, font=font.label, col=col.label)
}


#' Produce pie charts at specified locations on a map
#'
#' This function draw a pie chart at specified locations on a map. 
#'
#' @param x x coordinate of where the pie will be plotted
#' @param y y coordinate of where the pie will be plotted
#' @param values the values that will become the percentages of the pie.
#' @param sizes the size of each pie.
#' @param edges the number of increments used to describe the circle of the pie.
#' @param clockwise whether the divisions are in clockwise order.
#' @param init.angle where the pie starts te divisions
#' @param col pie colors (length must be the same as ncol values)
#' @param border border color(s) of the pies
#' @param ... passed to default plot function
#' 
#' @examples
#' data(SainteMarie85)
#' grains$SymSizes <- grains$SiloCount*8 #specify symbol sizes in map units
#' plotmat <- arrangeSymbols(grains, grains$SymSizes) #alternative pie locations
#' flines <- flylines(coordinates(grains)[,1],coordinates(grains)[,2],plotmat[,1],plotmat[,2]) #fly-lines from real locations to pie locations.
#' plot(boghar, col="grey75", border=NA, axes=TRUE, main="Barley (brown) and wheat (yellow)\nin 19th grain silos (Boghar, Algeria)") #background map
#' lines(flines, lwd=0.5) #fly-lines
#' pieSymbols(x=plotmat[,1], y=plotmat[,2], values=grains@data[,c("BarleyHL","WheatHL")], sizes=grains$SymSizes, col=c("saddlebrown","yellow")) #pies
#' points(grains, pch=19, cex=0.3, col="black") #real locations
#' @export
#' 
pieSymbols <- function (x, y, values, sizes=NULL, edges=360, clockwise=FALSE, init.angle=if(clockwise) 90 else 0, density=NULL, angle=45, col=NULL, border=NULL, lty=NULL, main=NULL, ...){
    ## This is a modification of the basic pie() function to allow overplotting on map data
    t2xy <- function(t){
        t2p <- twopi * t + init.angle * pi/180
        list(x=r * cos(t2p), y=r * sin(t2p))
    }
    for (a in 1:nrow(values)){
        v <- as.matrix(values)[a,]
        if (is.null(sizes)){
            sizes <- rep((max(y)-min(y))/3,length(sizes))
        } else if (length(sizes)==1){
            r <- rep(sizes,length(sizes))
        } else {
            r <- sizes[a]
        }
        if (!is.numeric(v) | any(is.na(v) | v < 0)){ 
            stop("values must be positive and numeric.")
        }
        v <- c(0, cumsum(v)/sum(v))
        dv <- diff(v)
        nv <- length(dv)
        if (!is.null(border)){ border <- rep_len(border, nv) }
        if (!is.null(lty)){ 
            lty <- rep_len(lty, nv)
            angle <- rep(angle, nv)
        }
        if (!is.null(density)){ density <- rep_len(density, nv) }
        if (clockwise){ twopi <- -2 * pi } else { twopi <- 2 * pi }
        for (i in 1L:nv) {
            n <- max(2, floor(edges * dv[i]))
            P <- t2xy(seq.int(v[i], v[i + 1], length.out=n))
            polygon(c(x[a]+P$x,x[a]), c(y[a]+P$y,y[a]), density=density[i], angle=angle[i], 
                    border=border[i], col=col[i], lty=lty[i])
        }
    }
}

#' Rearrange the location of those symbols that overalap with each other on a map.
#'
#' This function takes the orginal plot locations and suggests new locations that avoid each other, given a stated plotting size for each symbol. 
#'
#' @param x a SpatialPoints* object
#' @param plotsizes a vector of sizes (typically radii or perhaps half a bounding box) for the plotted symbols 
#' 
#' @return matrix of new suggested x and y locations for the centre of each plot symbol.
#' @export
#' 
arrangeSymbols <- function(x, plotsizes, method="buffers1"){
    sitebuffs <- gBuffer(x, width=plotsizes, byid=TRUE)
    polyids <- sapply(slot(sitebuffs, "polygons"), function(x) slot(x, "ID"))
    names(plotsizes) <- polyids
    nogozone <- gBuffer(x, width=max(plotsizes))    
    overlapcheck <- rowSums(gOverlaps(sitebuffs, byid=TRUE)) > 0
    overlaps <- x[overlapcheck,]
    nonoverlaps <- x[!overlapcheck,]
    nooverlapsbuff <- sitebuffs[!overlapcheck,]
    overlapssizes <- plotsizes[overlapcheck]
    overlaps <- overlaps[order(overlapssizes),]
    overlapssizes <- overlapssizes[order(overlapssizes)]
    res <- nooverlapsbuff
    sizeorder <- names(sort(overlapssizes, decreasing=TRUE))
    for (a in 1:length(sizeorder)){
        buff <- gBuffer(overlaps[sizeorder[a],], width=overlapssizes[sizeorder[a]], byid=TRUE)
        nonoverlapcheck <- is.null(gIntersection(nogozone,buff))
        if (nonoverlapcheck){
            buff <- spChFIDs(buff,sizeorder[a])
            res <- spRbind(res,buff)
            nogozone <- spRbind(nogozone,buff)
        } else {
            tmp1 <- gBuffer(nogozone, width=overlapssizes[sizeorder[a]])
            tmp2 <- gBuffer(tmp1, width=overlapssizes[sizeorder[a]])
            newareas <- gDifference(tmp2,tmp1)
            set.seed(123)
            pts <- spsample(newareas,1000, type="random")
            nearest <- which.min(gDistance(overlaps[sizeorder[a],], pts, byid=TRUE))
            newpt <- pts[nearest,]
            newbuff <- gBuffer(newpt, width=overlapssizes[sizeorder[a]])
            newbuff <- spChFIDs(newbuff,sizeorder[a])
            res <- spRbind(res,newbuff)
            nogozone <- spRbind(nogozone,newbuff)
        }
    }
    res <- coordinates(res)
    res <- as.data.frame(res[order(as.numeric(row.names(res))), ], stringsAsFactors=FALSE)
    return(res)
}

#' Produce flyout lines for where mapped pies and other charts need to avoid overlapping.
#'
#' This function takes the orginal plot locations and the suggested new plot locations created by arrangeSymbols and creates lines between them that can be plotted. 
#'
#' @param x a vector of x coordinates for the original locations
#' @param y a vector of y coordinates for the original locations
#' @param xplot a vector of x coordinates for where the charts will actually be plotted.
#' @param yplot a vector of y coordinates for where the charts will actually be plotted.
#' 
#' @return SpatialLines
#' @export
#' 
flylines <- function(x, y, xplot, yplot){
    fl <- NULL
    for (a in 1:length(x)){
        xmid <- mean(c(x[a],xplot[a]))
        ymid <- mean(c(y[a],yplot[a]))
        tmp <- cbind(c(x[a], xmid, xplot[a]),c(y[a], ymid, yplot[a]))
        newSL <- SpatialLines(list(Lines(list(Line(tmp)),a)))
        newSL <- spChFIDs(newSL, as.character(a))
        if (a!=1 & !is.null(fl)){
            fl <- spRbind(fl,newSL)
        } else {
            fl <- newSL
        }
    }
    return(fl)
}

#' Add a custom colour ramp to a plot
#'
#' This function gives greater flexibility for the plotting of colour ramps for maps and other plots.
#'
#' @export
#' 
colRibbon <- function(lut, breaks, xlim, ylim, nticks=3, ticksize=0.3, labels=TRUE, at=NULL, vertical=TRUE, border.box=NULL, col.ticks="black", col.labs="black", lwd.box=1, lwd.ticks=1, lty.box="solid", lty.ticks="solid", offset.labs=0.5, cex.labs=0.7, font.labs=1, ...){
    ##inspired by https://gist.github.com/johncolby/993830#file-colorbar-r
    ys <- seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/length(lut))
    for (i in 1:(length(lut))) {
        rect(xlim[1],ys[i],xlim[2],ys[i+1], col=lut[i], border=NA)
    }
    if (labels[1]!=FALSE){
        yscale <- (ylim[2]-ylim[1])/(breaks[length(breaks)]-breaks[1])
        if (is.null(at)){
            labels <- seq(breaks[1], breaks[length(breaks)], len=nticks)
            at <- ((labels - labels[1])*yscale)+ylim[1]
        } else {
            if (class(at)!="numeric"){
                stop("at must be NULL or a numeric vector of the same length as 'labels'.")
            }
            if (!is.numeric(labels[1])){
                labels <- at
            }
            if (length(at) != length(labels)){
                stop("at must be NULL or a numeric vector of the same length as 'labels'.")
            }
            at <- ((at - at[1])*yscale)+ylim[1]
        }
        ticksize <- ticksize*(xlim[2]-xlim[1])
        if (!is.null(border.box)){
            rect(xlim[1], ylim[1], xlim[2], ylim[2], border=border.box, lwd=lwd.box, lty=lty.box)
        }
        segments(x0=rep(xlim[2],length(at)), x1=rep(xlim[2]+ticksize,length(at)), y0=at, col=col.labs)
        text(x=rep(xlim[2]+ticksize,length(at)), y=at, labels=labels, cex=cex.labs, pos=4, offset=offset.labs)
    }
}
