
#' @title Convert an alphahull to a SpatialPolygonsDataFrame
#'
#' @description Convert an ahull object from the alphahull package to a SpatialPolygonsDataFramefor wider use.
#' 
#' @param x An object of class ahull from the alphahull package
#' @param increment A single numeric value for the number of increments with which to describe the circumference of a circle.
#' @param rnd A single numeric value controlling the precision of the estimate of the arcs (the default is usually sensible).
#' @param proj4string Projection parameters if the alphahull is based on geographic locations.
#' 
#' @return An object of class SpatialPolygonsDataFrame
#' @examples
#' ##Example point set and its alphahull
#' n <- 300
#' theta <- runif(n, 0,2 * pi)
#' r <- sqrt(runif(n, 0.25^2, 0.5^2))
#' x <- cbind(0.5 + r*cos(theta), 0.5 + r*sin(theta))
#' ahullobj <- ahull(x, alpha=0.1)
#' ahobj <- ah2sp(ahullobj)
#' ## Plot
#' par(mfrow=c(1,2), pty="s")
#' plot(ahullobj, main="ahull object")
#' plot(ahobj, col="lightgrey", pbg="white", axes=TRUE, main="sp object")
#' points(x)
#' @import alphahull
#' @import maptools
#' @export
#' 
ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){

  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  xdf <- as.data.frame(x$arcs)     
  xdf <- subset(xdf,xdf$r > 0)
  res <- NULL  
  if (nrow(xdf) > 0){
    linesj <- list()
    prevx<-NULL
    prevy<-NULL
    j<-1
    for (i in 1:nrow(xdf)){
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      #Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)
      angles <- anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length=ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      if (is.null(prevx)){
        prevx <- x
        prevy <- y
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
          if (i == nrow(xdf)){
            prevx <- append(prevx, x[2:ipoints])
            prevy <- append(prevy, y[2:ipoints])
            prevx[length(prevx)] <- prevx[1]
            prevy[length(prevy)] <- prevy[1]
            coordsj <- cbind(prevx, prevy)
            colnames(coordsj) <- NULL
            linej <- Line(coordsj)
            linesj[[j]] <- Lines(linej, ID = as.character(j))
          } else {
            prevx <- append(prevx, x[2:ipoints])
            prevy <- append(prevy, y[2:ipoints])
          }
    } else {
      prevx[length(prevx)] <- prevx[1]
      prevy[length(prevy)] <- prevy[1]
      coordsj<-cbind(prevx, prevy)
      colnames(coordsj) <- NULL
      linej <- Line(coordsj)
      linesj[[j]] <- Lines(linej, ID=as.character(j))
      j <- j+1
      prevx <- NULL
      prevy <- NULL
    }
  }
  lspl <- SpatialLines(linesj)
  lns <- slot(lspl, "lines")
  polys <- sapply(lns, function(x) { 
    crds <- slot(slot(x, "Lines")[[1]], "coords")
    identical(crds[1, ], crds[nrow(crds), ])
  }) 
  polyssl <- lspl[polys]
  list_of_Lines <- slot(polyssl, "lines")
  sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID="1")), proj4string=proj4string)
  hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
  areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
  df <- data.frame(hid, areas)
  names(df) <- c("HID", "Area")
  rownames(df) <- df$HID
  res <- SpatialPolygonsDataFrame(sppolys, data=df)
  res <- res[which(res@data$Area > 0),]
}  
  return(res)
}


##
 
