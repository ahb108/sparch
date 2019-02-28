
#' @export
removeIslands <- function(x){
  # Remove holes form polygon (based on function in wild1 library)
  require(sp)
  isles <-  unlist(lapply(x@Polygons, function(p) p@hole))
  p <- Polygons(x@Polygons[!isles], ID=x@ID)
  return(p)
}

##

#' @export
spdfClip <- function(toclip, clipby){

  #Class checks
  if(!class(toclip)[1] %in% c("SpatialLines","SpatialPolygons","SpatialLinesDataFrame","SpatialPolygonsDataFrame")){stop("'toclip' is not a SpatialLines* or SpatialPolygons* objects.")}
  if(!inherits(clipby, "SpatialPolygons")){stop("'clipby' is not a SpatialPolygons* object.")}
  #Projection compatibility
  if (!identical(proj4string(toclip),proj4string(clipby))){stop("'toclip' and 'clipby' are not the same CRS.")}
  #Set up up dataframe to fill if necessary
  if (class(toclip)[1]=="SpatialLinesDataFrame"){
    tc <- as(toclip,"SpatialLines")
    df <- toclip@data[1, ,drop=FALSE]
    df <- df[-1, ,drop=FALSE]
    df <- cbind(df,data.frame(RN=numeric(0)))
  } else if (class(toclip)[1]=="SpatialPolygonsDataFrame"){
    tc <- as(toclip,"SpatialPolygons")
    #df <- toclip@data[1,]
    df <- toclip@data[1, ,drop=FALSE]
    df <- df[-1, ,drop=FALSE]
    df <- cbind(df,data.frame(RN=numeric(0)))
  } else {
    tc <- toclip
  }
  #Loop through and clip (preserving dataframe links if necessary)
  clippedlist <- vector("list", length(tc))
  dflist <- vector("list", length(tc))
  for (i in 1:length(tc)){
    clippedtmp <- gIntersection(tc[i,],clipby)
    if (!is.null(clippedtmp)){
      if (class(tc)[1]=="SpatialLines"){
        RN <- lapply(slot(tc[i,],"lines"), function(x) {slot(x, "ID")})[[1]][1]
      } else {
        RN <- lapply(slot(tc[i,],"polygons"), function(x) {slot(x, "ID")})[[1]][1]
      }
      clippedtmp <- spChFIDs(clippedtmp, RN)
      clippedlist[[i]] <- clippedtmp      
      if (class(toclip)[1]=="SpatialLinesDataFrame"||class(toclip)[1]=="SpatialPolygonsDataFrame"){
        dfi <- toclip@data[i, ,drop=FALSE]
        dfi$RN <- RN
        dflist[[i]] <- dfi
      }
    }
  }
  clippedlist <- clippedlist[!sapply(clippedlist, is.null)]
  clippedall <- do.call("rbind", clippedlist)
  dflist <- dflist[!sapply(dflist, is.null)]
  df <- do.call("rbind", dflist) 
  #Recombine the results with the dataframe if necessary
  if (!is.null(clippedall)){
    if (class(toclip)[1]=="SpatialLinesDataFrame"){
      rownames(df) <- df$RN
      clippedall <- SpatialLinesDataFrame(clippedall, data=df)
    } else if (class(toclip)[1]=="SpatialPolygonsDataFrame"){
      rownames(df) <- df$RN
      clippedall <- SpatialPolygonsDataFrame(clippedall, data=df)
    }
  }
  return(clippedall)
}

##

#' @export
gExplode <- function(x) {
  # a modified version of something by Josh O'Brien
    a <- deparse(substitute(x))
    if(class(x)[[1]] == "SpatialLines" | class(x)[[1]] == "SpatialLinesDataFrame"){
        b <- unlist(lapply(x@lines, function(c) c@Lines))
        SpatialLines(lapply(seq_along(b), function(c) Lines(b[c], ID=paste0(a,c))), proj4string=CRS(proj4string(x)))
    } else if (class(x)[[1]] == "SpatialPolygons" | class(x)[[1]] == "SpatialPolygonsDataFrame"){
        b <- unlist(lapply(x@polygons, function(c) c@Polygons))
        SpatialPolygons(lapply(seq_along(b), function(c) Polygons(b[c], ID=paste0(a,c))), proj4string=CRS(proj4string(x)))
    } else {
        stop("Input must be of class SpatialLines* or SpatialPolygons*")
    }
}

##

#' @export
buffers <- function(x, bands, xIds=NULL, rings=TRUE, bMerge=FALSE, ...){
    for(i in 1:length(bands)){
        if (bMerge){
            pids <- paste("",bands[i], sep="_")
        } else if (!is.null(xIds)){
            pids <- paste(xIds,bands[i],sep="_")
        } else {
            pids <- paste(1:length(x),bands[i],sep="_")
        }
        b <- gBuffer(x, byid=TRUE,id=NULL, width=bands[i], ...)
        if (bMerge){ b <- gUnaryUnion(b) }
        if (i == 1){
            b <- spChFIDs(b, pids)
            bufs <- b
            a <- b
        } else {
            d <- b
            # clip out previous buffer
            if (rings){
                d <- gDifference(d, a, byid=TRUE, id=NULL)
                if (!bMerge){
                    d <- d[seq(1,length(d),sqrt(length(d))+1)]
                }
            }
            d <- spChFIDs(d, pids)
            bufs <- spRbind(bufs,d)
            a <- b
        }    
    }
    if (!is.null(xIds)){
        # Convert to SpatialPolygonsDataFrame
        pids <- sapply(slot(bufs, "polygons"), function(x) slot(x, "ID"))
        df <- data.frame(PID=pids, row.names=pids)
        bufs <- SpatialPolygonsDataFrame(bufs, data=df)
        if (!bMerge){
            bufs$xID <- sapply(strsplit(as.character(bufs$PID),"_"), "[", 1)
        }
        bufs$Buffer <- as.numeric(sapply(strsplit(as.character(bufs$PID),"_"), "[", 2))
    }
    return(bufs)
}

##

#' Simple and weighted Voronoi tesselation
#'
#' Function to produce a Voronoi tesselation (aka Dirichelet set or Thiessen polygons) from a set of point, without or without weights.
#'
#' @param x an object of class SpatialPoints*
#' @param ids identifiers for each point as lables for the resulting polygons
#' @param weights optional weight values for the calculation (e.g. measure of relative `influence' or `size' of each point)
#' @param win an object of class SpatialPolygons* for the exterior of the study area.
#' @param maxdist A cut-off (the radius of a circle) for the extent of a point's tesselation.
#' @param nsteps How many increments to use for defining the cruvature of a circle.
#' @return An object of class SpatialPolygonsDataFrame
#' @examples
#' w <- SpatialPolygons(list( Polygons(list(Polygon(cbind(c(0,0,100,100,0), c(0,100,100,0,0)))), "1")), 1:1)
#' pts <- data.frame(X=runif(3)*100,Y=runif(3)*100,ptID=letters[1:3], Size=c(1,2,4))
#' coordinates(pts) <- ~X+Y
#' vor <- voronoi(pts, ids=pts$ptID, win=w)
#' vorw <- voronoi(pts, ids=pts$ptID, weights=pts$Size, win=w)
#' plot(vor, lty="dotted", border="red")
#' plot(vorw, add=TRUE)
#' points(pts)
#' @import deldir
#' @export
voronoi <- function(x, ids, weights=NULL, win=NULL, maxdist=NULL, nsteps=36, ...){

    if (is.null(weights)){
        if (is.null(win)){ rw1 <- NULL } else { rw1 <- as.vector(t(bbox(win))) }
        a <- deldir(coordinates(x)[,1], coordinates(x)[,2], rw=rw1, ...)
        w <- tile.list(a)
        polys <- vector(mode='list', length=length(w))
        for (i in seq(along=polys)) {
            pcrds <- cbind(w[[i]]$x, w[[i]]$y)
            pcrds <- rbind(pcrds, pcrds[1,])
            polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
        }
        sp1 <- SpatialPolygons(polys)
        spdf1 <- SpatialPolygonsDataFrame(sp1, data=data.frame(PtID=ids, row.names=sapply(slot(sp1, 'polygons'), function(x) slot(x, 'ID'))))
        if (!is.null(win)){
            proj4string(spdf1) <- CRS(proj4string(win))
            res <- spdfClip(spdf1, win)
        }
        return(res)
    } else {
        for (a in 1:nrow(x)){
            da <- NULL
            for (b in 1:nrow(x)){
                if (a != b){
                    if (weights[a] <= weights[b]){
                        p1 <- coordinates(x)[a,]
                        p2 <- coordinates(x)[b,]
                        w1 <- weights[a]
                        w2 <- weights[b]
                    } else {
                        p1 <- coordinates(x)[b,]
                        p2 <- coordinates(x)[a,]
                        w1 <- weights[b]
                        w2 <- weights[a]
                    }
                    d12 <- sqrt(sum((p1 - p2) ^ 2)) #euclidean distance
                    c1 <- (w2^2*p1-w1^2*p2) / (w2^2-w1^2)
                    r1 <- (w1*w2*d12) / (w1^2-w2^2)
                    if (weights[a] == weights[b]){
                        tmpts <- rbind(x[a,],x[b,])
                        ac <- voronoi(tmpts,c(w1,w2),win=win,nsteps=nsteps)
                        ac <- as(ac[1,],"SpatialPolygons")
                    } else {
                        ac <- spEllipse(c1[1],c1[2],r1, nsteps=nsteps)
                    }
                    if (weights[a] > weights[b]){
                        ac <- gDifference(win,ac)
                    }
                    if (is.null(da)){
                        da <- gIntersection(ac,win)
                    } else {
                        da <- gIntersection(ac,da)
                    }
                }
            }
            if (!is.null(maxdist)){
                if (!is.numeric(maxdist) | maxdist <= 0){
                    stop("Maximum distance must be a postive number.")
                } else {
                    b <- gBuffer(x[a,], width=maxdist, quadsegs=round(nsteps/4,0))
                    da <- gIntersection(da,b)
                }
            }
            da <- spChFIDs(da, paste(ids[a],sep=""))
            da <- SpatialPolygonsDataFrame(da, data=data.frame(id=ids[a], weight=weights[a], row.names=paste(ids[a],sep="")))
            if (a==1){
                res <- da
            } else {
                res <- spRbind(res,da)
            }
        }
        proj4string(res) <- CRS(proj4string(pts))
        return(res)
    }
}

##

#' @export
polyAngle <- function(x, degsteps=1, type="azi", verbose=TRUE){
    if (!type %in% c("azi","sp")){
        stop("type argument must be either azi or sp.")
    }
    rotv <- vector(mode="numeric", length=nrow(x))
    if (length(x)>1 & verbose){
        print("Calculating dominant polygon angles...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(x), style=3)
    }
    for (b in 1:length(x)){
        if (length(x)>1 & verbose){ setTxtProgressBar(pb, b) }
        obj <- x[b,]
        rots <- seq(0,90,degsteps)
        ht <- abs(bbox(obj)[2,1] - bbox(obj)[2,2])
        wd <- abs(bbox(obj)[1,1] - bbox(obj)[1,2])
        area <- ht*wd
        minarea <- area
        for (a in 1:length(rots)){
            allcheck <- elide(obj, rotate=rots[a], center=apply(bbox(obj), 1, mean))
            thisht <- abs(bbox(allcheck)[2,1] - bbox(allcheck)[2,2])
            thiswd <- abs(bbox(allcheck)[1,1] - bbox(allcheck)[1,2])
            thisarea <- thisht*thiswd
            if (thisarea < minarea){
                minarea <- thisarea
                rot <- rots[a]
            }
        }
        obj1 <- elide(obj, rotate=rot, center=apply(bbox(obj), 1, mean))
        ranges <- apply(bbox(obj1),1,diff)
        if (ranges["x"]>ranges["y"]){
            if(type=="sp"){
                xy <- cbind(c(bbox(obj1)["x",1],bbox(obj1)["x",2]),c(bbox(obj1)["y",1]+(ranges["y"]/2),bbox(obj1)["y",1]+(ranges["y"]/2)))
            }
            rotv[b] <- 90-rot
        } else {
            if(type=="sp"){
                xy <- cbind(c(bbox(obj1)["x",1]+(ranges["x"]/2),bbox(obj1)["x",1]+(ranges["x"]/2)),c(bbox(obj1)["y",1],bbox(obj1)["y",2]))
            }
            rotv[b] <- 180-rot
        }
        if(type=="sp"){
            L <- Line(xy)
            lid <- sapply(slot(obj1, "polygons"), function(k) slot(k, "ID"))
            newSL <- SpatialLines(list(Lines(list(L),lid)))
            newSL <- SpatialLinesDataFrame(newSL,data.frame(LineID=lid, row.names=lid))
            res <- elide(newSL, rotate=-rot, center=apply(bbox(obj), 1, mean))
            proj4string(res) <- proj4string(obj)
            res <- gIntersection(res,obj)
            res <- spChFIDs(res, lid)
            if (b == 1){
                allres <- res
            } else {
                allres <- spRbind(allres,res)
            }
        }
    }
    if (length(x)>1 & verbose){ close(pb) }
    if (type=="azi"){ allres <- rotv }
    if (verbose){ print("Done.") }
    return(allres)
}
