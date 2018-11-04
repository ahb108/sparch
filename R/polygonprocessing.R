
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

#' @export
voronoi <- function(x,xID,win=NULL, ...){
    # wraps deldir for SpatialPolygons
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
    spdf1 <- SpatialPolygonsDataFrame(sp1, data=data.frame(PtID=xID, row.names=sapply(slot(sp1, 'polygons'), function(x) slot(x, 'ID'))))
    if (!is.null(win)){
        proj4string(spdf1) <- CRS(proj4string(win))
        spdf1 <- spdfClip(spdf1, win)
    }
    return(spdf1)
}

##

#' @export
mwvoronoi <- function(pts, ptID, ptWeights=rep(1,nrow(pts)), win=NULL, maxdist=NULL, nsteps=36, ...){

    # Multiplicatively weighted Voronoi tesselation
    for (a in 1:nrow(pts)){
        da <- NULL
        for (b in 1:nrow(pts)){
            if (a != b){
                if (ptWeights[a] <= ptWeights[b]){
                    p1 <- coordinates(pts)[a,]
                    p2 <- coordinates(pts)[b,]
                    w1 <- ptWeights[a]
                    w2 <- ptWeights[b]
                } else {
                    p1 <- coordinates(pts)[b,]
                    p2 <- coordinates(pts)[a,]
                    w1 <- ptWeights[b]
                    w2 <- ptWeights[a]
                }
                d12 <- sqrt(sum((p1 - p2) ^ 2)) #euclidean distance
                c1 <- (w2^2*p1-w1^2*p2) / (w2^2-w1^2)
                r1 <- (w1*w2*d12) / (w1^2-w2^2)
                if (ptWeights[a] == ptWeights[b]){
                    tmpts <- rbind(pts[a,],pts[b,])
                    ac <- voronoi(tmpts,c(w1,w2),win=win,nsteps=nsteps)
                    ac <- as(ac[1,],"SpatialPolygons")
                } else {
                    ac <- spEllipse(c1[1],c1[2],r1, nsteps=nsteps)
                }
                if (ptWeights[a] > ptWeights[b]){
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
                b <- gBuffer(pts[a,], width=maxdist, quadsegs=round(nsteps/4,0))
                da <- gIntersection(da,b)
            }
        }
        da <- spChFIDs(da, paste(pts$ptID[a],sep=""))
        da <- SpatialPolygonsDataFrame(da, data=data.frame(ptID=pts$ptID[a], weight=ptWeights[a], row.names=paste(pts$ptID[a],sep=""))) # convert to SpatialPolygonsDataFrame
        if (a==1){
            res <- da
        } else {
            res <- spRbind(res,da)
        }
    }
    proj4string(res) <- CRS(proj4string(pts))
    return(res)
}

##
