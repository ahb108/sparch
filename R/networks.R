#' @title Build a geographical network from a network edge-list.
#' @description Build a geographical network from a network edge-list.  
#' @param sites To Add.
#' @param siteid To Add.
#' @param edges To Add.
#' @param edgew To Add.
#' @param cutoff To Add.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details This function takes a point dataset that can be mapped in geogrpahical space (e.g. a SpatialPointsDataFrame) and a basic set of from-and-to edges (optionally with weights) that might used in a standard network model, and it builds line data that can be mapped corrctly as a geographical layout. 
#' 
#' @return An object of class SpatialLines.
#' @export
spnetBuild <- function(sites, siteid, edges, edgew, cutoff=0, verbose=TRUE){  
    if (verbose){
        print("Creating spatial lines...")
        flush.console()
        pb <- txtProgressBar(min=1, max=nrow(sites), style=3)
    }
    edges$Path <- paste("from_",as.character(edges$From),"_to_",as.character(edges$To),sep="")
    check <- TRUE
    for (i in 1:nrow(sites)){
        if (verbose){ setTxtProgressBar(pb,i) }
        xyi <- coordinates(sites[i,])
        tos <- sites[-i,]
        for (j in 1:nrow(tos)){
            lid <- paste("from_",as.character(sites@data[i,siteid]),"_to_",as.character(tos@data[j,siteid]),sep="")
            wt <- edges[edges$Path == lid,edgew]
            if (length(wt)>0){
                if (wt >= cutoff){
                    xyj <- coordinates(tos[j,])
                    L <- Line(rbind(t(matrix(xyi)),t(matrix(xyj))))
                    Ls <- Lines(list(L),j)
                    newSL0 <- SpatialLines(list(Ls))
                    newSL <- spChFIDs(newSL0, lid)
                    df <- data.frame(Weight=wt, row.names=lid)
                    newSL <- SpatialLinesDataFrame(newSL,df)
                    if (check){
                        sitenet <- newSL
                        check <- FALSE
                    } else {
                        sitenet <- spRbind(sitenet,newSL)
                    }
                }
            }
        }
    }
    if (verbose){ close(pb) }
    return(sitenet)
}

