#'
#' @title
#' Dual Mesh Weights
#'
#' @description
#' Computes the dual mesh weights for the integration of the point process
#'
#'
#' @param dmesh Dual mesh
#' @param region A region
#'
#'
#' @details
#' none
#'
#' @return
#' A list with element weights containing a vector of areas for each cell within the region area
#'
#' @references
#' none
#'
#'
#' @importFrom sf st_as_sf
#' @importFrom sf st_intersects
#' @importFrom sf st_intersection
#' @importFrom sf st_area
#'
#' @export
#'
#'
dmesh_weights <- function(dmesh, region){

  #--------------------------------------------------------------
  ### Find the intersection between the polygons in the dual mesh
  ### and the location domain
  #--------------------------------------------------------------

  # An intersect and a within are used to first find polygons that will
  # be cut to reduce the duration of the intersection (about twice as fast)
  if(nrow(region)>1){
    region<-st_union(region)
  }

  dm<-dmesh$dmesh

  ### Calculate weight
  #dm$id <- 1:nrow(dm)
  weights <- numeric(nrow(dm))
  overlaps <- lengths(st_intersects(dm, region))
  #within <- lengths(st_within(dm, region))
  within <- as.integer((1:nrow(dm))%in%(st_contains(region, dm)[[1]])) # much faster than the previous line https://github.com/r-spatial/sf/issues/1261
  o <- overlaps > 0L & !within > 0L
  suppressWarnings(
    cuts <- st_intersection(dm[o, ],region)[,names(dm)]
  )
  st_geometry(cuts)<-"geometry"
  dmeshcuts <- rbind(cuts, dm[within > 0L, ])
  #dmeshcuts <- dmeshcuts[order(dmeshcuts$id), ]
  #w <- which(as.logical(overlaps))
  a <- as.numeric(st_area(dmeshcuts))
  areas <- data.frame(id=dmeshcuts$id,area=a)
  areas <- aggregate(area~id,data=areas,FUN=sum)
  #if(length(w)!=length(areas)){
  #  stop("Number of resulting polygons different from the number of touching cells")
  #}
  weights[areas$id] <- areas$area

  ### Check to make sure there are integration points with 0 weights
  if(all(weights > 0)){
    stop("There needs to be some weights that are of 0")
  }
  dmesh[["weights"]]<-weights
  dmesh[["dmeshcuts"]]<-dmeshcuts
  dmesh

}
