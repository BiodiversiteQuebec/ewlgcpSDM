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
#'
#'
#' @return
#'
#'
#' @references
#'
#'
#' @examples
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

  dm<-dmesh$dmesh

  ### Calculate weight
  dm$id <- 1:nrow(dm)
  weights <- numeric(nrow(dm))
  overlaps <- lengths(st_intersects(dm, region))
  within <- lengths(st_within(dm, region))
  o <- overlaps > 0L & !within > 0L
  suppressWarnings(
    cuts <- st_intersection(region, dm[o, ])
  )
  st_geometry(cuts)<-"geometry"
  dmeshcuts <- rbind(cuts, dm[within > 0L, ])
  dmeshcuts <- dmeshcuts[order(dmeshcuts$id), ]
  w <- which(as.logical(overlaps))
  areas <- as.numeric(st_area(dmeshcuts))
  if(length(w)!=length(areas)){
    stop("Number of resulting polygons different from the number of touching cells")
  }
  weights[w] <- areas

  ### Check to make sure there are integration points with 0 weights
  if(all(weights > 0)){
    stop("There needs to be some weights that are of 0")
  }
  dmesh[["weights"]]<-weights
  dmesh[["dmeshcuts"]]<-dmeshcuts
  dmesh

}
