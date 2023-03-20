#'
#' @title
#' Summarize effort to dual mesh
#'
#' @description
#' This function summarizes effort to dual mesh.
#'
#' @param dmesh A dual mesh sf object
#' @param obs An sf data.frame with observations of the target species
#' @param background An sf data.frame with observations of the target group
#' @param adjust Whether to adjust effort to be species-specific
#' @param buffer A sf polygon to be used as a buffer around locations to prevent extrapolation outside of the species range. Dual mesh cells without any effort outside of this buffer will be assigned an effort value to force model predictions toward 0.
#' @param nsimeff Effort value to assign to cells outside of the buffer (\code{integer} representing an umber of background observations).
#' @param \dots Arguments passed to \code{inla}
#'
#'
#' @details
#' Either points or a raster can be used. If a raster, empty cells will be filled with 0.
#'
#' @returns
#'
#' @references
#' Simpson, D. Illian, J. B., Lindgren, F. Sørbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 \url{https://doi.org/10.1093/biomet/asv064}
#'
#' @keywords
#'
#' @importFrom sf st_intersects st_transform st_crs
#' @importFrom exactextractr exact_extract
#'
#' @export
#'
#'
dmesh_effort<-function(dmesh,obs,background,adjust=FALSE,buffer=NULL, nsimeff=20){

  dm<-dmesh$dmesh
  nobs<-lengths(st_intersects(dm,obs))
  if(inherits(background,"SpatRaster")){
    nbackground<-exact_extract(
      background,st_transform(dm,st_crs(background)),
      fun=function(values,coverage_fractions){
        sum(values*coverage_fractions,na.rm=TRUE)/sum(coverage_fractions,na.rm=TRUE)
      }
    )
  }else{
    nbackground<-lengths(st_intersects(dm,background))
  }
  miss<-is.na(nbackground) | is.nan(nbackground)
  if(any(miss)){
    nbackground[miss]<-0
  }

  if(adjust){ # to complete

    o<-st_intersects(background,dm)
    splist<-data.frame(species=background$species,id=dm$id[unlist(o)])
    res<-aggregate(species~id,data=splist,fun=function(i){length(unique(i))})

    pres<-as.integer(nobs>0)
    vals<-nobs*scales::rescale(pres/nbsp,to=c(1,max(nbsp)^1))
  }



  if(!is.null(buffer)){
    o<-!as.logical(lengths(st_intersects(dmesh$dmesh,buff)))
    nbackground<-ifelse(o & nbackground==0L,nsimeff,nbackground)
  }

  dmesh[["effort"]]<-data.frame(nobs,nbackground)
  dmesh
  #nbobs<-obs[,.(nbobs=.N),by=dmesh]
  #nbsp<-obs[,.(nbsp=length(unique(species))),by=.(dmesh)]
  #dmesh$dmesh<-dmesh$id
  #dmesh<-merge(dmesh,nbobs,all.x=TRUE)
  #dmesh$nbobs[is.na(dmesh$nbobs)]<-0
  #dmesh<-merge(dmesh,nbsp,all.x=TRUE)
  #dmesh$nbsp[is.na(dmesh$nbsp)]<-0

  ### Species effort
  #temp<-spobs[,.(spobs=.N),by=.(dmesh)]
  #dmesh<-merge(dmesh,temp,all.x=TRUE)
  #dmesh$spobs[is.na(dmesh$spobs)]<-0
  #dmesh$pres<-as.integer(dmesh$spobs>0)
  #vals<-dmesh$nbobs*scales::rescale(dmesh$pres/dmesh$nbsp,to=c(1,max(dmesh$nbsp)^1))
  #vals<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
  #dmesh$effoccs<-vals

}
