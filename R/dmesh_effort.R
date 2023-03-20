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

  # https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
  # rescale values to the a-b range
  rescale_ab<-function(x,a,b){
    ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
  }

  dm<-dmesh$dmesh
  nobs<-lengths(st_intersects(dm,obs))
  if(inherits(background,"SpatRaster")){
    if(adjust){
      warning("Species-specific adjustments ignored when background is SpatRaster")
    }
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

  if(adjust && inherits(background,"sf")){
    o<-st_intersects(background,dm)
    splist<-data.frame(species=background$species,id=dm$id[unlist(o)])
    res<-aggregate(species~id,
                   data=splist,FUN=function(i){
                     length(unique(i))
                   }
    )
    nsp<-res$species[match(dm$id,res$id)]
    nsp[is.na(nsp)]<-0
    pres<-as.integer(nobs>0)
    vals<-nobs*rescale_ab(pres/nsp,a=1,b=max(nsp))
    nbackgroundadjusted<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
  }else{
    nbackgroundadjusted<-nbackground
  }

  if(!is.null(buffer)){
    o<-!as.logical(lengths(st_intersects(dmesh$dmesh,buffer)))
    nbackgroundadjusted<-ifelse(o & nbackgroundadjusted==0L,nsimeff,nbackgroundadjusted)
  }

  if(adjust && inherits(background,"sf")){
    eff<-data.frame(nobs,nbackground,pres,nsp,nbackgroundadjusted)
  }else{
    eff<-data.frame(nobs,nbackground,nbackgroundadjusted)
  }

  dmesh[["effort"]]<-eff
  dmesh

}
