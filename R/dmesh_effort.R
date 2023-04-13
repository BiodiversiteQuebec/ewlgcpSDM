#'
#' @title
#' Summarize effort to dual mesh
#'
#' @description
#' This function summarizes effort to dual mesh.
#'
#' @param dmesh A dual mesh sf object
#' @param obs A sf spatial data.frame with observations of the target species
#' @param background A sf data.frame with observations of the target group or a terra raster with the sum of observations of the target group for each pixel.
#' @param adjust Whether to adjust effort to be species-specific. Default to \code{FALSE}. If \code{adjust = TRUE}, a column named \"species\" with species name must be present in the \code{background data.frame}. Currently ignored for when a species column is not given in \code{background} or when it is a raster.
#' @param buffer A sf polygon to be used as a buffer around locations to prevent extrapolation outside of the species range. Dual mesh cells without any effort outside of this buffer will be assigned an effort value to force model predictions toward 0.
#' @param nsimeff Effort value to assign to cells outside of the buffer. An \code{integer} representing a number of background observations.
#' @param \dots Arguments passed to \code{inla}
#'
#'
#' @details
#' Either an \code{sf} spatial object with points or a raster can be used. If a raster, empty cells will be filled with 0.
#'
#' @returns
#' A list with element effort with a \code{data.frame} summarizing the number of observations and the various effort measures for each cell of the dual mesh.
#'
#' Depending on the options chosen, the \code{data.frame} will contain some or all of the following:
#'
#' \itemize{
#'   \item\code{nobs}{ : number of observations}
#'   \item\code{nbackground}{ : number of background observations from the target group}
#'   \item\code{npres}{ : whether the species is present in a cell or not (1 = present, 0 = absent)}
#'   \item\code{nsp}{ : number of species in a cell}
#'   \item\code{nbackgroundwithbuff}{ : nbackground to which fictious observations have been added using the effort buffer}
#'   \item\code{nbackgroundspadjusted}{ : nbackground with species-specific adjustment}
#'   \item\code{nbackgroundspadjustedwithbuff}{ : nbackground with species-specific adjustment and to which fictious observations have been added using the effort buffer}

#' }
#'
#' @references
#' Simpson, D. Illian, J. B., Lindgren, F. Sørbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 \url{https://doi.org/10.1093/biomet/asv064}
#'
#' @importFrom sf st_intersects st_transform st_crs
#' @importFrom exactextractr exact_extract
#'
#' @export
#'
#'
dmesh_effort<-function(
    dmesh,
    obs,
    background,
    adjust=FALSE,
    buffer=NULL,
    nsimeff=20
){

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

  leff<-list()
  leff[["nobs"]]<-nobs
  leff[["nbackground"]]<-nbackground

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
    vals<-nbackground*rescale_ab(pres/nsp,a=1,b=max(nsp))
    nbackgroundspadjusted<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
    leff[["nsp"]]<-nsp
    leff[["pres"]]<-pres
    leff[["nbackgroundspadjusted"]]<-nbackgroundspadjusted
  }

  if(!is.null(buffer)){
    o<-!as.logical(lengths(st_intersects(dmesh$dmesh,buffer)))
    if(any("nbackground"==names(leff))){
      leff[["nbackgroundwithbuff"]]<-ifelse(o & nbackground==0L,nsimeff,nbackground)
    }
    if(any("nbackgroundspadjusted"==names(leff))){
      leff[["nbackgroundspadjustedwithbuff"]]<-ifelse(o & nbackgroundspadjusted==0L,nsimeff,nbackgroundspadjusted)
    }
  }

  #if(adjust && inherits(background,"sf")){
  #  eff<-data.frame(nobs,nbackground,pres,nsp,nbackgroundspadjusted)
  #}else{
  #  if(!is.null(buffer)){
  #    eff<-data.frame(nobs,nbackground,nbackgroundspadjusted)
  #  }else{
  #    eff<-data.frame(nobs,nbackground)
  #  }
  #}

  #eff<-data.frame(leff)

  dmesh[["effort"]]<-data.frame(leff)
  dmesh

}
