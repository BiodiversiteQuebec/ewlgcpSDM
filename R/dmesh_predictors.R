#'
#' @title
#' Summarize predictors to dual mesh
#'
#' @description
#' This function summarizes predictors to dual mesh. Predictors are assumed to be numerical values.
#'
#' @param dmesh A dual mesh
#' @param predictors A SpatRaster with predictor values. If processing in parallel, the SpatRastre needs to be wrapped using \link{\code{terra::wrap}}.
#'
#'
#' @details
#' Parallelization can be enabled by choosing a plan from the future package and the mesh will be splitted in chunks detemined by the number of workers. Don't forget to reset the plan to sequential when done.
#'
#' When the predictors do not cover the extent of the dual mesh, empty cells will be filled with the value of the nearest neighbour cell.
#'
#' @returns
#' A list to which an element predictors is added with a \code{data.frame} of predictor values.
#'
#' @examples
#' none
#'
#' @importFrom exactextractr exact_extract
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_apply
#' @importFrom sf st_centroid st_coordinates st_transform st_crs
#' @importFrom terra unwrap
#' @importFrom data.table as.data.table
#' @importFrom FNN knnx.index
#'
#' @export
#'
#'
#'
dmesh_predictors<-function(dmesh,predictors){
  cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
  if(cores>1){
    if(!inherits(predictors,"PackedSpatRaster")){
      stop("The SpatRaster of predictors needs to be wrapped for parallel processing. See ?terra::wrap.")
    }
    predictors<-unwrap(predictors)
  }
  dm<-dmesh$dmesh
  if(!(st_crs(dm)==st_crs(predictors))){
    stop("Dual mesh and predictors have different crs")
  }
  chunks <- split(1:nrow(dm), rep(1:cores, each=ceiling(nrow(dm)/cores))[1:nrow(dm)])
  options(future.globals.maxSize = 1000 * 1024 ^ 2)
  res<-future_lapply(chunks,function(chunksi){
    res<-exact_extract(predictors,
                    dm[chunksi,],
                    fun = function(values, coverage_fraction){
                      vals<-as.matrix(values)
                      covs<-matrix(rep(coverage_fraction,ncol(vals)),ncol=ncol(vals))
                      covs[is.na(vals)]<-NA
                      colSums(vals*covs,na.rm=TRUE)/colSums(covs,na.rm=TRUE)
                    },
                    force_df = FALSE,
                    progress = TRUE)
    if(is.null(dim(res))){
      res<-matrix(res,nrow=1)
      dimnames(res)[[1]]<-names(predictors)
    }
    t(res)
  })
  #plan(sequential)
  res<-as.data.frame(do.call("rbind",res))

  #dm<-cbind(dm,res)
  #plot(dm["tmean"])
  #plot(predictors$tmean)


  ### Fill mesh where values are missing
  xy<-st_coordinates(st_centroid(dm))
  xy<-cbind(xy,res)
  xynotna<-xy |> as.data.table() |> na.omit() |> as.matrix()
  nn<-knnx.index(xynotna[,1:2],xy[,1:2],k=1)
  res<-as.data.frame(xynotna[nn,-(1:2),drop=FALSE])

  #dm$tmean<-res2[,"elevation"]
  #plot(dm["tmean"])
  
  ### do not scale dummy variables that don't vary
  w<-which(apply(res,2,sd)!=0)
  if(any(w)){
    res[w]<-lapply(res[w],scale)
  }
  
  trans<-as.matrix(res)^2
  colnames(trans)<-paste0(colnames(trans),"2")
  res<-cbind(res,as.data.frame(trans))
  res<-res[,order(colnames(res))]


  dmesh[["predictors"]]<-res
  dmesh
}
