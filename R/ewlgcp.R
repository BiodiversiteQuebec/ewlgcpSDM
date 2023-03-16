#'
#' @title
#' Effort-Weighted Log Gaussian Cox Process SDM
#'
#' @description
#' Runs the effort-weighted Log Gaussian Cox Process species distribution model
#'
#' @param formula A formula of the form \code{y ~ x1 + x2}.
#' @param dmesh A dual mesh as a \code{MULTIPOLYGONS} \code{sf} object
#' @param effort Logical. Whether to adjust the model for effort. Default \code{TRUE}.
#' @param adjust Logical. Whether to adjust the effort for being species specific. Default \code{TRUE}.
#' @param buffer Logical. Whether to add a fictious effort to help reducing prediction outside of the species range. Default \code{TRUE}.
#' @param orthogonal Logical. Whether to make the spatial field orthogonal to the predictors? Default \code{TRUE}.
#' @param prior.beta x
#' @param prior.range x
#' @param prior.sigma x
#' @param smooth x
#' @param \dots x
#'
#' @details
#'
#'
#' @returns A model of class \code{INLA}.
#'
#' @references
#' Simpson, D. Illian, J. B., Lindgren, F. Sørbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 \url{https://doi.org/10.1093/biomet/asv064}
#'
#'
#' @examples
#'
#' add(10, 1)
#'
#'
#' @importFrom terra xyFromCell
#' @importFrom terra extract
#' @importFrom terra values
#' @importFrom sf st_coordinates
#' @importFrom INLA inla
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom INLA inla.spde.make.A
#' @importFrom INLA inla.spde.make.index
#' @importFrom INLA inla.stack
#' @importFrom INLA inla.stack.data
#' @importFrom INLA inla.stack.A
#' @importFrom Matrix Diagonal
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' @importFrom exactextractr exact_extract
#'
#'
#' @export
#'
#'
#'
ewlgcp <- function(formula, dmesh, effort = TRUE, adjust = FALSE, buffer = TRUE, orthogonal = TRUE, prior.beta = NULL, prior.range = NULL, prior.sigma = NULL, smooth = 2, ...) {

  ### Params
  #f<-formula(paste("y~",paste(vars,collapse="+")))
  #formula<-f
  #obs<-st_as_sf(spobs,coords=c("x","y"),crs=crsr)
  #explanaMesh<-explana
  bpriors<-list(prec=list(default=1/(0.5)^2,Intercept=1/(20)^2),mean=list(default=0,Intercept=0))


  vars<-all.vars(formula[[3]])

  prior.range<-c(50,0.1)
  prior.sigma<-c(1,0.1)
  prior.beta<-NULL
  smooth<-3/2
  num.threads<-1:1
  blas.num.threads<-1
  control.inla<-list(strategy="adaptive",int.strategy="eb",huge=TRUE) # adaptive, eb
  inla.mode<-"experimental"
  control.fixed<-NULL
  control.compute<-list(config=TRUE,openmp.strategy="pardiso.parallel")
  verbose<-TRUE


  #==============
  # Basic objects
  #==============
  #nobs <- length(obs)
  #nobs <- nrow(obs)

  #============
  # Define SPDE
  #============
  SPDE <- inla.spde2.pcmatern(mesh=dmesh$mesh,
                              alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.sigma,
                              constr=TRUE)

  #====================================================
  # Rescale weights if sampling bias offset is included
  #====================================================

  #dmesh$areas<-as.numeric(st_area(dmesh))
  #dmesh$weights<-dmesh$areas
  eff<-dmesh$effort$nbackground

  if(adjust){
    # adjust effort here for species specific
  }

  if(effort){
    k <- dmesh$weights > 0
    e <- eff[k]
    dmesh$weights[k]<-dmesh$weights[k] * ((e/dmesh$weights[k])/max(e/dmesh$weights[k]))
  }


  #========================
  # Define response objects
  #========================

  #spaceAgg <- mapSpecies:::aggData(obs, meshSpace = explanaMesh$meshSpace, meshDual = attributes(ppWeight)$dmesh)
  #spaceAgg <- data.frame(space = 1:nrow(dmesh), Freq = dmesh$spobs)
  # Pseudo-absences are the number of edges on the mesh
  # Occurences are the number of points
  yPP <- dmesh$effort$nobs
  # weight associated to pseudo-absences (w) and occurrences (0)
  ePP <- dmesh$weights#[spaceAgg$space]

  XEst<-dmesh$predictors[,vars,drop=FALSE]
  #XEst<-apply(XEst,2,function(i){scales::rescale(i,0:1)})
  XPred <- XEst


  #=====================
  # Fix given predictors
  #=====================
  if(FALSE){ # to fix predictors with arg fix
    m <- match(fix, colnames(XPred))
    v <- apply(XPred[, m, drop=FALSE], 2, max, na.rm = TRUE)
    XPred[, m] <- rep(v, each = nrow(XPred))
  }

  #================================================
  # Define projection matrix and build stack object
  #================================================
  # Projection matrix
  ProjInfer <- inla.spde.make.A(dmesh$mesh,loc = dmesh$mesh$loc)
  IDSpace <- inla.spde.make.index('i', dmesh$mesh$n)

  # Build stack objects
  StackEst <- inla.stack(data = list(y = yPP, e = ePP),
                         A = list(1, ProjInfer),
                         effects = list(c(list(Intercept = 1),
                                          asplit(XEst, 2)),
                                        IDSpace),
                         tag = "est")

  StackPred <- inla.stack(data = list(y = NA, e = NA),
                          A = list(1, ProjInfer),
                          effects = list(c(list(Intercept = 1),
                                           asplit(XPred, 2)),
                                         IDSpace),
                          tag = "pred")

  Stack <- inla.stack(StackEst, StackPred)


  X <- paste(colnames(XEst),collapse=" + ")
  fixed <- paste("y ~ 0 + Intercept +",X)
  formule <- formula(paste(fixed,"f(i, model=SPDE)", sep=" + "))

  if(TRUE){
    # build constraints
    XX = cbind(rep(1, nrow(XEst)), XEst)
    Q = qr.Q(qr(XX))
    AA <- as.matrix(t(Q)%*%ProjInfer)
    ee <- rep(0, ncol(XX))
    formule <- formula(paste(fixed,
                             "f(i, model=SPDE, extraconstr = list(A = AA, e = ee))",
                             sep=" + "))
  }




  verbose<-TRUE
  model <- inla(formule,
                family = "poisson",
                data = inla.stack.data(Stack),
                control.predictor = list(A = inla.stack.A(Stack),link = 1),
                E = inla.stack.data(Stack)$e,
                num.threads=2:2,
                blas.num.threads=2,
                control.inla=list(strategy="adaptive",int.strategy="eb",huge=TRUE,
                                  control.vb=list(enable=TRUE, verbose=verbose)),
                inla.mode="experimental",
                control.fixed=bpriors,
                control.compute=list(config=TRUE,openmp.strategy="pardiso"),
                verbose=verbose
  )


  nameRes <- names(model)

  #prog(message = sprintf("Model done %s", species[j]))

  #=============
  # Return model
  #=============

  attributes(model) <- list(formula = formula,
                            obs = obs,
                            XEst = XEst,
                            XPred = XPred,
                            meshSpace = dmesh$mesh,
                            Stack = Stack)

  names(model) <- nameRes



  #class(model) <- "ppSpace"

  model


}

















