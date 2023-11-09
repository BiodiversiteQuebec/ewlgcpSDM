#'
#' @title
#' Effort-Weighted Log Gaussian Cox Process SDM
#'
#' @description
#' Runs the effort-weighted Log Gaussian Cox Process species distribution model
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...}.
#' @param dmesh A list with elements produced by the different \code{dmesh_} functions.
#' @param effort Logical. Whether to adjust the model for effort. Default \code{TRUE}.
#' @param adjust Logical. Whether to adjust the effort for being species specific. Default \code{TRUE}.
#' @param buffer Logical. Whether to add the effort buffer to help reducing prediction outside of the species range. Default \code{TRUE}.
#' @param orthogonal Logical. Whether to make the spatial field orthogonal to the predictors? Default \code{TRUE}.
#' @param prior.beta Normal priors for the betas of the fixed effects coefficients as required by \code{INLA}. Default is \code{list(prec=list(default=1/(1)^2,Intercept=1/(20)^2),mean=list(default=0,Intercept=0))} which means a prior with \code{mean = 0} and \code{sd = 1} for all coefficients and a prior with \code{mean = 0} and \code{sd = 20} for the model intercept. See \code{\link{?control.fixed}}.
#' @param prior.range Penalized complexity prior for the range of the spatial field. A vector of length two giving the probability that the range is inferior to a given value. The default is \code{prior.range = c(50, 0.01)} which represents a 1% chance that the range is inferior to 50 (in the units of the crs used).
#' @param prior.sigma Penalized complexity prior for the standard deviation (sd) of the spatial field. A vector of length two giving the probability that the sd is superior to a given value. The default is \code{prior.sigma = c(1, 0.01)} which represents a 1% chance that the range is superior to 1.
#' @param smooth x
#' @param \dots Further arguments to pass to \code{inla}
#'
#' @details
#' none
#'
#' @returns
#' A model of class \code{INLA}.
#'
#' @encoding
#' UTF-8
#'
#' @references
#' Simpson, D. Illian, J. B., Lindgren, F. S\u00f8rbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 \url{https://doi.org/10.1093/biomet/asv064}
#'
#' Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. 2019 Constructing Priors that Penalize the Complexity of Gaussian Random Fields. Journal of the American Statistical Association, 114(525): 445-452 \url{https://doi.org/10.1080/01621459.2017.1415907}
#'
#' @examples
#' none
#'
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
ewlgcp <- function(
    formula,
    dmesh,
    effort = TRUE,
    adjust = TRUE,
    buffer = TRUE,
    orthogonal = TRUE,
    prior.beta = NULL,
    prior.range=c(50,0.1),
    prior.sigma=c(1,0.1),
    smooth = 3/2,
    ...
){

  ### Params
  #f<-formula(paste("y~",paste(vars,collapse="+")))
  #formula<-f
  #obs<-st_as_sf(spobs,coords=c("x","y"),crs=crsr)
  #explanaMesh<-explana


  #model.arguments<-as.list(environment())
  dots<-list(...)

  vars<-all.vars(formula[[3]])

  if(is.null(prior.beta)){
    prior.beta<-list(prec=list(default=1/(1)^2,Intercept=1/(20)^2),mean=list(default=0,Intercept=0))
  }

  inla.defaults<-list(
    num.threads=2:2,
    #blas.num.threads=2,
    control.inla=list(
      strategy="adaptive",
      int.strategy="eb",
      huge=TRUE,
      control.vb=list(
        enable=TRUE,
        verbose=TRUE
      )
    ),# adaptive, eb
    inla.mode="experimental",
    control.fixed=prior.beta,
    control.compute=list(config=TRUE,openmp.strategy="pardiso"),
    verbose=TRUE
  )

  inla.arguments<-c(dots,inla.defaults)
  inla.arguments<-inla.arguments[!duplicated(names(inla.arguments))]



  #miss<-which(!names(inla.arguments)%in%names(arguments))
  #if(any(miss)){
  #  arguments<-c(arguments,inla.arguments[w])
  #}



  #num.threads=2:2,
  #blas.num.threads=2,
  #control.inla=list(strategy="adaptive",int.strategy="eb",huge=TRUE,
  #                  control.vb=list(enable=TRUE, verbose=verbose)),
  #inla.mode="experimental",
  #control.fixed=prior.beta,
  #control.compute=list(config=TRUE,openmp.strategy="pardiso"),
  #verbose=verbose

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

  if((adjust || buffer) && !effort){
    warning("When effort = FALSE, adjust and buffer are ignored")
  }

  if(effort){
    if(adjust && buffer){
      if(is.null(dmesh$effort$nbackgroundspadjustedwithbuff)){
        stop("Missing effort info from output of dmesh_effort")
      }else{
        eff<-dmesh$effort$nbackgroundspadjustedwithbuff
      }
    }
    if(adjust && !buffer){
      if(is.null(dmesh$effort$nbackgroundspadjusted)){
        stop("Missing effort info from output of dmesh_effort")
      }else{
        eff<-dmesh$effort$nbackgroundspadjusted
      }
    }
    if(!adjust && !buffer){
      if(is.null(dmesh$effort$nbackground)){
        stop("Missing effort info from output of dmesh_effort")
      }else{
        eff<-dmesh$effort$nbackground
      }
    }
    if(!adjust && buffer){
      if(is.null(dmesh$effort$nbackgroundwithbuff)){
        stop("Missing effort info from output of dmesh_effort")
      }else{
        eff<-dmesh$effort$nbackgroundwithbuff
      }
    }
    k<-dmesh$weights > 0
    e<-eff[k]
    dmesh$weights[k]<-dmesh$weights[k]*((e/dmesh$weights[k])/max(e/dmesh$weights[k]))
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

  model.arguments<-list(
    formula=formule,
    family = "poisson",
    data = inla.stack.data(Stack),
    control.predictor = list(A = inla.stack.A(Stack),link = 1),
    E = inla.stack.data(Stack)$e
  )

  #verbose<-TRUE
  #model <- inla(formula = formule,
  #              family = "poisson",
  #              data = inla.stack.data(Stack),
  #              control.predictor = list(A = inla.stack.A(Stack),link = 1),
  #              E = inla.stack.data(Stack)$e,
  #              num.threads=2:2,
  #              blas.num.threads=2,
  #              control.inla=list(strategy="adaptive",int.strategy="eb",huge=TRUE,
  #                                control.vb=list(enable=TRUE, verbose=verbose)),
  #              inla.mode="experimental",
  #              control.fixed=prior.beta,
  #              control.compute=list(config=TRUE,openmp.strategy="pardiso"),
  #              verbose=verbose
  #)


  model<-do.call("inla",c(model.arguments,inla.arguments))

  nameRes <- names(model)

  #prog(message = sprintf("Model done %s", species[j]))

  #=============
  # Return model
  #=============

  attributes(model) <- list(formula = formula,
                            obs = dmesh$effort,
                            XEst = XEst,
                            XPred = XPred,
                            meshSpace = dmesh$mesh,
                            Stack = Stack)

  names(model) <- nameRes



  #class(model) <- "ppSpace"

  model


}

















