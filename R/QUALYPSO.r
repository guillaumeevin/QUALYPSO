###===============================###===============================###
### Guillaume Evin
### 11/02/2019, Grenoble
###  IRSTEA
### guillaume.evin@irstea.fr
###
### These functions use data augmentation and Bayesian techniques for the assessment
### of single-member and incomplete ensembles of climate projections.
### It provides unbiased estimates of climate change responses of all
### simulation chains and of all uncertainty variables. It additionally propagates
### uncertainty due to missing information in the estimates.
###
### REFERENCES
###
### Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
### Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
### Journal of Climate. https://doi.org/10.1175/JCLI-D-18-0606.1.
###
###
###===============================###===============================###


#==============================================================================
#' get.Qstar.mat
#'
#' Provide matrix containing Helmert contrasts (see Eq. A7 in Evin et al., 2019).
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie. (2019) <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @param p integer
#'
#' @return \item{matrix}{p x (p-1) matrix containing Helmert contrasts}
#'
#' @author Guillaume Evin
get.Qstar.mat = function(p){
  # initialize matrix
  M = matrix(NA,nrow=p,ncol=(p-1))

  for(i in 1:p){ # by row
    for(j in 1:(p-1)){ # by column
      if(i<(p-j)){
        M[i,j] = 0
      }else if((i+j)==p){
        M[i,j] = -j
      }else{
        M[i,j] = 1
      }
    }
  }

  # return matrix
  return(M)
}


#==============================================================================
#' get.Qmat
#'
#' Provide matrix Q derived from a matrix of Helmert contrasts: \deqn{ Q = Q^* (Q^{*T} Q^*)^{-1/2}}
#' See Eq. A6 in Evin et al., 2019.
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie. (2019) <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @param p integer
#'
#' @return \item{matrix}{p x p matrix}
#'
#' @author Guillaume Evin
get.Qmat = function(p){
  # get Qstar
  Qstar = get.Qstar.mat(p)

  # get Q: Eq. A6 in Evin et al. (2019)
  Q = Qstar %*% MASS::ginv(expm::sqrtm(t(Qstar) %*% Qstar))

  return(Q)
}


#==============================================================================
#' fit.climate.response
#'
#' fit trends for each simulation chain of an ensemble of \code{nS} projections. Each projection can be a time series
#' of \code{nY} time steps (e.g. number of years).
#'
#' @param Y matrix ofsimulation chains: \code{nS} x \code{nY}
#' @param parSmooth smoothing parameter \code{spar} in \link[stats]{smooth.spline}: varies in [0,1]
#' @param indexReferenceYear index of the control year (e.g. 1990)
#' @param typeChangeVariable type of change variable: "abs" (absolute) or "rel" (relative)
#'
#' @return list with the following fields:
#' \itemize{
#'   \item \strong{phiStar}: climate change response
#'   \item \strong{etaStar}: internal variability
#'   \item \strong{phi}: raw trend obtained using \link[stats]{smooth.spline}
#'   \item \strong{climateResponse}: output from \link[stats]{smooth.spline}
#'   \item \strong{varInterVariability}: scalar, internal variability component of the MME
#' }
#'
#' @author Guillaume Evin
fit.climate.response = function(Y, parSmooth=1, indexReferenceYear=1, typeChangeVariable="abs"){

  # number of simulation chains
  nS = nrow(Y)

  # length of the simulation chains
  nY = ncol(Y)

  # prepare outputs
  phiStar = etaStar = phi = matrix(nrow=nS,ncol=nY)
  climateResponse = list()

  for(iS in 1:nS){
    # projection for this simulation chain
    Ys = Y[iS,]

    # fit a smooth signal (smooth cubic splines)
    zz = !is.na(Ys)
    phiS = rep(NA,nY)
    smooth.spline.out<-stats::smooth.spline((1:nY)[zz],Ys[zz],spar=parSmooth)
    climateResponse[[iS]] = smooth.spline.out
    phiS[zz] = smooth.spline.out$y

    # store climate response for this simulation chain
    phi[iS,] = phiS

    # climate response of the control year
    phiC = phiS[indexReferenceYear]

    # Climate change response: phiStar, and internal variability expressed as a change: etaStar
    if(typeChangeVariable=='abs'){
      # Eq. 5
      phiStar[iS,] = phiS-phiC
      etaStar[iS,] = Ys-phiS
    }else if(typeChangeVariable=='rel'){
      # Eq. 6
      phiStar[iS,] = phiS/phiC-1
      etaStar[iS,] = (Ys-phiS)/phiC
    }else{
      stop("fit.climate.response: argument type.change.var must be equal to 'abs' (absolute changes) or 'rel' (relative changes)")
    }
  }

  # Variance related to the internal variability: considered constant over the time period
  # (see Eq. 22 and 23 in Hingray and Said, 2014). We use a direct empirical estimate
  # of the variance of eta for each simulation chain and take the mean, see Eq. 19
  varInterVariability = mean(apply(etaStar,2,var),na.rm=T)

  # return objects
  return(list(phiStar=phiStar,etaStar=etaStar,phi=phi,climateResponse=climateResponse,varInterVariability=varInterVariability))
}




#==============================================================================
#' QUALYPSO.ANOVA.i
#'
#' Partition sources of uncertainty in climate change responses for one lead time 't' or grid point 'i'.
#'
#' @param phiStar.i vector of climate change response at time 't' or for grid point 'i': \code{nS x 1}
#' @param nMCMC number of MCMC simulation
#' @param listScenarioInput list containing specifications, provided by \code{QUALYPSO.process.scenario}
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{mu}: vector of length \code{nMCMC}, mean climate change response
#'   \item \strong{sigma2}: vector of length \code{nMCMC}, variance of the residual terms
#'   \item \strong{effect}: list with \code{nTypeEff} elements, where each element corresponds to a different type of effect (e.g. alpha, beta, gamma in Eq. 7)
#'   \item \strong{empEff}: list with \code{nTypeEmpEff} elements, where each element corresponds to an empirical effect
#'   Each element is a matrix \code{nMCMC} x \code{nMaineff}, and \code{nMaineff} is the number of main effects (e.g. number of GCMs, RCMs, etc.)
#' }
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. https://doi.org/10.1175/JCLI-D-18-0606.1.
#'
#' @author Guillaume Evin
QUALYPSO.ANOVA.i = function(phiStar.i, nMCMC, listScenarioInput){
  #============= retrieve objects related to the scenarios =====
  listEff=listScenarioInput$listEff
  scenAvail=listScenarioInput$scenAvail
  scenComp=listScenarioInput$scenComp
  nEff=listScenarioInput$nEff
  nTypeEff=listScenarioInput$nTypeEff
  nComp=listScenarioInput$nComp
  isMissing=listScenarioInput$isMissing
  nMissing=listScenarioInput$nMissing
  iMatchScen=listScenarioInput$iMatchScen
  indexEffInCompScen=listScenarioInput$indexEffInCompScen
  Qmat=listScenarioInput$Qmat
  objEmpEff=listScenarioInput$objEmpEf


  #=============  initialize effects =================
  effect.POST = list()
  if(all(phiStar.i==0)){
    # when i corresponds to the control period, phiStar.i = 0 and no ANOVA can be performed
    # return 0 instead of NA for representations
    mu.POST = sigma2.POST = rep(0,nMCMC)
    for(i.eff in 1:nEff) effect.POST[[i.eff]] = matrix(0,nrow=nMCMC,ncol=nTypeEff[i.eff])
    return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST))
  }else{
    # phi to complete
    phi2comp = phiStar.i[iMatchScen]
    # initialise matrix and arrays
    mu.POST = sigma2.POST = vector(length=nMCMC)
    for(i.eff in 1:nEff) effect.POST[[i.eff]] = matrix(nrow=nMCMC,ncol=nTypeEff[i.eff])
  }


  #=============  treat effects for which only empirical effects are computed if needed =================
  if(!is.null(listScenarioInput$objEmpEff)){
    nEmpEff = objEmpEff$nEmpEff
    nTypeEmpEff= objEmpEff$nTypeEmpEff
    indexEmpEffInCompScen = objEmpEff$indexEmpEffInCompScen
    empEff.POST = list()
    for(i.empEff in 1:nEmpEff)  empEff.POST[[i.empEff]] = matrix(data = 0, nrow=nMCMC,ncol=nTypeEmpEff[i.empEff])
  }else{
    empEff.POST = NULL
  }


  ###################################################################
  # Bayesian sampling is done with Gibbs algorithm. Details about the choice of the prior distributions,
  # the hyperparameters, and full conditional posterior distributions are given in Appendix in Evin et al. (2019)

  #============= hyper-parameters: see Appendix A, section g. =============
  mu.mu = mean(phiStar.i,na.rm=T)
  largevar = 16*var(as.vector(phiStar.i),na.rm=T)
  sig2.mu = largevar
  lam.sig = 0.5
  X.diff = phiStar.i-mean(phiStar.i,na.rm=T)
  for(i.eff in 1:nEff){
    eff.hat = aggregate(x = X.diff, by = list(scenAvail[,i.eff]), FUN = "mean")
    X.diff = X.diff - eff.hat$x[match(scenAvail[,i.eff],eff.hat$Group.1)]
  }
  nu.sig = 0.5*var(as.vector(X.diff),na.rm=T)
  s2eff = largevar


  #============ first iteration ===============

  # we sample from the prior distributions

  # error variance
  s2 = 1/rgamma(n=1, shape=lam.sig, rate=nu.sig) # Eq. A3
  sigma2.POST[1] = s2

  # grand mean
  mu = rnorm(n=1, mean=mu.mu, sd=sqrt(sig2.mu)) # Eq. A1
  mu.POST[1] = mu

  # main effects
  mat.eff = matrix(nrow=nComp,ncol=nEff)
  for(i.eff in 1:nEff){
    n = nTypeEff[i.eff]
    eff.star = rnorm(n=(n-1),sd=sqrt(s2eff)) # Eq A8
    eff = Qmat[[i.eff]] %*% eff.star # linear transform, see par. in Appedix A, section c. for the parameter beta
    mat.eff[,i.eff] = eff[indexEffInCompScen[,i.eff]]
    effect.POST[[i.eff]][1,] = eff
  }

  # missing values: sample from posterior
  mean.X = mu + Rfast::rowsums(mat.eff)
  sd.X = sqrt(s2)
  phi2comp[isMissing] = rnorm(n=nMissing, mean=mean.X[isMissing], sd=rep(sd.X,nMissing)) # Eq A14

  #============ iteration 2,... ===============
  for(i.MCMC in 2:nMCMC){

    #_________________ error variance ____________________
    sum.diff2 = sum((phi2comp - mu - Rfast::rowsums(mat.eff))^2)
    s2 = 1/rgamma(n=1, shape=nComp/2+lam.sig, rate=sum.diff2/2+nu.sig) # Eq. A5
    sigma2.POST[i.MCMC] = s2


    #__________________ grand mean _______________________
    mu.V = (1/s2)*sum(phi2comp) + mu.mu/sig2.mu
    mu.PSi = 1/(nComp/s2 + 1/sig2.mu)
    mu = rnorm(n=1, mean=mu.V*mu.PSi, sd=sqrt(mu.PSi)) # Eq. A2
    mu.POST[i.MCMC] = mu


    #________________ main effects _______________________
    # remove global mean: main effects + residual term
    X.shift = phi2comp-mu

    for(i.eff in 1:nEff){
      # number of factors for each type of effect (i.e. number of GCMs, RCMs)
      n = nTypeEff[i.eff]
      # marginal sums
      mar.diff = Rfast::group.sum(X.shift,indexEffInCompScen[,i.eff])
      # arguments of full conditional posterior distribution
      V = as.numeric((1/s2)*(t(Qmat[[i.eff]]) %*% mar.diff))
      Psi = 1/(nComp/(s2*n)+1/s2eff)
      eff.star = Psi*V+rnorm(n=(n-1))*sqrt(Psi) # Eq A9
      eff = Qmat[[i.eff]] %*% eff.star
      mat.eff[,i.eff] = eff[indexEffInCompScen[,i.eff]]
      effect.POST[[i.eff]][i.MCMC,] = eff
    }


    #_____ missing values: sample from posterior ________
    mean.X = mu + Rfast::rowsums(mat.eff)
    sd.X = sqrt(s2)
    phi2comp[isMissing] = rnorm(n=nMissing, mean=mean.X[isMissing], sd=rep(sd.X,nMissing)) # Eq A14


    #_____ empirical effects ________
    if(!is.null(objEmpEff)){
      # xi
      xi = phi2comp - mu - Rfast::rowsums(mat.eff)
      for(i.empEff in 1:nEmpEff){
        for(iType in 1:nTypeEmpEff[i.empEff]){
          selEmpEff = which(objEmpEff$indexEmpEffInAvailScen[,i.empEff] == iType)
          empEff.POST[[i.empEff]][i.MCMC,iType] = mean(xi[match(selEmpEff,iMatchScen)])
        }
      }

    }
  }


  # return MCM samples from the posterior
  return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST,empEff=empEff.POST))
}

#==============================================================================
#' QUALYPSO.process.scenario
#'
#' Process input scenarios.
#'
#' @param scenAvail matrix of available combinations
#' @param computeEmpEff vector of column indices in scenAvail corresponding to effects which are estimated empirically
#'
#' @return list of objects
#'
#' @author Guillaume Evin
QUALYPSO.process.scenario = function(scenAvail,computeEmpEff){
  # number of scenarios
  nS = nrow(scenAvail)

  # empirical interactions
  if(!is.null(computeEmpEff)){
    # check input
    if(!computeEmpEff%in%1:nS) stop('computeEmpEff must be a vector of column indices in scenAvail')

    # number of effects for which we compute empirical interactions
    nEmpEff = length(computeEmpEff)

    # prepare length of types of effect
    nTypeEmpEff = vector(length=nEmpEff)

    # list of unique effects
    listEmpEff = list()

    # find indices of the effects in available scenarios
    indexEmpEffInAvailScen = matrix(nrow=nS,ncol=nEmpEff)


    for(i.eff in 1:nEmpEff){
      indexEmpEff = computeEmpEff[i.eff]
      vecEmpEff = scenAvail[,indexEmpEff]
      # unique effects
      UniqueEmpEff  = unique(vecEmpEff)
      listEmpEff[[i.eff]] = UniqueEmpEff
      nTypeEmpEff[i.eff] = length(UniqueEmpEff)
      # find indices of the effects in available scenarios
      indexEmpEffInAvailScen[,i.eff] = match(vecEmpEff,UniqueEmpEff)
    }

    # modify scenAvail
    scenAvail = scenAvail[,-computeEmpEff]

    # return object
    objEmpEff = list(computeEmpEff=computeEmpEff,nEmpEff=nEmpEff,nTypeEmpEff=nTypeEmpEff,listEmpEff=listEmpEff,
                     indexEmpEffInAvailScen=indexEmpEffInAvailScen)
  }else{
    objEmpEff = NULL
  }

  # list of effects
  nEff = ncol(scenAvail)
  listEff = list()
  for(i in 1:nEff) listEff[[i]] = unique(scenAvail[,i])
  nTypeEff = unlist(lapply(listEff,length))

  # possible combinations of main effects
  scenComp = expand.grid(listEff)
  nComp = nrow(scenComp)


  #########  Missing scenarios   #########
  vScenComp <- apply(scenComp, 1, paste, collapse='.')
  vScenAvail <- apply(scenAvail, 1, paste, collapse='.')
  isMissing = !vScenComp%in%vScenAvail
  nMissing = sum(isMissing)


  #########   vector of projections to complete with data augmentation   #########
  iMatchScen = match(vScenComp,vScenAvail)


  #########   matrix of effect index: for each main effect, index of 'scenComp' (combinations of scenarios) related to listEff   #########
  indexEffInCompScen = matrix(nrow=nComp,ncol=nEff)
  for(i.eff in 1:nEff){
    indexEffInCompScen[,i.eff] = match(scenComp[,i.eff],listEff[[i.eff]])
  }


  #########   transformation matrices Q  #########
  Qmat = list()
  for(i.eff in 1:nEff){
    Qmat[[i.eff]] = get.Qmat(nTypeEff[i.eff])
  }


  return(list(listEff=listEff,
              scenAvail=scenAvail, scenComp=scenComp,
              nEff=nEff, nTypeEff=nTypeEff, nComp=nComp,
              isMissing=isMissing, nMissing=nMissing,
              iMatchScen=iMatchScen,
              indexEffInCompScen=indexEffInCompScen,
              Qmat=Qmat, objEmpEff=objEmpEff))
}


#==============================================================================
#' QUALYPSO.check.option
#'
#' Check possible options.
#'
#' @param listOption list of options
#'
#' @return listOption
#'
#' @author Guillaume Evin
QUALYPSO.check.option = function(listOption){
  if(is.null(listOption)){
    listOption = list()
  }

  # parSmooth
  if('parSmooth' %in% names(listOption)){
    parSmooth = listOption[['parSmooth']]
    if(!(is.numeric(parSmooth)&(parSmooth>0&parSmooth<=1))) stop('parSmooth must be in ]0,1]')
  }else{
    listOption[['parSmooth']] = 1
  }

  # typeChangeVariable
  if('typeChangeVariable' %in% names(listOption)){
    typeChangeVariable = listOption[['typeChangeVariable']]
    if(!(typeChangeVariable%in%c('abs','rel'))) stop("typeChangeVariable must be equal to 'abs' or 'rel'")
  }else{
    listOption[['typeChangeVariable']] = 'abs'
  }

  # nBurn
  if('nBurn' %in% names(listOption)){
    nBurn = listOption[['nBurn']]
    if(!(is.numeric(nBurn)&(nBurn>=0))) stop('wrong value for nBurn')
  }else{
    listOption[['nBurn']] = 1000
  }

  # nKeep
  if('nKeep' %in% names(listOption)){
    nKeep = listOption[['nKeep']]
    if(!(is.numeric(nKeep)&(nKeep>=0))) stop('wrong value for nKeep')
  }else{
    listOption[['nKeep']] = 2000
  }

  # nMCMC
  listOption$nMCMC = listOption$nKeep+listOption$nBurn


  # nCluster
  if('nCluster' %in% names(listOption)){
    nCluster = listOption[['nCluster']]
    if(!(is.numeric(nCluster)&(nCluster>=0))) stop('wrong value for nCluster')
  }else{
    listOption[['nCluster']] = 1
  }

  # doCompress
  if('doCompress' %in% names(listOption)){
    doCompress = listOption[['doCompress']]
    if(!(is.logical(doCompress))) stop('wrong value for doCompress')
  }else{
    listOption[['doCompress']] = TRUE
  }

  # computeEmpEff
  if('computeEmpEff' %in% names(listOption)){
    computeEmpEff = listOption[['computeEmpEff']]
    if(!(is.vector(computeEmpEff)&is.numeric(computeEmpEff))){
      stop('computeEmpEff must be a vector of column indices in scenAvail')
    }
  }else{
    listOption[['computeEmpEff']] = NULL
  }

  # doCompress
  if('quantileCompress' %in% names(listOption)){
    quantileCompress = listOption[['quantileCompress']]
    if(!(is.numeric(quantileCompress))) stop('wrong value for quantileCompress')
  }else{
    listOption[['quantileCompress']] = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)
  }

  # Version
  listOption[['version']] = 'v1.0.0'

  return(listOption)
}


#==============================================================================
#' QUALYPSO.ANOVA
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach
#'
#' @param phiStar array of climate response: \code{nS} x \code{nY}
#' @param scenAvail matrix of available combinations
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000)
#'   \item \strong{nKeep}: number of kept samples (default: 2000)
#'   \item \strong{nCluster}: number of clusters used for the parellelization
#'   \item \strong{doCompress}: logical, indicates if all the samples from the posterior distributions are stored (=FALSE)
#'   or if only quantiles are retrieved (=TRUE)
#'   \item \strong{computeEmpEff}: vector of column indices in scenAvail corresponding to effects which are estimated empirically
#' }
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{mu}: mean climate change response
#'   \item \strong{sigma2}: variance of the residual terms
#'   \item \strong{effect}: list with \code{nTypeEff} elements, where each element corresponds to a different type of effect (e.g. alpha, beta, gamma in Eq. 7)
#'   Each element contains the main effects (e.g. number of GCMs, RCMs, etc.)
#'   \item \strong{varEffect}: matrix \code{nTypeEff} x \code{nY} of variances related to the main effects
#'   \item \strong{varResidualEffect}: vector of length \code{nY} of variances of residual effects
#'   \item \strong{years}: number of years (nY)
#'   \item \strong{scenAvail}: matrix of available combinations given as inputs
#'   \item \strong{listEff}: list of effects contained in scenAvail
#'   \item \strong{listOption}: list of options used to obtained these results
#' }
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. https://doi.org/10.1175/JCLI-D-18-0606.1.
#'
#' @author Guillaume Evin
QUALYPSO.ANOVA = function(phiStar,scenAvail,listOption=NULL){
  #########  process input #########
  if(any(is.na(phiStar))){
    stop('na values are not allowed in phiStar')
  }

  # number of grid points / years
  n = dim(phiStar)[2]

  # Check list of options
  listOption = QUALYPSO.check.option(listOption)

  # number of MCMC samples
  nKeep = listOption$nKeep
  vec.keep = (listOption$nBurn+1):listOption$nMCMC

  # Process scenarios data.frame to get different objects
  listScenarioInput = QUALYPSO.process.scenario(scenAvail = scenAvail,computeEmpEff = listOption$computeEmpEff)
  nEff = listScenarioInput$nEff
  nTypeEff = listScenarioInput$nTypeEff

  #########  Apply ANOVA for each time step: parallelisation over time steps #########
  cl <- parallel::makeCluster(listOption$nCluster)
  doParallel::registerDoParallel(cl)
  i = NULL # avoid warning during check
  Anova.POST <- foreach(i=1:n, .multicombine=TRUE,
                        .export=c("QUALYPSO.ANOVA.i")) %dopar% {
                          POST.out = QUALYPSO.ANOVA.i(phiStar.i=phiStar[,i], nMCMC=listOption$nMCMC,listScenarioInput = listScenarioInput)
                          return(POST.out)
                        }
  parallel::stopCluster(cl)

  # number of years for which QE-ANOVA is applied
  y.POST = length(Anova.POST)

  #============================
  # POSTERIOR
  #============================

  # grand mean
  mu.POST = matrix(nrow=y.POST,ncol=nKeep)
  for(y in 1:y.POST) mu.POST[y,]=Anova.POST[[y]]$mu[vec.keep]

  # uncertainty climate response / interactions
  sigma2.POST = matrix(nrow=y.POST,ncol=nKeep)
  for(y in 1:y.POST) sigma2.POST[y,]=Anova.POST[[y]]$sigma2[vec.keep]

  # main effects
  eff.POST = list()
  for(i.eff in 1:nEff){
    # trim posterior distributions of the main effects
    eff.POST.i = array(dim=c(y.POST,nKeep,nTypeEff[i.eff]))
    for(y in 1:y.POST) eff.POST.i[y,,]=Anova.POST[[y]]$effect[[i.eff]][vec.keep,]
    eff.POST[[i.eff]] = eff.POST.i
  }

  POST = list(mu=mu.POST,sigma2=sigma2.POST,effect=eff.POST)


  #============================
  # VARIANCES
  #============================

  # variance of residual effects, With the Bayesian estimation, this variance
  # is part of the inference (sigma^2) and we take the mean of the posterior as the point estimate
  varResidualEffect = apply(sigma2.POST,1,mean)

  # variance of the main effects
  varEffect = matrix(nrow=nEff,ncol=y.POST)

  for(i.eff in 1:nEff){
    # predictive variance: mean of variances, which correspond to the mean of squares since the mean is 0 by constraint (see Eq 16, 17 and 18)
    varEffect[i.eff,] = apply(eff.POST[[i.eff]]^2,1,mean)
  }


  #============================
  # QUANTILES
  #============================

  # mean and s2
  mu.QUANT = t(apply(mu.POST,1,quantile,probs=listOption$quantileCompress))
  sigma2.QUANT = t(apply(sigma2.POST,1,quantile,probs=listOption$quantileCompress))

  # effects
  eff.QUANT = list()
  for(i.eff in 1:nEff){
    eff.quant.i.eff = apply(eff.POST[[i.eff]],c(1,3),quantile,probs=listOption$quantileCompress)
    eff.QUANT[[i.eff]] = aperm(eff.quant.i.eff, c(2,1,3))
  }

  QUANT = list(mu=mu.QUANT,sigma2=sigma2.QUANT,effect=eff.QUANT)


  #============================
  # MEAN
  #============================

  # mean main effects
  muMean = apply(mu.POST,1,mean)
  sigmaMean = apply(sigma2.POST,1,mean)
  effMean = list()
  for(i.eff in 1:nEff){
    effMean[[i.eff]] = apply(eff.POST[[i.eff]],c(1,3),mean)
  }

  MEAN = list(mu=muMean,sigma2=sigmaMean,effect=effMean)


  #============================
  # MAIN CHANGE BY EFFECT (MU + EFFECT)
  #============================

  meanChangeByEffect = quantChangebyEffect = list()
  eff.QUANT = list()
  for(i.eff in 1:nEff){
    eff.POST.i = eff.POST[[i.eff]]
    for(j in 1:dim(eff.POST.i)[3]){
      eff.POST.i[,,j] = eff.POST.i[,,j] + mu.POST
    }

    # mean
    meanChangeByEffect[[i.eff]] = apply(eff.POST.i,c(1,3),mean)

    # quantile
    eff.quant.i.eff = apply(eff.POST.i,c(1,3),quantile,probs=listOption$quantileCompress)
    quantChangebyEffect[[i.eff]] = aperm(eff.quant.i.eff, c(2,1,3))
  }

  MEAN$ChangeByEffect = meanChangeByEffect
  QUANT$ChangeByEffect = quantChangebyEffect


  #============================
  # EMPIRICAL EFFECTS
  #============================

  if(!is.null(listOption$computeEmpEff)){
    objEmpEff = listScenarioInput$objEmpEff
    nEmpEff = objEmpEff$nEmpEff
    nTypeEmpEff = objEmpEff$nTypeEmpEff

    #========= mean interactions =========
    empEff.POST = list()
    for(i.empEff in 1:nEmpEff){
      # trim posterior distributions of the main effects
      empEff.POST.i = array(dim=c(y.POST,nKeep,nTypeEmpEff[i.empEff]))
      for(y in 1:y.POST) empEff.POST.i[y,,]=Anova.POST[[y]]$empEff[[i.empEff]][vec.keep,]
      empEff.POST[[i.empEff]] = empEff.POST.i
    }
    POST$empEff = empEff.POST

    #========= quantiles =========
    empEff.QUANT = list()
    for(i.empEff in 1:nEmpEff){
      empEff.quant.i = apply(empEff.POST[[i.empEff]],c(1,3),quantile,probs=listOption$quantileCompress)
      empEff.QUANT[[i.empEff]] = aperm(empEff.quant.i, c(2,1,3))
    }
    QUANT$empEff = empEff.QUANT

    #========= mean =========
    empEff.Mean = list()
    for(i.empEff in 1:nEmpEff){
      empEff.Mean[[i.empEff]] = apply(empEff.POST[[i.empEff]],c(1,3),mean)
    }
    MEAN$empEff = empEff.Mean
  }

  # if doCompress is true, we do not return all the samples from the posterior distributions
  # which saves memory
  if(listOption$doCompress) POST=NULL

  # return results
  return(list(POSTERIOR=POST,QUANT=QUANT,MEAN=MEAN,varEffect=varEffect,varResidualEffect=varResidualEffect,
              years=y.POST,scenAvail=scenAvail,listOption=listOption,listEff=listScenarioInput$listEff))

}



#==============================================================================
#' QUALYPSO
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach
#'
#' @param Y array of climate response: \code{nS} x \code{nY}
#' @param scenAvail matrix of available combinations
#' @param vecYears vector of years corresponding to the projections
#' @param indexReferenceYear index corresponding to the control year (e.g. 1990)
#' @param futureYear index corresponding to a future year (e.g. 2090)
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{parSmooth}: smoothing parameter \code{spar} in \link[stats]{smooth.spline}: varies in [0,1]
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute) or "rel" (relative)
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000)
#'   \item \strong{nKeep}: number of kept samples (default: 2000)
#'   \item \strong{nCluster}: number of clusters used for the parellelization
#'   \item \strong{doCompress}: logical, indicates if all the samples from the posterior distributions are stored (=FALSE)
#'   or if only quantiles are retrieved (=TRUE)
#'   \item \strong{computeEmpEff}: vector of column indices in scenAvail corresponding to effects which are estimated empirically
#'
#' }
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{mu}: mean climate change response
#'   \item \strong{sigma2}: variance of the residual terms
#'   \item \strong{effect}: list with \code{nTypeEff} elements, where each element corresponds to a different type of effect (e.g. alpha, beta, gamma in Eq. 7)
#'   Each element contains the main effects (e.g. number of GCMs, RCMs, etc.)
#'   \item \strong{varEffect}: matrix \code{nTypeEff} x \code{nY} of variances related to the main effects
#'   \item \strong{varResidualEffect}: vector of length \code{nY} of variances of residual effects
#'   \item \strong{years}: number of years (nY)
#'   \item \strong{scenAvail}: matrix of available combinations given as inputs
#'   \item \strong{listEff}: list of effects contained in scenAvail
#'   \item \strong{listOption}: list of options used to obtained these results
#' }
#' @examples
#' scenAvail = cbind(c('GCM1','GCM1','GCM2'),c('RCM1','RCM2','RCM1'))
#' listOption = list(nBurn=100,nKeep=100)
#' QUALYPSOOUT = QUALYPSO(Y=dataQUALYPSO, scenAvail=scenAvail, listOption=listOption)
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. https://doi.org/10.1175/JCLI-D-18-0606.1.
#'
#' @export
#' @author Guillaume Evin
QUALYPSO = function(Y,scenAvail,vecYears=NULL,indexReferenceYear=NULL,futureYear=NULL,listOption=NULL){
  ######### Check inputs and assign default values ##########

  ######## Check list of options ########
  listOption = QUALYPSO.check.option(listOption)

  # dimensions
  d = dim(Y)

  if(length(d)==3){
    # Y is an array: GridPoints x Scenario x Time
    nG = d[1]
    nS = d[2]
    nY = d[3]

    # parallelization over space (e.g. on a grid)
    paralType = 'Grid'
    if(is.null(futureYear)) stop("Since Y is an array of dimension 3: 'futureYear' must be provided'")
  }else{
    # Y is a matrix: Scenario x Time
    nS = nrow(Y)
    nY = ncol(Y)

    # parallelization over time
    paralType = 'Time'
  }

  # vecYears
  if(!is.null(vecYears)){
    if(nY!=length(vecYears)) stop('length of vecYears different from the number of columns of Y')
    if(!all(is.numeric(vecYears))) stop('vecYears must be a vector of years')
  }else{
    vecYears = 1:nY
  }

  # indexReferenceYear
  if(!is.null(indexReferenceYear)){
    if(!any(indexReferenceYear==1:nY)) stop('wrong value for indexReferenceYear')
  }else{
    indexReferenceYear = 1
  }

  # futureYear
  if(paralType == 'Grid'){
    if(!any(futureYear==1:nY)) stop('wrong value for futureYear')
  }

  ##############################################
  # extract climate response
  if(paralType == 'Time'){
    climResponse = fit.climate.response(Y, parSmooth=listOption$parSmooth,indexReferenceYear=indexReferenceYear,
                                        typeChangeVariable=listOption$typeChangeVariable)

    # extract quantities from these fits
    phiStarAllTime = climResponse$phiStar
    etaStarAllTime = climResponse$etaStar
    phiAllTime = climResponse$phi

    # years when some simulation chains are missing (actually we could apply QUALYPSO to a different set of scenarios
    # each year but some configurations may fail if few scenarios are available, and it is difficult to identify in advance)
    hasNa = apply(phiStarAllTime,2,function(x) any(is.na(x)))
    # apply ANOVA after the reference year and when there is no na
    filterYearsAnova = 1:nY>indexReferenceYear&!hasNa
    if(!any(filterYearsAnova)) stop('Error in QUALYPSO: check indexReferenceYear and NAs in Y')

    # vector of years for ANOVA
    vecYearsANOVA = vecYears[filterYearsAnova]

    # phiStar for ANOVA
    phiStar = phiStarAllTime[,filterYearsAnova]

    # internal variability
    varInterVariability = rep(climResponse$varInterVariability,length(vecYearsANOVA))
  }else if(paralType == 'Grid'){
    climResponse = list()
    phiStarAllTime = etaStarAllTime = phiAllTime = array(dim=d)
    varInterVariability = vector(length=nG)
    for(g in 1:nG){
      climResponse[[g]] = fit.climate.response(Y[g,,], parSmooth=listOption$parSmooth,indexReferenceYear=indexReferenceYear,
                                               typeChangeVariable=listOption$typeChangeVariable)

      # extract quantities from these fits
      phiStarAllTime[g,,] = climResponse[[g]]$phiStar
      etaStarAllTime[g,,] = climResponse[[g]]$etaStar
      phiAllTime[g,,] = climResponse[[g]]$phi
      varInterVariability[g] = climResponse[[g]]$varInterVariability
    }

    # phiStar for ANOVA
    phiStar = t(phiStarAllTime[,,futureYear])

    # vector of years for ANOVA
    vecYearsANOVA = vecYears[futureYear]
  }


  #############################################
  # ANOVA on phiStar
  anova = QUALYPSO.ANOVA(phiStar = phiStar, scenAvail = scenAvail, listOption = listOption)


  #############################################
  # format variance
  ANOVAVARIANCE = list()
  ANOVAVARIANCE$eff = anova$varEffect
  ANOVAVARIANCE$ResidualEffect = anova$varResidualEffect
  ANOVAVARIANCE$InterVariability = varInterVariability
  ANOVAVARIANCE$TotalVar = Rfast::colsums(anova$varEffect) + anova$varResidualEffect + varInterVariability


  # return results
  return(list(CLIMATEESPONSE=list(phiStar=phiStarAllTime,etaStar=etaStarAllTime,phi=phiAllTime),
              ANOVAPOST=anova$POSTERIOR,ANOVAQUANT=anova$QUANT,ANOVAMEAN=anova$MEAN,
              ANOVAVARIANCE=ANOVAVARIANCE,vecYears=listOption$vecYears,vecYearsANOVA=vecYearsANOVA,
              Y=Y,scenAvail=scenAvail,listEff=anova$listEff,listOption=anova$listOption))
}
