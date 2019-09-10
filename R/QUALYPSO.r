###===============================###===============================###
### Guillaume Evin
### 02/05/2019, Grenoble
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
### references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
### Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
### Journal of Climate. J. Climate, 32, 2423–2440. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
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
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
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
#' Provide matrix Q derived from a matrix Q* of Helmert contrasts: \deqn{Q = Q^* (Q^{*T} Q^*)^{-1/2}}
#' See Eq. A6 in Evin et al., 2019.
#'
#' @param p integer
#'
#' @return \item{matrix}{p x p matrix}
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
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
#' Fit trends for each simulation chain of an ensemble of \code{nS} projections. Each simulation chain is a time series
#' of \code{nY} time steps (e.g. number of years).
#'
#' @param Y matrix of simulation chains: \code{nS} x \code{nY}
#' @param parSmooth smoothing parameter \code{spar} in \code{\link[stats]{smooth.spline}}: varies in [0,1]
#' @param indexReferenceYear index of the reference year
#' @param typeChangeVariable type of change variable: "abs" or "rel"
#'
#' @return list with the following fields for each simulation chain:
#' \itemize{
#'   \item \strong{phiStar}: climate change response
#'   \item \strong{etaStar}: internal variability
#'   \item \strong{phi}: raw trend obtained using \link[stats]{smooth.spline}
#'   \item \strong{climateResponse}: output from \link[stats]{smooth.spline}
#'   \item \strong{varInterVariability}: scalar, internal variability component of the MME
#' }
#'
#' @details
#' See \code{\link{QUALYPSO}} for further information on arguments \code{indexReferenceYear} and \code{typeChangeVariable}.
#'
#' @export
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
fit.climate.response = function(Y, parSmooth, indexReferenceYear, typeChangeVariable){

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

    # climate response of the reference year
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
#' Partition sources of uncertainty in climate change responses for one lead time or one grid point.
#'
#' @param phiStar.i vector of \code{nS} climate change response for one lead time or for one grid point: \code{nS x 1}
#' @param nMCMC number of MCMC simulation required
#' @param listScenarioInput list containing specifications, provided by \code{\link{QUALYPSO.process.scenario}}
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
#' Journal of Climate. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
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

  if(any(is.na(phiStar.i))){
    # when i corresponds to the reference period, phiStar.i = 0 and no ANOVA can be performed
    # return 0 instead of NA for representations
    mu.POST = sigma2.POST = rep(0,nMCMC)
    for(i.eff in 1:nEff) effect.POST[[i.eff]] = matrix(0,nrow=nMCMC,ncol=nTypeEff[i.eff])
    return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST))
  }

  if(all(phiStar.i==0)){
    # when i corresponds to the reference period, phiStar.i = 0 and no ANOVA can be performed
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
    eff.star = rnorm(n=(n-1),sd=sqrt(s2eff)) # Eq. A8
    eff = Qmat[[i.eff]] %*% eff.star # linear transform, see Appendix A for the parameter beta
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
#' @param scenAvail matrix of available combinations \code{nS} x \code{nEff}
#' @param computeEmpEff vector of column indices in \code{scenAvail} corresponding to effects which are estimated empirically
#'
#' @return list of preprocessed objects (\code{listEff, scenAvail, scenComp, nEff, nTypeEff, nComp, isMissing, nMissing, iMatchScen,
#' indexEffInCompScen, Qmat, EmpEff})
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
#' Check if input options provided in \code{\link{QUALYPSO}} are valid and assigned default values if missing.
#'
#' @param listOption list of options
#'
#' @return List containing the complete set of options.
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

  # quantileCompress
  if('quantileCompress' %in% names(listOption)){
    quantileCompress = listOption[['quantileCompress']]
    if(!(is.numeric(quantileCompress))) stop('wrong value for quantileCompress')
  }else{
    listOption[['quantileCompress']] = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)
  }

  # Version
  listOption[['version']] = 'v1.1'

  return(listOption)
}


#==============================================================================
#' QUALYPSO.ANOVA
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach.
#'
#' @param phiStar matrix of climate change responses (absolute or relative changes): \code{nS} x \code{n}. \code{n} can be the number of time steps or the number of grid points
#' @param scenAvail matrix of available combinations \code{nS} x \code{nEff}
#' @param listOption list of options (see \code{\link{QUALYPSO}})
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{POSTERIOR}: list of MCMC samples representing the posterior distributions of inferred quantities. \code{=NULL} if\code{listOption$doCompress=T}
#'   \item \strong{QUANT}: list of quantiles from the posterior distributions of inferred quantities
#'   \item \strong{MEAN}: list of mean of the posterior distributions of inferred quantities
#'   \item \strong{varEffect}: matrix \code{nEff} x \code{n} of variances related to the main effects
#'   \item \strong{varResidualEffect}: vector of length \code{n} of variances of residual effects
#'   \item \strong{listOption}: list of options used to obtained these results (obtained from \code{\link{QUALYPSO.check.option}})
#'   \item \strong{listScenarioInput}: list of scenario characteristics (obtained from \code{\link{QUALYPSO.process.scenario}})
#' }
#'
#' @export
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
QUALYPSO.ANOVA = function(phiStar,scenAvail,listOption=NULL){
  #########  process input #########
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

  # variance of individual effects
  contribEffect = list()
  for(i.eff in 1:nEff){
    contribEffect[[i.eff]] = matrix(nrow=nTypeEff[i.eff],ncol=y.POST)
    for(iIndEff in 1:nTypeEff[i.eff]){
      contribEffect[[i.eff]][iIndEff,] = apply(eff.POST[[i.eff]][,,iIndEff]^2,1,mean)/(varEffect[i.eff,]*nTypeEff[i.eff])
    }
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
  return(list(POSTERIOR=POST,QUANT=QUANT,MEAN=MEAN,
              varEffect=varEffect,contribEffect=contribEffect,
              varResidualEffect=varResidualEffect,
              listOption=listOption,listScenarioInput=listScenarioInput))

}



#==============================================================================
#' QUALYPSO
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach.
#'
#' @param Y matrix \code{nS} x \code{nY} or array \code{nG} x \code{nS} x \code{nY} of climate projections
#' @param scenAvail matrix of available combinations \code{nS} x \code{nEff}. The number of characteristics \code{nEff} corresponds to the
#' number of main effects which will be included in the ANOVA model.
#' @param vecYears (optional) vector of years corresponding to the projections (e.g. \code{vecYears=2001:2100}. Optional,
#' mainly used for records. By default, a vector \code{1:nY} is created.
#' @param indexReferenceYear (optional) index in \code{vecYears} corresponding to the control year. For example, if \code{vecYears=1980:2100}
#' and we want to specify a control year equals to 1990, we indicate \code{indexReferenceYear=11} or, equivalently \code{indexReferenceYear=which(vecYears==1990)}
#' \strong{if} \code{vecYears} is already available in the workspace
#' @param indexFutureYear index in \code{indexFutureYear} corresponding to a future year (similarly to \code{indexReferenceYear}). This index is necessary when
#' \code{Y} is an array \code{nG} x \code{nS} x \code{nY} available for \code{nG} grid points. Indeed, in this case, we run QUALYPSO only for one future year.
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{parSmooth}: smoothing parameter \code{spar} in \link[stats]{smooth.spline}: varies in [0,1]
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative)
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000). If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: number of kept samples (default: 2000). If \code{nKeep} is too small, MCMC samples might not be represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{nCluster}: number of clusters used for the parallelization (default: 1). When \code{nCluster} is greater than one, parallelization is used to
#'   apply \code{QUALYPSO} over multiple time steps or grid points simultaneously.
#'   \item \strong{doCompress}: logical, indicates if all the samples from the posterior distributions are stored (if FALSE)
#'   or if only quantiles are retrieved (if TRUE). Equals TRUE by default
#'   \item \strong{computeEmpEff}: vector of column indices in \code{scenAvail} corresponding to effects which are estimated empirically (e.g. interactions) when
#'   the number of available runs is not sufficient to identify / estimate these additional effects.
#'   \item \strong{quantileCompress}: vector of probabilities (in [0,1]) for which we compute the quantiles from the posterior distributions
#'    \code{quantileCompress = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)} by default
#' }
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{CLIMATEESPONSE}: list of climate change responses and corresponding internal variability. Contains \code{phiStar} (climate change responses),
#'   \code{etaStar} (deviation from the climate change responses as a result of internal variability), and \code{phi} (fitted climate responses)
#'   \item \strong{ANOVAPOST}: list of MCMC samples representing the posterior distributions of inferred quantities. \code{=NULL} if\code{listOption$doCompress=T}
#'   \item \strong{ANOVAQUANT}: list of quantiles from the posterior distributions of inferred quantities
#'   \item \strong{ANOVAMEAN}: list of mean of the posterior distributions of inferred quantities
#'   Each element contains the main effects (e.g. number of GCMs, RCMs, etc.)
#'   \item \strong{ANOVAVARIANCE}: matrix \code{nTypeEff} x \code{nY} of variances related to the main effects
#'   \item \strong{vecYears}: vector of years
#'   \item \strong{vecYearsANOVA}: vector of years for the ANOVA decomposition (start at \code{indexReferenceYear})
#'   \item \strong{Y}: matrix of available combinations given as inputs
#'   \item \strong{listOption}: list of options used to obtained these results (obtained from \code{\link{QUALYPSO.check.option}})
#'   \item \strong{listScenarioInput}: list of scenario characteristics (obtained from \code{\link{QUALYPSO.process.scenario}})
#' }
#'
#' @examples
#' ##########################################################################
#' # SYNTHETIC SCENARIOS
#' ##########################################################################
#' # create nS=3 fictive climate scenarios with 2 GCMs and 2 RCMs, for a period of nY=20 years
#' n=20
#' t=1:n/n
#'
#' # GCM effects (sums to 0 for each t)
#' effGCM1 = t*2
#' effGCM2 = t*-2
#'
#' # RCM effects (sums to 0 for each t)
#' effRCM1 = t*1
#' effRCM2 = t*-1
#'
#' # These climate scenarios are a sum of effects and a random gaussian noise
#' scenGCM1RCM1 = effGCM1 + effRCM1 + rnorm(n=n,sd=0.5)
#' scenGCM1RCM2 = effGCM1 + effRCM2 + rnorm(n=n,sd=0.5)
#' scenGCM2RCM1 = effGCM2 + effRCM1 + rnorm(n=n,sd=0.5)
#' Y = rbind(scenGCM1RCM1,scenGCM1RCM2,scenGCM2RCM1)
#'
#' # Here, scenAvail indicates that the first scenario is obtained with the combination of the
#' # GCM "GCM1" and RCM "RCM1", the second scenario is obtained with the combination of
#' # the GCM "GCM1" and RCM "RCM2" and the third scenario is obtained with the combination
#' # of the GCM "GCM2" and RCM "RCM1".
#' scenAvail = data.frame(GCM=c('GCM1','GCM1','GCM2'),RCM=c('RCM1','RCM2','RCM1'))
#'
#' ##########################################################################
#' # RUN QUALYPSO
#' ##########################################################################
#' # call main QUALYPSO function: two arguments are mandatory:
#' # - Y: Climate projections for nS scenarions and nY time steps. if Y is a matrix nS x nY, we
#' # run QUALYPSO nY times, for each time step. If Y is an array nG x nS x nY, for nG grid points,
#' # we run QUALYPSO nG times, for each grid point, for one time step specified using the argument
#' # indexFutureYear.
#' # - scenAvail: matrix or data.frame of available combinations nS x nEff. The number of
#' # characteristics nEff corresponds to the number of main effects which will be included in the
#' # ANOVA model. In the following example, we have nEff=2 main effects corresponding to the GCMs
#' # and RCMs.
#'
#' # Many options can be specified in the argument "listOption". Here, we change the default values
#' # for nBurn and nKeep in order to speed up computation time for this small example. However, it must
#' # be noticed that convergence and sampling of the posterior distributions often require higher
#' # values for these two parameters.
#' listOption = list(nBurn=100,nKeep=100,quantileCompress=c(0.025,0.5,0.975))
#'
#' # run QUALYPSO
#' QUALYPSOOUT = QUALYPSO(Y=Y, scenAvail=scenAvail, vecYears=2001:2020, listOption=listOption)
#'
#' ##########################################################################
#' # SOME PLOTS
#' ##########################################################################
#' # plot grand mean
#' plotQUALYPSOgrandmean(QUALYPSOOUT)
#'
#' # plot main GCM effects
#' plotQUALYPSOeffect(QUALYPSOOUT,iEff=1)
#'
#' # plot main RCM effects
#' plotQUALYPSOeffect(QUALYPSOOUT,iEff=2)
#'
#' # plot fraction of total variance for the differences sources of uncertainty
#' plotQUALYPSOTotalVarianceDecomposition(QUALYPSOOUT)
#'
#' # plot mean prediction and total variance with the differences sources of uncertainty
#' # for one scenario (e.g. a RCP scenario)
#' plotQUALYPSOTotalVarianceByScenario(QUALYPSOOUT,iEff=1,nameScenario='GCM1')
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. \url{https://doi.org/10.1175/JCLI-D-18-0606.1}.
#'
#' @export
#'
#' @author Guillaume Evin
QUALYPSO = function(Y,scenAvail,vecYears=NULL,indexReferenceYear=NULL,indexFutureYear=NULL,listOption=NULL){
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
    if(is.null(indexFutureYear)) stop("Since Y is an array of dimension 3: 'indexFutureYear' must be provided'")
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
    if(!any(indexReferenceYear==1:nY)){
      printMessageDimension(Y,scenAvail,vecYears)
      stop('wrong value for indexReferenceYear')
    }
  }else{
    indexReferenceYear = 1
  }

  # indexFutureYear
  if(paralType == 'Grid'){
    if(!any(indexFutureYear==1:nY)){
      printMessageDimension(Y,scenAvail,vecYears)
      stop('wrong value for indexFutureYear')
    }
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

    # if there are some runs which have only na values
    if(all(hasNa)){
      warning("Some runs leads to undefined variable changes (e.g. relative differences with a reference equals to 0)")
      iRunNA = which(!apply(phiStarAllTime,1,function(x) all(is.na(x))))
      if(any(iRunNA)){
        print('try to run with runs on well-defined variable changes')
        return(QUALYPSO(Y = Y[iRunNA,],scenAvail = scenAvail[iRunNA,],vecYears=vecYears,indexReferenceYear=indexReferenceYear,
                        indexFutureYear=indexFutureYear,listOption=listOption))
      }else{
        warning('Error in QUALYPSO: check NAs (or computation of relative diff. with zeros) in Y')
        return(NULL)
      }
    }

    # apply ANOVA after the reference year and when there is no na
    filterYearsAnova = 1:nY>indexReferenceYear&!hasNa
    if(!any(filterYearsAnova)){
      warning('Error in QUALYPSO: check indexReferenceYear and NAs in Y')
      return(NULL)
    }

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
    phiStar = t(phiStarAllTime[,,indexFutureYear])

    # vector of years for ANOVA
    vecYearsANOVA = vecYears[indexFutureYear]
  }


  #############################################
  # ANOVA on phiStar
  anova = QUALYPSO.ANOVA(phiStar = phiStar, scenAvail = scenAvail, listOption = listOption)


  #############################################
  # format variance
  ANOVAVARIANCE = list()
  ANOVAVARIANCE$eff = anova$varEffect
  ANOVAVARIANCE$contribEffect = anova$contribEffect
  ANOVAVARIANCE$ResidualEffect = anova$varResidualEffect
  ANOVAVARIANCE$InterVariability = varInterVariability
  ANOVAVARIANCE$TotalVar = Rfast::colsums(anova$varEffect) + anova$varResidualEffect + varInterVariability


  # return results
  return(list(CLIMATEESPONSE=list(phiStar=phiStarAllTime,etaStar=etaStarAllTime,phi=phiAllTime),
              ANOVAPOST=anova$POSTERIOR,ANOVAQUANT=anova$QUANT,ANOVAMEAN=anova$MEAN,paralType=paralType,
              ANOVAVARIANCE=ANOVAVARIANCE,vecYears=listOption$vecYears,vecYearsANOVA=vecYearsANOVA,
              Y=Y,listOption=anova$listOption,listScenarioInput=anova$listScenarioInput))
}



#==============================================================================
#' plotQUALYPSOgrandmean
#'
#' Plot prediction of grand mean ensemble. By default, we plot the credible interval corresponding to a probability 0.95.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param CIlevel probabilities for the credible intervals, default is equal to \code{c(0.025,0.975)}
#' @param lim y-axis limits (default is NULL)
#' @param col color for the overall mean and the credible interval
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOgrandmean = function(QUALYPSOOUT,CIlevel=c(0.025,0.975),lim=NULL,
                                   col='black',xlab="Years",ylab="Grand mean",addLegend=T,
                                   ...){
  # vector of years
  vecYears = QUALYPSOOUT$vecYearsANOVA

  # find index quantiles
  iMedian = which(QUALYPSOOUT$listOption$quantileCompress==0.5)
  iBinf = which(QUALYPSOOUT$listOption$quantileCompress==CIlevel[1])
  iBsup = which(QUALYPSOOUT$listOption$quantileCompress==CIlevel[2])
  if((length(iBinf)!=1)|(length(iBsup)!=1)) stop(paste0('Quantiles ',CIlevel, "are not available, check argument CIlevel"))

  # retrieve median and limits
  medRaw = QUALYPSOOUT$ANOVAQUANT$mu[,iMedian]
  binfRaw = QUALYPSOOUT$ANOVAQUANT$mu[,iBinf]
  bsupRaw = QUALYPSOOUT$ANOVAQUANT$mu[,iBsup]

  # get smooth limits
  med = predict(loess(medRaw~vecYears))
  binf = predict(loess(binfRaw~vecYears))
  bsup = predict(loess(bsupRaw~vecYears))

  # colors polygon
  colPoly = adjustcolor(col,alpha.f=0.2)

  # initiate plot
  if(is.null(lim)) lim = range(c(binf,bsup),na.rm=TRUE)
  plot(-100,-100,xlim=range(vecYears),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  # add confidence interval
  polygon(c(vecYears,rev(vecYears)),c(binf,rev(bsup)),col=colPoly,lty=0)

  # add median
  lines(vecYears,med,lwd=3,col=col)

  # legend
  if(addLegend){
    pctCI = plotQUALYPSOgetCI(QUALYPSOOUT,iBinf,iBsup)
    legend('topleft',bty='n',fill=c(NA,colPoly),lwd=c(2,NA),lty=c(1,NA),
           border=c(NA,col),col=c(col,NA),legend=c('Median',paste0(pctCI,'%CI')))
  }
}


#==============================================================================
#' plotQUALYPSOeffect
#'
#' Plot prediction of ANOVA effects for one main effect. By default, we plot we plot the credible intervals corresponding to a probability 0.95.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param iEff index of the main effect to be plotted in \code{QUALYPSOOUT$listScenarioInput$listEff}
#' @param includeMean if TRUE, the grand mean is added to the main effect in the plot
#' @param CIlevel probabilities for the credible intervals, default is equal to \code{c(0.025,0.975)}
#' @param lim y-axis limits (default is NULL)
#' @param col colors for each effect
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOeffect = function(QUALYPSOOUT,iEff,includeMean=FALSE,CIlevel=c(0.025,0.975),lim=NULL,
                                   col=1:20,xlab="Years",ylab="Effect",addLegend=TRUE,
                                   ...){
  # vector of years
  vecYears = QUALYPSOOUT$vecYearsANOVA

  # retrieve effects
  if(includeMean){
    QUANTEff = QUALYPSOOUT$ANOVAQUANT$ChangeByEffect[[iEff]]
  }else{
    QUANTEff = QUALYPSOOUT$ANOVAQUANT$effect[[iEff]]
  }
  nEff = dim(QUANTEff)[3]

  # find index quantiles
  iMedian = which(QUALYPSOOUT$listOption$quantileCompress==0.5)
  iBinf = which(QUALYPSOOUT$listOption$quantileCompress==CIlevel[1])
  iBsup = which(QUALYPSOOUT$listOption$quantileCompress==CIlevel[2])
  if((length(iBinf)!=1)|(length(iBsup)!=1)) stop(paste0('Quantiles ',CIlevel, "are not available, check argument CIlevel"))

  # retrieve median and limits
  medRaw = QUANTEff[,iMedian,]
  binfRaw = QUANTEff[,iBinf,]
  bsupRaw = QUANTEff[,iBsup,]

  # get smooth limits
  med = apply(medRaw,2,function(x) predict(loess(x~vecYears)))
  binf = apply(binfRaw,2,function(x) predict(loess(x~vecYears)))
  bsup = apply(bsupRaw,2,function(x) predict(loess(x~vecYears)))

  # initiate plot
  if(is.null(lim)) lim = range(c(binf,bsup),na.rm=TRUE)
  plot(-100,-100,xlim=range(vecYears),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  for(i in 1:nEff){
    # colors polygon
    colPoly = adjustcolor(col[i],alpha.f=0.2)

    # add confidence interval
    polygon(c(vecYears,rev(vecYears)),c(binf[,i],rev(bsup[,i])),col=colPoly,lty=0)

    # add median
    lines(vecYears,med[,i],lwd=3,col=col[i])
  }

  # legend
  if(addLegend){
    pctCI = plotQUALYPSOgetCI(QUALYPSOOUT,iBinf,iBsup)
    legend('topleft',bty='n',fill=c(NA,'black'),lwd=c(2,NA),lty=c(1,NA),
           border=c(NA,'black'),col=c('black',NA),legend=c('Median',paste0(pctCI,'%CI')))

    legend('bottomleft',bty='n',lwd=2,lty=1,col=col,
           legend=QUALYPSOOUT$listScenarioInput$listEff[[iEff]])
  }
}

#==============================================================================
# plotQUALYPSOgetCI
#
# return confidence level given indices in QUALYPSOOUT$listOption$quantileCompress
plotQUALYPSOgetCI = function(QUALYPSOOUT,iBinf,iBsup){
  vecq = QUALYPSOOUT$listOption$quantileCompress
  return(round((vecq[iBsup]-vecq[iBinf])*100))
}

#==============================================================================
# printMessageDimension
printMessageDimension = function(Y, scenAvail, vecYears){
  # Y
  dimY = dim(Y)
  if(length(dimY)==2){
    print(paste0('Y is a matrix of dimension nS=',dimY[1],' projections x nY=',dimY[2],' years'))
  }else if(length(dimY)==3){
    print(paste0('Y is an array of dimension nG=',dimY[1],' grid points x nS=',dimY[2],
                 'projections x nY=',dimY[3],' years'))
  }else(stop('Y must be a matrix nS x nY or a 3-dimension array  nG x nS x nY'))

  # scenAvail
  dimscenAvail = dim(scenAvail)
  if(length(dimscenAvail)==2){
    print(paste0('dimscenAvail is a matrix of dimension nS=',dimscenAvail[1],
                 ' projections x nEff=',dimscenAvail[2],' effects'))
  }else(stop('scenAvail must be a matrix nS x nEff'))

  # vecYears
  print("VecYears must be a vector of years or indices corresponding to the dimension 'time':")
  print(vecYears)
}


#==============================================================================
#' plotQUALYPSOTotalVarianceDecomposition
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param vecEff vector of indices corresponding to the main effects (NULL by default), so that the order of appearance in the plot can be modified
#' @param col colors for each source of uncertainty, the first two colors corresponding to internal variability and residual variability, respectively
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOTotalVarianceDecomposition = function(QUALYPSOOUT,vecEff=NULL,
                                                     col=c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1"),
                                                     xlab="Years",ylab="% Total Variance",addLegend=TRUE,...){
  # number of years
  vecYears = QUALYPSOOUT$vecYearsANOVA
  nY = length(vecYears)

  # Variance decomposition
  VV = QUALYPSOOUT$ANOVAVARIANCE


  # number of main effects
  if(is.null(vecEff)){
    vecEff = 1:nrow(VV$eff)
  }

  # select effects
  Veff = VV$eff[vecEff,]
  nEff = nrow(Veff)

  # concatenate variances
  Vbind = rbind(Veff,VV$ResidualEffect,VV$InterVariability)
  Vtot = colSums(Vbind)
  Vnorm = Vbind/t(replicate(n = nrow(Vbind), Vtot))



  # figure
  col = col[1:(nEff+2)]
  cum=cum.smooth=rep(0,nY)
  plot(-1,-1,xlim=range(vecYears),ylim=c(0,1),xaxs="i",yaxs="i",las=1,
       xlab=xlab,ylab=ylab,...)
  for(i in 1:(nEff+2)){
    cumPrevious = cum.smooth
    cum = cum + Vnorm[i,]
    cum.smooth = predict(loess(cum~vecYears))
    cum.smooth[cum.smooth<0] = 0
    polygon(c(vecYears,rev(vecYears)),c(cumPrevious,rev(cum.smooth)),col=rev(col)[i],lty=1)
  }
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    if(is.null(colnames(QUALYPSOOUT$listScenarioInput$scenAvail))){
      namesEff = paste0("Eff",vecEff)
    }else{
      namesEff = colnames(QUALYPSOOUT$listScenarioInput$scenAvail)[vecEff]
    }

    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}


#==============================================================================
#' plotQUALYPSOTotalVarianceByScenario
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param iEff index in \code{scenAvail} corresponding to the scenarios (e.g. RCP scenarios)
#' @param nameScenario name of the scenario to be plotted (as provided in \code{scenAvail})
#' @param probCI probability for the dredible interval, =0.9 by default
#' @param col colors for each source of uncertainty, the first two colors corresponding to internal variability and residual variability, respectively
#' @param ylim y-axis limits
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOTotalVarianceByScenario = function(QUALYPSOOUT,iEff,nameScenario,probCI=0.9,col=NULL,ylim=NULL,
                                               xlab="Years",ylab="Change variable",addLegend=TRUE,...){
  # number of years
  vecYears = QUALYPSOOUT$vecYearsANOVA
  nY = length(vecYears)

  # which scenario
  iScenario = which(QUALYPSOOUT$listScenarioInput$listEff[[iEff]] == nameScenario)

  # mean prediction
  meanPred = QUALYPSOOUT$ANOVAMEAN$ChangeByEffect[[iEff]][,iScenario]

  # Variance decomposition
  VV = QUALYPSOOUT$ANOVAVARIANCE

  # remove effect corresponding to the scenarios from the total variance
  Veff = VV$eff[-iEff,]

  # concatenate variances
  Vbind = rbind(Veff,VV$ResidualEffect,VV$InterVariability)
  nEff = nrow(Vbind)-2
  Vtot = colSums(Vbind)
  Vnorm = Vbind/t(replicate(n = nrow(Vbind), Vtot))

  # reverse
  vNormRev = apply(Vnorm,2,rev)


  # compute the lower bound if the distribution is gaussian
  binf = qnorm(p = (1-probCI)/2, mean = meanPred, sd = sqrt(Vtot))
  bsup = qnorm(p = 0.5+probCI/2, mean = meanPred, sd = sqrt(Vtot))

  # figure
  if(is.null(col)){
    default.col = c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1")
    col = default.col[1:(nEff+2)]
  }

  # obtain limits of the intervals, proportion corresponds to the part of the variance, lower and upper than the mean
  limIntInf = limIntSup = matrix(nrow=nEff+2,ncol=nY)
  limIntInf[1,] = predict(loess(binf~vecYears))
  limIntSup[1,] = predict(loess(bsup~vecYears))
  for(i in 1:(nEff+1)){
    binfi = limIntInf[i,]+vNormRev[i,]*(meanPred-binf)
    limIntInf[i+1,] = predict(loess(binfi~vecYears))
    bsupi = limIntSup[i,]-vNormRev[i,]*(bsup-meanPred)
    limIntSup[i+1,] = predict(loess(bsupi~vecYears))
  }

  # figure
  if(is.null(ylim)) ylim = c(min(binf),max(bsup))
  plot(-1,-1,xlim=range(vecYears),ylim=ylim,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1,...)
  for(i in 1:(nEff+2)){
    polygon(c(vecYears,rev(vecYears)),c(limIntInf[i,],rev(limIntSup[i,])),col=col[i],lty=1)
  }
  lines(vecYears,predict(loess(meanPred~vecYears)),col="white",lwd=1)

  # add horizontal lines
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    if(is.null(colnames(QUALYPSOOUT$listScenarioInput$scenAvail))){
      namesEff = paste0("Eff",1:nEff)
    }else{
      namesEff = colnames(QUALYPSOOUT$listScenarioInput$scenAvail)[-iEff]
    }

    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}
