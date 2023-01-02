###===============================###===============================###
### Guillaume Evin
### 30/06/2022, Grenoble
###  INRAE
### guillaume.evin@inrae.fr
###
### These functions use data augmentation and Bayesian techniques for the assessment
### of single-member and incomplete ensembles of climate projections.
### It provides unbiased estimates of climate change responses of all
### simulation chains and of all uncertainty variables. It additionally propagates
### uncertainty due to missing information in the estimates.
###
### REFERENCES
###
### references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
### Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
### Journal of Climate. J. Climate, 32, 2423–2440. <doi:10.1175/JCLI-D-18-0606.1>.
###
###
###===============================###===============================###


#==============================================================================
#' get.Qstar.mat
#'
#' Provide matrix containing Helmert contrasts (see Eq. A7 in Evin et al., 2019).
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020) <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @param p integer
#'
#' @return \item{matrix}{p x (p-1) matrix containing Helmert contrasts}
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. <doi:10.1175/JCLI-D-18-0606.1>.
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
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. <doi:10.1175/JCLI-D-18-0606.1>.
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
#' @param spar smoothing parameter \code{spar} in \code{\link[stats]{smooth.spline}}: varies in [0,1]
#' @param Xmat matrix of predictors corresponding to the projections, e.g. time or global temperature.
#' @param Xfut values of the predictor over which the ANOVA will be applied.
#' @param typeChangeVariable type of change variable: "abs" (absolute, value by default) or "rel" (relative)
#'
#' @return list with the following fields for each simulation chain:
#' \itemize{
#'   \item \strong{phiStar}: \code{nS x nF}, climate change responses
#'   \item \strong{etaStar}: \code{nS x nY}, deviation from the climate change response
#'   due to the internal variability, for \code{Xmat}
#'   \item \strong{phi}: \code{nS x nF}, raw trends obtained using \link[stats]{smooth.spline}
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
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. J. Climate, 32, 2423–2440. <doi:10.1175/JCLI-D-18-0606.1>.
fit.climate.response = function(Y, spar, Xmat, Xfut, typeChangeVariable){

  # number of simulation chains
  nS = nrow(Y)

  # length of the simulation chains
  nY = ncol(Y)

  # number of future time/global.tas
  nFut = length(Xfut)

  # Xref is the reference value for X used to compute absolute or relative changes
  # usually indicating the reference (control) period of the current climate
  # For simplicity, we consider that it is the first Xfut value (i.e. first year
  # in Xfut or minimum global temperature)
  Xref = Xfut[1]

  # prepare outputs
  phiStar = phi = matrix(nrow=nS,ncol=nFut)
  etaStar = YStar = matrix(nrow=nS,ncol=nY)
  climateResponse = list()

  for(iS in 1:nS){
    # projection for this simulation chain
    Ys = Y[iS,]
    Xs = Xmat[iS,]

    # fit a smooth signal (smooth cubic splines)
    zz = !is.na(Ys)
    smooth.spline.out<-stats::smooth.spline(Xs[zz],Ys[zz],spar=spar)

    # store spline object
    climateResponse[[iS]] = smooth.spline.out

    # fitted responses at the points of the fit (for etaStar)
    phiY = predict(smooth.spline.out, Xs)$y

    # fitted responses at unknown points ("Xfut")
    phiS = predict(smooth.spline.out, Xfut)$y

    # climate response of the reference/control time/global tas
    phiC = predict(smooth.spline.out, Xref)$y

    # store climate response for this simulation chain
    phi[iS,] = phiS

    # Climate change response: phiStar, and internal variability expressed as a change: etaStar
    if(typeChangeVariable=='abs'){
      # Eq. 5
      phiStar[iS,] = phiS-phiC
      etaStar[iS,] = Ys-phiY
      YStar[iS,] = Ys-phiC
    }else if(typeChangeVariable=='rel'){
      # Eq. 6
      phiStar[iS,] = phiS/phiC-1
      etaStar[iS,] = (Ys-phiY)/phiC
      YStar[iS,] = (Ys-phiC)/phiC
    }else{
      stop("fit.climate.response: argument type.change.var must be equal to 'abs' (absolute changes) or 'rel' (relative changes)")
    }
  }

  # Variance related to the internal variability: considered constant over the time period
  # (see Eq. 22 and 23 in Hingray and Said, 2014). We use a direct empirical estimate
  # of the variance of eta for each simulation chain and take the mean, see Eq. 19
  varInterVariability = mean(apply(etaStar,2,var,na.rm=T),na.rm=T)

  # return objects
  return(list(phiStar=phiStar,etaStar=etaStar,YStar=YStar,phi=phi,
              climateResponse=climateResponse,varInterVariability=varInterVariability))
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
#'   Each element is a matrix \code{nMCMC} x \code{nMaineff}, and \code{nMaineff} is the number of main effects (e.g. number of GCMs, RCMs, etc.)
#' }
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
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
  }


  # return MCM samples from the posterior
  return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST))
}

#==============================================================================
#' QUALYPSO.process.scenario
#'
#' Process input scenarios.
#'
#' @param scenAvail data.frame \code{nS} x \code{nEff} with the \code{nEff} characteristics (e.g. type of GCM) for each of the \code{nS} x \code{nS} scenarios
#'
#' @return list of preprocessed objects (\code{listEff, scenAvail, scenComp, nEff, nTypeEff, nComp, isMissing, nMissing, iMatchScen,
#' indexEffInCompScen, Qmat})
#'
#' @author Guillaume Evin
QUALYPSO.process.scenario = function(scenAvail){
  # number of scenarios
  nS = nrow(scenAvail)

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
              Qmat=Qmat))
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

  # spar
  if('spar' %in% names(listOption)){
    spar = listOption[['spar']]
    if(!(is.numeric(spar)&(spar>0&spar<=10))) stop('spar must be in ]0,10]')
  }else{
    listOption[['spar']] = 1
  }

  # ANOVAmethod
  if('ANOVAmethod' %in% names(listOption)){
    ANOVAmethod = listOption[['ANOVAmethod']]
    if(!(ANOVAmethod%in%c('QUALYPSO','lm'))) stop("ANOVAmethod must be equal to 'QUALYPSO' or 'lm'")
  }else{
    listOption[['ANOVAmethod']] = 'lm'
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

  # probCI
  if('probCI' %in% names(listOption)){
    probCI = listOption[['probCI']]
    if(!(is.numeric(probCI))) stop('wrong value for probCI')
    if(probCI<0|probCI>1) stop('wrong value for probCI: outside [0,1]')
  }else{
    listOption[['probCI']] = 0.9
  }

  # quantilePosterior
  if('quantilePosterior' %in% names(listOption)){
    quantilePosterior = listOption[['quantilePosterior']]
    if(!(is.numeric(quantilePosterior))) stop('wrong value for quantilePosterior')
  }else{
    listOption[['quantilePosterior']] = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)
  }

  # Version
  listOption[['version']] = 'v1.2'

  return(listOption)
}


#==============================================================================
#' QUALYPSO.ANOVA
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach.
#'
#' @param phiStar matrix of climate change responses (absolute or relative changes): \code{nS} x \code{n}.
#' \code{n} can be the number of time steps or the number of grid points
#' @param scenAvail data.frame \code{nS} x \code{nEff} with the \code{nEff} characteristics (e.g. type of GCM) for each of the \code{nS} x \code{nS} scenarios
#' @param listOption list of options (see \code{\link{QUALYPSO}})
#' @param namesEff names of the main effects
#'
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{GRANDMEAN}: List of estimates for the grand mean:
#'   \itemize{
#'      \item {strong}: MEAN: vector of length \code{n} of posterior means
#'      \item {strong}: SD: vector of length \code{n} of posterior standard dev.
#'      \item {strong}: CI: matrix \code{n} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'      \item {strong}: QUANT: matrix \code{n} x \code{nQ} of quantiles related
#'      to the probabilities \code{quantilePosterior} given in \code{listOption}
#'   }
#'   \item \strong{RESIDUALVAR}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item {strong}: MEAN: vector of length \code{n} of posterior means
#'      \item {strong}: SD: vector of length \code{n} of posterior standard dev.
#'      \item {strong}: CI: matrix \code{n} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'      \item {strong}: QUANT: matrix \code{n} x \code{nQ} of quantiles related
#'      to the probabilities \code{quantilePosterior} given in \code{listOption}
#'   }
#'   \item \strong{MAINEFFECT}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item {strong}: MEAN: matrix \code{n} x \code{nTypeEff} of posterior means
#'      \item {strong}: SD: matrix \code{n} x \code{nTypeEff} of posterior standard dev.
#'      \item {strong}: CI: array \code{n} x 2 x \code{nTypeEff} of credible
#'      intervals of probability \code{probCI} given in \code{listOption}.
#'      \item {strong}: QUANT: array \code{n} x \code{nQ} x \code{nTypeEff} of
#'      quantiles related to the probabilities \code{quantilePosterior} given in
#'      \code{listOption}
#'   }
#'   \item \strong{CHANGEBYEFFECT}: For each main effect, list of estimates for
#'   the mean change by main effect, i.e. mean change by scenario (RCP4.5). For
#'   each main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item {strong}: MEAN: matrix \code{n} x \code{nTypeEff} of posterior means
#'      \item {strong}: SD: matrix \code{n} x \code{nTypeEff} of posterior standard dev.
#'      \item {strong}: CI: array \code{n} x 2 x \code{nTypeEff} of credible
#'      intervals of probability \code{probCI} given in \code{listOption}.
#'      \item {strong}: QUANT: array \code{n} x \code{nQ} x \code{nTypeEff} of
#'      quantiles related to the probabilities \code{quantilePosterior} given in
#'      \code{listOption}
#'   }
#'   \item \strong{EFFECTVAR}: variability related to the main effects (i.e.
#'   variability between the different RCMs, GCMs,..). Matrix \code{n} x
#'   \code{nTypeEff}
#'   \item \strong{CONTRIB_EACH_EFFECT}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{n} x \code{nTypeEff}
#'   \item \strong{listOption}: list of options used to obtained these results
#'    (obtained from \code{\link{QUALYPSO.check.option}})
#'   \item \strong{listScenarioInput}: list of scenario characteristics
#'    (obtained from \code{\link{QUALYPSO.process.scenario}})
#' }
#'
#' @export
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
QUALYPSO.ANOVA = function(phiStar,scenAvail,listOption=NULL,namesEff){
  #########  process input #########
  # number of grid points / years
  n = dim(phiStar)[2]

  # Check list of options
  listOption = QUALYPSO.check.option(listOption)

  # number of MCMC samples
  nKeep = listOption$nKeep
  vec.keep = (listOption$nBurn+1):listOption$nMCMC

  # Process scenarios data.frame to get different objects
  listScenarioInput = QUALYPSO.process.scenario(scenAvail = scenAvail)
  nEff = listScenarioInput$nEff
  nTypeEff = listScenarioInput$nTypeEff

  #########  Apply ANOVA for each time step: parallel computation over time steps #########
  # apply parallel computation only if more than one cluster has been indicated
  if(listOption$nCluster==1){
    Anova.POST = list()
    for(i in 1:n){
      Anova.POST[[i]] = QUALYPSO.ANOVA.i(phiStar.i=phiStar[,i], nMCMC=listOption$nMCMC,listScenarioInput = listScenarioInput)
    }
  }else{
    i = NULL # avoid warning during check
    Anova.POST <- foreach(i=1:n, .multicombine=TRUE,
                          .export=c("QUALYPSO.ANOVA.i")) %dopar% {
                            POST.out = QUALYPSO.ANOVA.i(phiStar.i=phiStar[,i], nMCMC=listOption$nMCMC,listScenarioInput = listScenarioInput)
                            return(POST.out)
                          }
  }


  #-------------------------------------------------------------------------
  # POSTPROCESS: EXTRACT SOME QUANTITIES
  #-------------------------------------------------------------------------
  qCI = c((1-listOption$probCI)/2,0.5+listOption$probCI/2)
  qPost = listOption$quantilePosterior


  #============================
  # GRANDMEAN
  #============================
  GRANDMEAN = list()

  # extract posterior
  mu.POST = matrix(nrow=n,ncol=nKeep)
  for(i in 1:n) mu.POST[i,]=Anova.POST[[i]]$mu[vec.keep]

  # mean
  GRANDMEAN$MEAN = apply(mu.POST,1,mean)

  # sd
  GRANDMEAN$SD = apply(mu.POST,1,sd)

  # CI
  GRANDMEAN$CI = t(apply(mu.POST,1,quantile,probs=qCI))

  # QUANT
  GRANDMEAN$QUANT = t(apply(mu.POST,1,quantile,probs=qPost))


  #============================
  # RESIDUALVAR
  #============================
  # variance of residual effects, With the Bayesian estimation, this variance
  # is part of the inference (sigma^2) and we take the mean of the posterior as the point estimate
  RESIDUALVAR = list()

  # extract posterior
  sigma2.POST = matrix(nrow=n,ncol=nKeep)
  for(i in 1:n) sigma2.POST[i,]=Anova.POST[[i]]$sigma2[vec.keep]

  # mean
  RESIDUALVAR$MEAN = apply(sigma2.POST,1,mean)

  # sd
  RESIDUALVAR$SD = apply(sigma2.POST,1,sd)

  # CI
  RESIDUALVAR$CI = t(apply(sigma2.POST,1,quantile,probs=qCI))

  # QUANT
  RESIDUALVAR$QUANT = t(apply(sigma2.POST,1,quantile,probs=qPost))


  #==============================================================
  # MAINEFFECT + CHANGEBYEFFECT + EFFECTVAR + CONTRIB_EACH_EFFECT
  #==============================================================
  # main effects
  MAINEFFECT = list()

  # change by effect
  CHANGEBYEFFECT = list()

  # variance of the main effects
  EFFECTVAR = matrix(nrow=n,ncol=nEff)

  # contribution to the variance of individual effects
  CONTRIB_EACH_EFFECT = list()

  # main effects
  for(i.eff in 1:nEff){
    eff = namesEff[i.eff]
    #============
    # MAINEFFECT
    #============
    # trim posterior distributions of the main effects
    eff.POST = array(dim=c(n,nKeep,nTypeEff[i.eff]))
    for(i in 1:n) eff.POST[i,,]=Anova.POST[[i]]$effect[[i.eff]][vec.keep,]

    # retrieve estimates for the main effects
    lHatMain =  list()
    lHatMain$MEAN = apply(eff.POST,c(1,3),mean)
    lHatMain$SD = apply(eff.POST,c(1,3),sd)
    eff.ci = apply(eff.POST,c(1,3),quantile,probs=qCI)
    lHatMain$CI = aperm(eff.ci, c(2,1,3))
    eff.quant = apply(eff.POST,c(1,3),quantile,probs=qPost)
    lHatMain$QUANT = aperm(eff.quant, c(2,1,3))
    MAINEFFECT[[eff]] = lHatMain

    #===============
    # CHANGEBYEFFECT
    #===============
    # trim posterior distributions of the main effects
    meanChange.POST = array(dim=c(n,nKeep,nTypeEff[i.eff]))
    for(j in 1:nTypeEff[i.eff]) meanChange.POST[,,j]=eff.POST[,,j] + mu.POST

    # retrieve estimates for the mean change by effect
    lHatchange =  list()
    lHatchange$MEAN = apply(meanChange.POST,c(1,3),mean)
    lHatchange$SD = apply(meanChange.POST,c(1,3),sd)
    eff.ci = apply(meanChange.POST,c(1,3),quantile,probs=qCI)
    lHatchange$CI = aperm(eff.ci, c(2,1,3))
    eff.quant = apply(meanChange.POST,c(1,3),quantile,probs=qPost)
    lHatchange$QUANT = aperm(eff.quant, c(2,1,3))
    CHANGEBYEFFECT[[eff]] = lHatchange

    #============
    # EFFECTVAR
    #============
    # predictive variance: mean of variances, which correspond to the mean of squares since the mean is 0 by constraint (see Eq 16, 17 and 18)
    EFFECTVAR[,i.eff] = apply(eff.POST^2,1,mean)

    #============
    # CONTRIB_EACH_EFFECT
    #============
    matContrib = matrix(nrow=n,ncol=nTypeEff[i.eff])
    for(iIndEff in 1:nTypeEff[i.eff]){
      # individual effect squared
      indEff2 = apply(eff.POST[,,iIndEff]^2,1,mean)
      # variance of this main effect x number of factors
      infEff2Tot = (EFFECTVAR[,i.eff]*nTypeEff[i.eff])
      # contribution of this individual effect, as a percentage
      matContrib[,iIndEff] = indEff2/infEff2Tot
    }
    CONTRIB_EACH_EFFECT[[eff]] = matContrib
  }


  # return results
  return(list(GRANDMEAN=GRANDMEAN,
              RESIDUALVAR=RESIDUALVAR,
              MAINEFFECT=MAINEFFECT,
              CHANGEBYEFFECT=CHANGEBYEFFECT,
              EFFECTVAR=EFFECTVAR,
              CONTRIB_EACH_EFFECT=CONTRIB_EACH_EFFECT,
              listOption=listOption,
              listScenarioInput=listScenarioInput))
}


#==============================================================================
#' lm.ANOVA
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach.
#'
#' @param phiStar matrix of climate change responses (absolute or relative changes): \code{nS} x \code{n}. \code{n} can be the number of time steps or the number of grid points
#' @param scenAvail data.frame \code{nS} x \code{nEff} with the \code{nEff} characteristics (e.g. type of GCM) for each of the \code{nS} x \code{nS} scenarios
#' @param listOption list of options (see \code{\link{QUALYPSO}})
#' @param namesEff names of the main effects
#'
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{GRANDMEAN}: List of estimates for the grand mean:
#'   \itemize{
#'      \item {strong}: MEAN: vector of length \code{n} of means
#'      \item {strong}: SD: vector of length \code{n} of standard dev.
#'      \item {strong}: CI: matrix \code{n} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'   }
#'   \item \strong{RESIDUALVAR}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item {strong}: MEAN: vector of length \code{n}
#'   }
#'   \item \strong{MAINEFFECT}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item {strong}: MEAN: matrix \code{n} x \code{nTypeEff}
#'   }
#'   \item \strong{CHANGEBYEFFECT}: For each main effect, list of estimates for
#'   the mean change by main effect, i.e. mean change by scenario (RCP4.5). For
#'   each main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item {strong}: MEAN: matrix \code{n} x \code{nTypeEff}
#'   }
#'   \item \strong{EFFECTVAR}: variability related to the main effects (i.e.
#'   variability between the different RCMs, GCMs,..). Matrix \code{n} x
#'   \code{nTypeEff}
#'   \item \strong{CONTRIB_EACH_EFFECT}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{n} x \code{nTypeEff}
#'   \item \strong{listOption}: list of options used to obtained these results
#'    (obtained from \code{\link{QUALYPSO.check.option}})
#'   \item \strong{listScenarioInput}: list of scenario characteristics
#'    (obtained from \code{\link{QUALYPSO.process.scenario}})
#' }
#'
#' @export
#'
#' @author Guillaume Evin
lm.ANOVA = function(phiStar,scenAvail,listOption=NULL,namesEff){
  #########  process input #########
  # number of grid points / years
  n = dim(phiStar)[2]

  # Check list of options
  listOption = QUALYPSO.check.option(listOption)

  # Process scenarios data.frame to get different objects
  listScenarioInput = QUALYPSO.process.scenario(scenAvail = scenAvail)
  nEff = listScenarioInput$nEff
  nTypeEff = listScenarioInput$nTypeEff
  listEff = listScenarioInput$listEff

  #########  Apply lm for each time step #########
  # contrasts
  list.contrasts = list()
  for(j in 1:length(namesEff)){
    list.contrasts[[namesEff[j]]] = contr.sum
  }

  # formula
  formula = paste0("phiStar ~ ",paste0(namesEff,collapse = " + "))

  # build data for the call to the lm function
  lm.data = scenAvail

  lm.out = lm.sum = list()
  for(i in 1:n){
    lm.data$phiStar=phiStar[,i]
    lm.out[[i]] = lm(formula, lm.data,contrasts=list.contrasts)
    if(any(is.na(lm.out[[i]]$coefficients))){
      stop("Undefined effects probably due to an ill-posed ANOVA")
    }
    lm.sum[[i]] = summary(lm.out[[i]])
  }

  # list of effects
  listEff.LM = lm.out[[1]]$xlevels


  #-------------------------------------------------------------------------
  # POSTPROCESS: EXTRACT SOME QUANTITIES
  #-------------------------------------------------------------------------
  pCI = c((1-listOption$probCI)/2,0.5+listOption$probCI/2)
  qPost = listOption$quantilePosterior
  Qt <- qt(pCI, lm.out[[1]]$df.residual, lower.tail = TRUE)

  #============================
  # GRANDMEAN
  #============================
  GRANDMEAN = list()
  GRANDMEAN$MEAN = rep(0,n)
  GRANDMEAN$SD = rep(0,n)
  GRANDMEAN$CI = matrix(data=0,nrow = n, ncol = 2)

  # extract estimates for each time/grid point
  for(i in 1:n){
    lm.sum.i = lm.sum[[i]]
    GRANDMEAN$MEAN[i] = lm.sum.i$coefficients[1,1]
    GRANDMEAN$SD[i] = lm.sum.i$coefficients[1,2]
    GRANDMEAN$CI[i,] = GRANDMEAN$MEAN[i] + outer(GRANDMEAN$SD[i], Qt)
  }


  #============================
  # RESIDUALVAR
  #============================
  # variance of residual effects, With the Bayesian estimation, this variance
  # is part of the inference (sigma^2) and we take the mean of the posterior as the point estimate
  RESIDUALVAR = list()
  RESIDUALVAR$MEAN = rep(0,n)

  # extract estimates for each time/grid point
  for(i in 1:n){
    lm.sum.i = lm.sum[[i]]
    RESIDUALVAR$MEAN[i] = lm.sum.i$sigma^2
  }


  #============================
  # MAINEFFECT + CHANGEBYEFFECT
  #============================
  # main effects
  MAINEFFECT = list()

  # change by effect
  CHANGEBYEFFECT = list()

  # variance of the main effects
  EFFECTVAR = matrix(data=0,nrow=n,ncol=nEff)

  # contribution to the variance of individual effects
  CONTRIB_EACH_EFFECT = list()

  # main effects
  for(i.eff in 1:nEff){
    eff = namesEff[i.eff]
    #============
    # MAINEFFECT
    #============
    # retrieve estimates for the main effects
    lHat = list()
    lHat$MEAN = matrix(data = 0, nrow = n, ncol = nTypeEff[i.eff])
    for(i in 1:n){
      lm.out.i = lm.out[[i]]
      matchEff = match(listEff[[i.eff]],listEff.LM[[eff]])
      e.raw = lm.out.i$coefficients[lm.out.i$assign==i.eff]
      e.unsorted = c(e.raw,-sum(e.raw))
      lHat$MEAN[i,] = e.unsorted[matchEff]
    }
    MAINEFFECT[[eff]] = lHat


    #===============
    # CHANGEBYEFFECT
    #===============
    changeHat = list()
    changeHat$MEAN = MAINEFFECT[[eff]]$MEAN + replicate(nTypeEff[i.eff],GRANDMEAN$MEAN)
    CHANGEBYEFFECT[[eff]] = changeHat

    #============
    # EFFECTVAR
    #============
    # predictive variance: mean of variances, which correspond to the mean of squares since the mean is 0 by constraint (see Eq 16, 17 and 18)
    EFFECTVAR[,i.eff] = apply(MAINEFFECT[[eff]]$MEAN^2,1,mean)

    #============
    # CONTRIB_EACH_EFFECT
    #============
    matContrib = matrix(data=0,nrow=n,ncol=nTypeEff[i.eff])
    AllEff2 = MAINEFFECT[[eff]]$MEAN^2
    for(iIndEff in 1:nTypeEff[i.eff]){
      # individual effect squared
      indEff2 = AllEff2[,iIndEff]
      # variance of this main effect x number of factors
      infEff2Tot = (EFFECTVAR[,i.eff]*nTypeEff[i.eff])
      # contribution of this individual effect, as a percentage
      matContrib[,iIndEff] = indEff2/infEff2Tot
    }
    CONTRIB_EACH_EFFECT[[eff]] = matContrib
  }


  # return results
  list.output = list(GRANDMEAN=GRANDMEAN,
                     RESIDUALVAR=RESIDUALVAR,
                     MAINEFFECT=MAINEFFECT,
                     CHANGEBYEFFECT=CHANGEBYEFFECT,
                     EFFECTVAR=EFFECTVAR,
                     CONTRIB_EACH_EFFECT=CONTRIB_EACH_EFFECT,
                     listOption=listOption,
                     listScenarioInput=listScenarioInput)
}



#==============================================================================
#' QUALYPSO
#'
#' Partition uncertainty in climate responses using an ANOVA applied to climate change responses.
#'
#' @param Y matrix \code{nS} x \code{nY} or array \code{nG} x \code{nS} x \code{nY} of climate projections.
#' @param scenAvail data.frame \code{nS} x \code{nEff} with the \code{nEff} characteristics
#' (e.g. type of GCM) for each of the \code{nS} scenarios. The number of characteristics
#'  \code{nEff} corresponds to the number of main effects that will be included in the ANOVA model.
#' @param X (optional) predictors corresponding to the projections, e.g. time or global temperature.
#' It can be a vector if the predictor is the same for all scenarios (e.g. \code{X=2001:2100}) or
#' a matrix of the same size as Y if these predictors are different for the scenarios. By default,
#' a vector \code{1:nY} is created.
#' @param Xfut (optional) \code{nF} values of the predictor over which the ANOVA will be applied. It must be
#' a vector of values within the range of values of X. By default, it corresponds to X if X is a vector,
#' \code{1:nY} if X is \code{NULL} or a vector of 10 values equally spaced between the minimum and
#' maximum values of X if X is a matrix.
#' @param iFut index in \code{1:nF} corresponding to a future predictor value . This index is
#' necessary when \code{Y} is an array \code{nG} x \code{nS} x \code{nY} available for \code{nG} grid points.
#' Indeed, in this case, we run QUALYPSO only for one future predictor.
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{spar}: smoothing parameter \code{spar} in \link[stats]{smooth.spline}: typically (but not necessarily) in (0,1].
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative).
#'   \item \strong{ANOVAmethod}: ANOVA method: "QUALYPSO" applies the method described in Evin et al. (2020),
#'   "lm" applies a simple linear model to estimate the main effects.
#'   \item \strong{nBurn}: if \code{ANOVAmethod=="QUALYPSO"}, number of burn-in samples (default: 1000).
#'   If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: if \code{ANOVAmethod=="QUALYPSO"}, number of kept samples (default: 2000).
#'   If \code{nKeep} is too small, MCMC samples might not represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{nCluster}: number of clusters used for the parallelization (default: 1).
#'   When \code{nCluster} is greater than one, parallelization is used to
#'   apply \code{QUALYPSO} over multiple time steps or grid points simultaneously.
#'   \item \strong{probCI}: probability (in [0,1]) for the confidence intervals, \code{probCI = 0.9} by default.
#'   \item \strong{quantilePosterior}: vector of probabilities (in [0,1]) for which
#'   we compute the quantiles from the posterior distributions
#'    \code{quantilePosterior = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)} by default.
#'   \item \strong{climResponse}: NULL by default. If it is provided, it must correspond to the outputs
#'   of \code{\link{fit.climate.response}}, i.e. a list with \code{phiStar} [nS x nF], \code{etaStar} [nS x nY],
#'    \code{phi} [nS x nF] and \code{varInterVariability} [scalar] if \code{Y} is a matrix [nS x nY],
#'    or a list with \code{phiStar} [nG x nS x nF], \code{etaStar} [nG x nS x nY], \code{phi} [nG x nS x nF] and
#'     \code{varInterVariability} vector of length \code{nG} if \code{Y} is an array [nG x nS x nY].
#' }
#'
#' @return  List providing the results for each of the \code{n} values of \code{Xfut}
#' if \code{Y} is a matrix or for each grid point if \code{Y} is an array, with the following fields:
#' \itemize{
#'   \item \strong{CLIMATEESPONSE}: list of climate change responses and
#'  corresponding internal variability. Contains \code{phiStar} (climate change
#'  responses), \code{etaStar} (deviation from the climate change responses as
#'  a result of internal variability), \code{Ystar} (change variable from the
#'  projections),and \code{phi} (fitted climate responses).
#'   \item \strong{GRANDMEAN}: List of estimates for the grand mean:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{n} of means.
#'      \item \strong{SD}: vector of length \code{n} of standard dev.
#'      if \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{CI}: matrix \code{n} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{QUANT}: matrix \code{n} x \code{nQ} of quantiles of
#'      probability \code{quantilePosterior} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'   }
#'   \item \strong{MAINEFFECT}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item \strong{MEAN}: matrix \code{n} x \code{nTypeEff}
#'      \item \strong{SD}: matrix \code{n} x \code{nTypeEff} of standard dev.
#'      if \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{CI}: array \code{n} x 2 x \code{nTypeEff} of credible
#'      intervals of probability \code{probCI} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{QUANT}: array \code{n} x \code{nQ} x \code{nTypeEff} of
#'      quantiles of probability \code{quantilePosterior} given in
#'      \code{listOption} if \code{ANOVAmethod=="QUALYPSO"}.
#'   }
#'   \item \strong{CHANGEBYEFFECT}: For each main effect, list of estimates for
#'   the mean change by main effect, i.e. mean change by scenario. For
#'   each main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item \strong{MEAN}: matrix \code{n} x \code{nTypeEff}
#'      \item \strong{SD}: matrix \code{n} x \code{nTypeEff} of standard dev.
#'      if \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{CI}: array \code{n} x 2 x \code{nTypeEff} of credible
#'      intervals of probability \code{probCI} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{QUANT}: array \code{n} x \code{nQ} x \code{nTypeEff} of
#'      quantiles of probability \code{quantilePosterior} given in
#'      \code{listOption} if \code{ANOVAmethod=="QUALYPSO"}.
#'   }
#'   \item \strong{EFFECTVAR}: Matrix \code{n} x \code{nTypeEff} giving, for each
#'   time variability related to the main effects (i.e.
#'   variability between the different RCMs, GCMs,..).
#'   \item \strong{CONTRIB_EACH_EFFECT}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{n} x \code{nTypeEff}
#'   \item \strong{RESIDUALVAR}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{n}.
#'      \item \strong{SD}: vector of length \code{n} of standard dev.
#'      if \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{CI}: matrix \code{n} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'      \item \strong{QUANT}: matrix \code{n} x \code{nQ} of quantiles of
#'      probability \code{quantilePosterior} given in \code{listOption} if
#'      \code{ANOVAmethod=="QUALYPSO"}.
#'   }
#'   \item \strong{INTERNALVAR}: Internal variability (constant over time)
#'   \item \strong{TOTALVAR}: total variability, i.e. the sum of internal variability,
#'       residual variability and variability related to the main effects
#'   \item \strong{DECOMPVAR}: Decomposition of the total variability for each component
#'   \item \strong{RESERR}: differences between the climate change responses and the additive anova formula (grand mean + main effects)
#'   \item \strong{Xmat}: matrix of predictors
#'   \item \strong{Xfut}: future predictor values
#'   \item \strong{paralType}: type of parallelisation (Time or Grid)
#'   \item \strong{namesEff}: names of the main effects
#'   \item \strong{Y}: matrix of available combinations given as inputs
#'   \item \strong{listOption}: list of options used to obtained these results
#'   (obtained from \code{\link{QUALYPSO.check.option}})
#'   \item \strong{listScenarioInput}: list of scenario characteristics
#'   (obtained from \code{\link{QUALYPSO.process.scenario}})
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
#' Y.synth = rbind(scenGCM1RCM1,scenGCM1RCM2,scenGCM2RCM1)
#'
#' # Here, scenAvail indicates that the first scenario is obtained with the combination of the
#' # GCM "GCM1" and RCM "RCM1", the second scenario is obtained with the combination of
#' # the GCM "GCM1" and RCM "RCM2" and the third scenario is obtained with the combination
#' # of the GCM "GCM2" and RCM "RCM1".
#' scenAvail.synth = data.frame(GCM=c('GCM1','GCM1','GCM2'),RCM=c('RCM1','RCM2','RCM1'))
#'
#' ##########################################################################
#' # RUN QUALYPSO
#' ##########################################################################
#' # call main QUALYPSO function: two arguments are mandatory:
#' # - Y: Climate projections for nS scenarios and nY time steps. if Y is a matrix nS x nY, we
#' # run QUALYPSO nY times, for each time step. If Y is an array nG x nS x nY, for nG grid points,
#' # we run QUALYPSO nG times, for each grid point, for one time step specified using the argument
#' # iFut
#' # - scenAvail: matrix or data.frame of available combinations nS x nEff. The number of
#' # characteristics nEff corresponds to the number of main effects that will be included in the
#' # ANOVA model. In the following example, we have nEff=2 main effects corresponding to the GCMs
#' # and RCMs.
#'
#' # Many options can be specified in the argument "listOption". When ANOVAmethod=="QUALYPSO"
#' # a Bayesian inference is performed. Here, we change the default values for nBurn and nKeep
#' # in order to speed up computation time for this small example. However, it must be noticed
#' # that convergence and sampling of the posterior distributions often require higher values
#' #  for these two arguments.
#' listOption = list(nBurn=100,nKeep=100,ANOVAmethod="QUALYPSO",quantilePosterior=c(0.025,0.5,0.975))
#'
#' # run QUALYPSO
#' QUALYPSO.synth = QUALYPSO(Y=Y.synth, scenAvail=scenAvail.synth, X=2001:2020, listOption=listOption)
#'
#' ##########################################################################
#' # SOME PLOTS
#' ##########################################################################
#' # plot grand mean
#' plotQUALYPSOgrandmean(QUALYPSO.synth,xlab="Years")
#'
#' # plot main GCM effects
#' plotQUALYPSOeffect(QUALYPSO.synth,nameEff="GCM",xlab="Years")
#'
#' # plot main RCM effects
#' plotQUALYPSOeffect(QUALYPSO.synth,nameEff="RCM",xlab="Years")
#'
#' # plot fraction of total variance for the differences sources of uncertainty
#' plotQUALYPSOTotalVarianceDecomposition(QUALYPSO.synth,xlab="Years")
#'
#' # plot mean prediction and total variance with the differences sources of uncertainty
#' plotQUALYPSOMeanChangeAndUncertainties(QUALYPSO.synth,xlab="Years")
#'
#' #____________________________________________________________
#' # EXAMPLE OF QUALYPSO WHEN THE PREDICTOR IS TIME
#' #____________________________________________________________
#'
#' # list of options
#' listOption = list(typeChangeVariable='abs')
#'
#' # call QUALYPSO
#' QUALYPSO.time = QUALYPSO(Y=Y,scenAvail=scenAvail,X=X_time_vec,
#'                          Xfut=Xfut_time,listOption=listOption)
#'
#' # grand mean effect
#' plotQUALYPSOgrandmean(QUALYPSO.time,xlab="Years")
#'
#' # main GCM effects
#' plotQUALYPSOeffect(QUALYPSO.time,nameEff="GCM",xlab="Years")
#'
#' # main RCM effects
#' plotQUALYPSOeffect(QUALYPSO.time,nameEff="RCM",xlab="Years")
#'
#' # mean change and associated uncertainties
#' plotQUALYPSOMeanChangeAndUncertainties(QUALYPSO.time,xlab="Years")
#'
#' # variance decomposition
#' plotQUALYPSOTotalVarianceDecomposition(QUALYPSO.time,xlab="Years")
#'
#' #____________________________________________________________
#' # EXAMPLE OF QUALYPSO WHEN THE PREDICTOR IS THE GLOBAL TEMPERATURE
#' #____________________________________________________________
#'
#' # list of options
#' listOption = list(typeChangeVariable='abs')
#'
#' # call QUALYPSO
#' QUALYPSO.globaltas = QUALYPSO(Y=Y,scenAvail=scenAvail,X=X_globaltas,
#'                               Xfut=Xfut_globaltas,listOption=listOption)
#'
#' # grand mean effect
#' plotQUALYPSOgrandmean(QUALYPSO.globaltas,xlab="Global warming (Celsius)")
#'
#' # main GCM effects
#' plotQUALYPSOeffect(QUALYPSO.globaltas,nameEff="GCM",xlab="Global warming (Celsius)")
#'
#' # main RCM effects
#' plotQUALYPSOeffect(QUALYPSO.globaltas,nameEff="RCM",xlab="Global warming (Celsius)")
#'
#' # mean change and associated uncertainties
#' plotQUALYPSOMeanChangeAndUncertainties(QUALYPSO.globaltas,xlab="Global warming (Celsius)")
#'
#' # variance decomposition
#' plotQUALYPSOTotalVarianceDecomposition(QUALYPSO.globaltas,xlab="Global warming (Celsius)")
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @export
#'
#' @author Guillaume Evin
QUALYPSO = function(Y,scenAvail,X=NULL,Xfut=NULL,iFut=NULL,listOption=NULL){
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
  }else{
    # Y is a matrix: Scenario x Time
    nS = nrow(Y)
    nY = ncol(Y)

    # parallelization over time
    paralType = 'Time'
  }

  # X is the explanatory (or dependent) variable against the evolution of trajectory
  # are assessed. It is usually the years corresponding to the climate simulations
  # but can also be the global temperature.
  # X can be:
  # - NULL: in that case, we create a simple index corresponding to the size of Y
  # - a vector: the same depending var. is used for all the climate simulations,
  # its length must correspond to the number of columns of Y.
  # - a matrix: same size as Y, each line indicates the depending var. for this simulation.
  if(is.null(X)){
    Xmat = matrix(rep(1:nY,nS),byrow=T,nrow=nS,ncol=nY)
  }else if(is.vector(X)){
    if(nY!=length(X)){
      stop('if X is provided as a vector, its length must equal the number of columns of Y')
    }else{
      # repeat the vector to obtain a matrix
      Xmat = matrix(rep(X,nS),byrow=T,nrow=nS,ncol=nY)
    }
  }else if(is.matrix(X)){
    Xmat = X
    if(length(dim(Y))==2){
      if(any(dim(X)!=dim(Y))){
        stop('if X is provided as a matrix and Y as a grid, its size must match the size of Y')
      }
    }else if(length(dim(Y))==3){
      if(any(dim(X)!=dim(Y)[c(2,3)])){
        stop('if X is provided as a matrix, its first dimension must match the second
             dimension of Y (number of scenarios) and and its second dimension must match
             the third dimension of Y (number of future predictions)')
      }
    }
  }else{
    stop('X must be NULL, a vector or a matrix')
  }


  # Xfut are the values for X used to compute absolute or relative changes
  # usually indicating the future periods. When a grid is provided for Y, it must
  # be a single value. Indeed, in this case, we run QUALYPSO only for one future
  # year/global temperature. Otherwise, it can be a vector.
  # if Xfut is provided, we check that is a single value within the values of X
  if(is.null(Xfut)){
    if(is.null(X)){
      Xfut = 1:nY
    }else if(is.vector(X)){
      Xfut = X
    }else if(is.matrix(X)){
      Xfut = seq(from=min(X),to=max(X),length.out=10)
    }
  }else if(!(is.vector(Xfut)&is.numeric(Xfut))){
    stop('Xfut must be a vector of numerical values')
  }else if(length(Xfut)==1){
    stop('Xfut must be a vector with a length greater than 1')
  }else if(min(Xfut)<min(Xmat)|max(Xfut)>max(Xmat)){
    stop('Xfut must be within the range of X')
  }


  if(paralType == 'Grid'){
    # if Y is an array, iFut must be provided (we provide the ANOVA decomposition
    # only for one future time/global tas)
    if(is.null(iFut)){
      stop('if Y is an array, iFut must be provided')
    }else if(length(iFut)!=1|!is.numeric(iFut)){
      stop('iFut must be a single numeric value')
    }else if(!any(iFut==1:length(Xfut))){
      stop(paste0('iFut must be an index of the vector Xfut, i.e. an integer between 1 and ',length(Xfut)))
    }
  }

  # register clusters
  if(listOption$nCluster>1){
    cl <- parallel::makeCluster(listOption$nCluster)
    doParallel::registerDoParallel(cl)
  }


  ##############################################
  # check presence of NAs in climate projections
  if(paralType == 'Time'){
    # check is some simulation chains are entirely missing
    hasAllNa = apply(Y,1,function(x) all(is.na(x)))

    if(any(hasAllNa)){
      warning(paste0('Error in QUALYPSO: some projections have only NAs in Y: ',paste(which(hasAllNa),collapse = ',')))
      return(NULL)
    }
  }

  ##############################################
  # extract climate response
  if(paralType == 'Time'){
    if(is.null(listOption$climResponse)){
      climResponse = fit.climate.response(Y,
                                          spar=listOption$spar,
                                          Xmat=Xmat, Xfut=Xfut,
                                          typeChangeVariable=listOption$typeChangeVariable)
    }else{
      climResponse = listOption$climResponse
    }

    # extract quantities from these fits
    phiStar = climResponse$phiStar
    etaStar = climResponse$etaStar
    YStar = climResponse$YStar
    phi = climResponse$phi

    # internal variability
    varInterVariability = rep(climResponse$varInterVariability,length(Xfut))

    # phiStar for the ANOVA
    phiStar.ANOVA = phiStar

  }else if(paralType == 'Grid'){
    if(is.null(listOption$climResponse)){
      # initialise outputs
      phiStar = phi = array(dim=c(nG,nS,length(Xfut)))
      etaStar = YStar = array(dim=d)
      varInterVariability = vector(length=nG)

      # parallelize if required
      if(listOption$nCluster==1){
        for(g in 1:nG){
          # check is some simulation chains are entirely missing
          hasAllNa = apply(Y[g,,],1,function(x) all(is.na(x)))
          if(any(hasAllNa)){
            phiStar[g,,] = NA
            etaStar[g,,] = NA
            phi[g,,] = NA
            varInterVariability[g] = NA
          }else{
            climResponse = fit.climate.response(Y[g,,],
                                                spar=listOption$spar,
                                                Xmat=Xmat, Xfut=Xfut,
                                                typeChangeVariable=listOption$typeChangeVariable)

            # extract quantities from these fits
            phiStar[g,,] = climResponse$phiStar
            etaStar[g,,] = climResponse$etaStar
            YStar[g,,] = climResponse$YStar
            phi[g,,] = climResponse$phi
            varInterVariability[g] = climResponse$varInterVariability
          }
        }
      }else{
        g = NULL # avoid warning during check
        climResponse <- foreach(g=1:nG,.export=c("fit.climate.response")) %dopar% {
          # check is some simulation chains are entirely missing
          hasAllNa = apply(Y[g,,],1,function(x) all(is.na(x)))
          if(any(hasAllNa)){
            climResponse.g = list(phiStar = NA, etaStar = NA, phi = NA,
                                  varInterVariability = NA)
          }else{
            climResponse.g = fit.climate.response(Y[g,,],
                                                  spar=listOption$spar,
                                                  Xmat=Xmat, Xfut=Xfut,
                                                  typeChangeVariable=listOption$typeChangeVariable)
          }
          return(climResponse.g)
        }

        # fill the matrices
        for(g in 1:nG){
          phiStar[g,,] = climResponse[[g]]$phiStar
          etaStar[g,,] = climResponse[[g]]$etaStar
          phi[g,,] = climResponse[[g]]$phi
          varInterVariability[g] = climResponse[[g]]$varInterVariability
        }
      }
    }else{
      climResponse = listOption$climResponse
      YStar = climResponse$YStar
      phiStar = climResponse$phiStar
      etaStar = climResponse$etaStar
      phi = climResponse$phi
      varInterVariability = climResponse$varInterVariability
    }

    # phiStar for ANOVA
    phiStar.ANOVA = t(phiStar[,,iFut])
  }

  # names of the main effect
  if(is.null(colnames(scenAvail))){
    namesEff = paste0("Eff",1:ncol(scenAvail))
  }else{
    namesEff = colnames(scenAvail)
  }

  #====================================================================================================

  ################
  # Check for singularities

  # contrasts
  list.contrasts = list()
  for(j in 1:length(namesEff)){
    list.contrasts[[namesEff[j]]] = contr.sum
  }
  # formula
  formula = paste0("phiStar ~ ",paste0(namesEff,collapse = " + "))
  #lm
  lm.data = scenAvail
  lm.data$phiStar=phiStar.ANOVA[,ncol(phiStar.ANOVA)]
  lm.out = lm(formula, lm.data,contrasts=list.contrasts)
  if(any(is.na(lm.out$coefficients))){
    stop("singular fit encountered: the effects cannot be estimated (ill-posed problem)")
  }

  ##################
  # ANOVA on phiStar
  if(listOption$ANOVAmethod=="QUALYPSO"){
    anova = QUALYPSO.ANOVA(phiStar = phiStar.ANOVA, scenAvail = scenAvail, listOption = listOption, namesEff = namesEff)
  }else{
    anova = lm.ANOVA(phiStar = phiStar.ANOVA, scenAvail = scenAvail, listOption = listOption, namesEff = namesEff)
  }



  #====================================================================================================
  # further computation using point estimates

  # first retrieve point estimates and some quantities
  effhat = anova$MAINEFFECT
  muHat = anova$GRANDMEAN$MEAN
  listEff = anova$listScenarioInput$listEff
  nEff = length(effhat)
  nP = dim(phiStar.ANOVA)[2]

  # retrieve residual errors: differences between the climate change responses
  # and the additive effects (grand mean + main effect)
  mat.eff = matrix(nrow=nP,ncol=nEff)
  RESERR = matrix(nrow=nP,ncol=nS)
  for(iS in 1:nS){
    for(iE in 1:nEff){
      indE = which(scenAvail[iS,iE]==listEff[[iE]])
      mat.eff[,iE] = effhat[[namesEff[iE]]]$MEAN[,indE]
    }
    RESERR[,iS] = phiStar.ANOVA[iS,1:nP] - muHat - Rfast::rowsums(mat.eff)
  }

  # internal variability
  INTERNALVAR = varInterVariability


  # variance decomposition
  Vbind = cbind(anova$EFFECTVAR,anova$RESIDUALVAR$MEAN,INTERNALVAR)

  # Total variability
  TOTALVAR = Rfast::rowsums(Vbind)

  # Decomposition of the total uncertainty
  DECOMPVAR = Vbind/replicate(n = ncol(Vbind), TOTALVAR)
  colnames(DECOMPVAR) = c(namesEff,"ResidualVar","InternalVar")

  # stop clusters
  if(listOption$nCluster>1){
    parallel::stopCluster(cl)
  }


  #############################################
  # return results
  return(list(CLIMATEESPONSE=list(phiStar=phiStar,etaStar=etaStar,YStar=YStar,phi=phi),
              GRANDMEAN=anova$GRANDMEAN,
              MAINEFFECT=anova$MAINEFFECT,
              CHANGEBYEFFECT=anova$CHANGEBYEFFECT,
              EFFECTVAR=anova$EFFECTVAR,
              CONTRIB_EACH_EFFECT=anova$CONTRIB_EACH_EFFECT,
              RESIDUALVAR=anova$RESIDUALVAR,
              INTERNALVAR=INTERNALVAR,
              TOTALVAR=TOTALVAR,
              DECOMPVAR=DECOMPVAR,
              RESERR=RESERR,
              Xmat=Xmat,Xfut=Xfut,
              paralType=paralType,namesEff=namesEff,
              Y=Y,listOption=anova$listOption,
              listScenarioInput=anova$listScenarioInput))
}


#==============================================================================
#' plotQUALYPSOclimateResponse
#'
#' Plot the climate responses.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param lim y-axis limits (default is NULL)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOclimateResponse = function(QUALYPSOOUT,lim=NULL,xlab="X",ylab="Y",...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  phi = QUALYPSOOUT$CLIMATEESPONSE$phi

  # list of scenarios
  scenAvail = QUALYPSOOUT$listScenarioInput$scenAvail

  # Xmat and Y arguments
  Xmat = QUALYPSOOUT$Xmat
  Y = QUALYPSOOUT$Y

  # number of scenarios
  nS = nrow(Y)

  for(iS in 1:nS){
    Ys = Y[iS,]
    Xs = Xmat[iS,]
    phis = phi[iS,]

    plot(-1, -1, xlim = range(c(Xs,Xfut)), ylim = range(c(Ys,phis)),
         main=paste0(scenAvail[iS,],collapse = " / "),
         xlab = xlab, ylab = ylab, ...)

    # add lines of raw projection and climate projection
    lines(Xs, Ys, lwd = 1)
    lines(Xfut, phis, lwd = 3)

    # add legend
    legend("topleft",legend = c("Raw projection", "Climate response"),
           lty=1,lwd=c(1,3),bty="n")

    readline(prompt = "Press Enter")
  }
}


#==============================================================================
#' plotQUALYPSOclimateChangeResponse
#'
#' Plot climate change responses.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param lim y-axis limits (default is NULL)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOclimateChangeResponse = function(QUALYPSOOUT,lim=NULL,xlab="",
                                             ylab="Climate change response",...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  phiStar = QUALYPSOOUT$CLIMATEESPONSE$phiStar


  # initiate plot
  if(is.null(lim)) lim = range(phiStar)
  plot(-100,-100,xlim=range(Xfut),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  for(i in 1:nrow(phiStar)){
    lines(Xfut,phiStar[i,],lwd=3,col=i)
  }
}

#==============================================================================
#' plotQUALYPSOgrandmean
#'
#' Plot prediction of grand mean ensemble.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
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
plotQUALYPSOgrandmean = function(QUALYPSOOUT,lim=NULL,col='black',xlab="",
                                 ylab="Grand mean",addLegend=T,...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  meanPred = QUALYPSOOUT$GRANDMEAN$MEAN

  # retrieve limits
  binf = QUALYPSOOUT$GRANDMEAN$CI[,1]
  bsup = QUALYPSOOUT$GRANDMEAN$CI[,2]

  # colors polygon
  colPoly = adjustcolor(col,alpha.f=0.2)

  # initiate plot
  if(is.null(lim)) lim = range(c(binf,bsup),na.rm=TRUE)
  plot(-100,-100,xlim=range(Xfut),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  # add confidence interval
  polygon(c(Xfut,rev(Xfut)),c(binf,rev(bsup)),col=colPoly,lty=0)

  # add median
  lines(Xfut,meanPred,lwd=3,col=col)

  # legend
  if(addLegend){
    pctCI = round(QUALYPSOOUT$listOption$probCI*100)
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
#' @param nameEff name of the main effect to be plotted in \code{QUALYPSOOUT$namesEff}
#' @param includeMean if TRUE, the grand mean is added to the main effect in the plot
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
plotQUALYPSOeffect = function(QUALYPSOOUT,nameEff,includeMean=FALSE,lim=NULL,
                              col=1:20,xlab="",ylab="Effect",addLegend=TRUE,
                              ...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # index of this effect
  iEff = which(QUALYPSOOUT$namesEff==nameEff)
  if(length(iEff)==0) stop("wrong value for nameEff")

  # retrieve effects
  if(includeMean){
    EffHat = QUALYPSOOUT$CHANGEBYEFFECT[[nameEff]]
  }else{
    EffHat = QUALYPSOOUT$MAINEFFECT[[nameEff]]
  }
  nEff = dim(EffHat$MEAN)[2]

  # retrieve mean
  meanPred = EffHat$MEAN

  # add CI if the ANOVA method is QUALYPSO
  add.CI = QUALYPSOOUT$listOption$ANOVAmethod=="QUALYPSO"

  # initiate plot
  if(is.null(lim)){
    if(add.CI){
      lim = range(EffHat$CI,na.rm=TRUE)
    }else{
      lim = range(meanPred,na.rm=TRUE)
    }
  }
  plot(-100,-100,xlim=range(Xfut),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  for(i in 1:nEff){
    # colors polygon
    colPoly = adjustcolor(col[i],alpha.f=0.2)

    # add confidence interval
    if(add.CI){
      polygon(c(Xfut,rev(Xfut)),c(EffHat$CI[,1,i],rev(EffHat$CI[,2,i])),col=colPoly,lty=0)
    }

    # add median
    lines(Xfut,meanPred[,i],lwd=3,col=col[i])
  }

  # legend
  if(addLegend){
    pctCI = round(QUALYPSOOUT$listOption$probCI*100)
    if(QUALYPSOOUT$listOption$ANOVAmethod=="QUALYPSO"){
      legend('topleft',bty='n',fill=c(NA,'black'),lwd=c(2,NA),lty=c(1,NA),
             border=c(NA,'black'),col=c('black',NA),
             legend=c('Median',paste0(pctCI,'%CI')))
    }

    legend('bottomleft',bty='n',lwd=2,lty=1,col=col,
           legend=QUALYPSOOUT$listScenarioInput$listEff[[iEff]])
  }
}


#==============================================================================
# printMessageDimension
printMessageDimension = function(Y, scenAvail, X){
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

  # X
  print("X must be a vector/matrix of years or global temperatures corresponding to the dimension of Y:")
  print(X)
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
                                                  xlab="",ylab="% Total Variance",addLegend=TRUE,...){
  # future predictor values
  Xfut = QUALYPSOOUT$Xfut
  nFut = length(Xfut)
  nEff = QUALYPSOOUT$listScenarioInput$nEff

  # number of main effects
  if(is.null(vecEff)){
    vecEff = 1:nEff
  }

  # Variance decomposition
  VARDECOMP = QUALYPSOOUT$DECOMPVAR
  VARDECOMP[,1:nEff] = VARDECOMP[,vecEff]

  # figure
  col = col[1:(nEff+2)]
  cum=rep(0,nFut)
  plot(-1,-1,xlim=range(Xfut),ylim=c(0,1),xaxs="i",yaxs="i",las=1,xlab=xlab,ylab=ylab,...)
  for(i in 1:(nEff+2)){
    cumPrevious = cum
    cum = cum + VARDECOMP[,i]
    polygon(c(Xfut,rev(Xfut)),c(cumPrevious,rev(cum)),col=rev(col)[i],lty=1)
  }
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    legend('topleft',bty='n',cex=1.1, fill=rev(col),
           legend=c(QUALYPSOOUT$namesEff[vecEff],'Res. Var.','Int. Variab.'))
  }
}


#==============================================================================
#' plotQUALYPSOTotalVarianceByScenario
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param nameEff name of the main effect to be plotted in \code{QUALYPSOOUT$namesEff}
#' @param nameScenario name of the scenario to be plotted (as provided in \code{scenAvail})
#' @param col colors for each source of uncertainty, the first two colors
#' corresponding to internal variability and residual variability, respectively
#' @param ylim y-axis limits
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOTotalVarianceByScenario = function(QUALYPSOOUT,nameEff,nameScenario,
                                               col=NULL,ylim=NULL,xlab="",
                                               ylab="Change variable",
                                               addLegend=TRUE,...){
  # future predictor values
  Xfut = QUALYPSOOUT$Xfut
  nFut = length(Xfut)

  # which scenario
  iEff = which(QUALYPSOOUT$namesEff==nameEff)
  if(length(iEff)==0) stop("wrong value for nameEff")
  iScenario = which(QUALYPSOOUT$listScenarioInput$listEff[[iEff]] == nameScenario)

  # mean prediction
  meanPred = QUALYPSOOUT$CHANGEBYEFFECT[[nameEff]]$MEAN[,iScenario]

  # remove effect corresponding to the scenarios from the total variance
  Veff = QUALYPSOOUT$EFFECTVAR[,-iEff]

  # concatenate variances
  Vbind = rbind(t(Veff),QUALYPSOOUT$RESIDUALVAR$MEAN,QUALYPSOOUT$INTERNALVAR)
  nEff = nrow(Vbind)-2
  Vtot = colSums(Vbind)
  Vnorm = Vbind/t(replicate(n = nrow(Vbind), Vtot))

  # reverse
  vNormRev = apply(Vnorm,2,rev)

  # compute the lower bound if the distribution is gaussian
  probCI = QUALYPSOOUT$listOption$probCI
  binf = qnorm(p = (1-probCI)/2, mean = meanPred, sd = sqrt(Vtot))
  bsup = qnorm(p = 0.5+probCI/2, mean = meanPred, sd = sqrt(Vtot))

  # figure
  if(is.null(col)){
    default.col = c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1")
    col = default.col[1:(nEff+2)]
  }

  # obtain limits of the intervals, proportion corresponds to the part of the variance, lower and upper than the mean
  limIntInf = limIntSup = matrix(nrow=nEff+2,ncol=nFut)
  limIntInf[1,] = binf
  limIntSup[1,] = bsup
  for(i in 1:(nEff+1)){
    limIntInf[i+1,] = limIntInf[i,]+vNormRev[i,]*(meanPred-binf)
    limIntSup[i+1,] = limIntSup[i,]-vNormRev[i,]*(bsup-meanPred)
  }

  # figure
  if(is.null(ylim)) ylim = c(min(binf),max(bsup))
  plot(-1,-1,xlim=range(Xfut),ylim=ylim,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1,...)
  for(i in 1:(nEff+2)){
    polygon(c(Xfut,rev(Xfut)),c(limIntInf[i,],rev(limIntSup[i,])),col=col[i],lty=1)
  }
  lines(Xfut,meanPred,col="white",lwd=1)

  # add horizontal lines
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    namesEff = QUALYPSOOUT$namesEff[-iEff]

    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}


#==============================================================================
#' plotQUALYPSOMeanChangeAndUncertainties
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
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
plotQUALYPSOMeanChangeAndUncertainties = function(QUALYPSOOUT,col=NULL,ylim=NULL,
                                                  xlab="",ylab="Change variable",addLegend=TRUE,...){

  # probCI
  probCI = QUALYPSOOUT$listOption$probCI

  # future predictor values
  Xfut = QUALYPSOOUT$Xfut
  nFut = length(Xfut)

  # mean prediction
  meanPred = QUALYPSOOUT$GRANDMEAN$MEAN
  vNorm = t(QUALYPSOOUT$DECOMPVAR)
  Vtot = QUALYPSOOUT$TOTALVAR
  nEff = nrow(vNorm)-2

  # reverse
  vNormRev = apply(vNorm,2,rev)

  # compute the lower bound if the distribution is gaussian
  binf = qnorm(p = (1-probCI)/2, mean = meanPred, sd = sqrt(Vtot))
  bsup = qnorm(p = 0.5+probCI/2, mean = meanPred, sd = sqrt(Vtot))

  # figure
  if(is.null(col)){
    default.col = c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1")
    col = default.col[1:(nEff+2)]
  }

  # obtain limits of the intervals, proportion corresponds to the part of the variance, lower and upper than the mean
  limIntInf = limIntSup = matrix(nrow=nEff+2,ncol=nFut)
  limIntInf[1,] = binf
  limIntSup[1,] = bsup
  for(i in 1:(nEff+1)){
    limIntInf[i+1,] = limIntInf[i,]+vNormRev[i,]*(meanPred-binf)
    limIntSup[i+1,] = limIntSup[i,]-vNormRev[i,]*(bsup-meanPred)
  }

  # figure
  if(is.null(ylim)) ylim = c(min(binf),max(bsup))
  plot(-1,-1,xlim=range(Xfut),ylim=ylim,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1,...)
  for(i in 1:(nEff+2)){
    polygon(c(Xfut,rev(Xfut)),c(limIntInf[i,],rev(limIntSup[i,])),col=col[i],lty=1)
  }
  lines(Xfut,meanPred,col="white",lwd=1)

  # add horizontal lines
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    namesEff = QUALYPSOOUT$names
    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}


#==============================================================================
#' plotQUALYPSOMeanChangeAndUncertaintiesBetatest
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
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
plotQUALYPSOMeanChangeAndUncertaintiesBetatest = function(QUALYPSOOUT,col=NULL,ylim=NULL,
                                                          xlab="",ylab="Change variable",addLegend=TRUE,...){

  # probCI
  probCI = QUALYPSOOUT$listOption$probCI

  # future predictor values
  Xfut = QUALYPSOOUT$Xfut
  nFut = length(Xfut)

  # phiStar
  phiStar = QUALYPSOOUT$CLIMATEESPONSE$phiStar

  # Ystar: change variable from the projection (Y(t)-phi(c))
  YStar = QUALYPSOOUT$CLIMATEESPONSE$YStar

  # compute a multiplicative factor SDtot/SD(phi*) to match the standard deviation
  # obtained from QUALYPSO
  sd.qua = sqrt(QUALYPSOOUT$TOTALVAR-QUALYPSOOUT$INTERNALVAR)
  sd.emp = apply(phiStar,2,sd)
  sd.emp[sd.emp==0] = 1 # avoid NaN for the reference year
  sd.corr = sd.qua/sd.emp
  phiStar.corr = phiStar*t(replicate(nrow(phiStar),sd.corr))

  # mean prediction
  meanPred = QUALYPSOOUT$GRANDMEAN$MEAN

  # variance decomposition without the internal var.
  Vbind = cbind(QUALYPSOOUT$EFFECTVAR,QUALYPSOOUT$RESIDUALVAR$MEAN)
  Vtot = Rfast::rowsums(Vbind)
  DECOMPVAR = Vbind/replicate(n = ncol(Vbind), Vtot)
  vNorm = t(DECOMPVAR)
  nEff = nrow(vNorm)-1
  vNormRev = apply(vNorm,2,rev)

  # compute the lower bound if the distribution is gaussian
  binf = apply(phiStar.corr,2,quantile,probs = (1-probCI)/2)
  bsup = apply(phiStar.corr,2,quantile,probs = 0.5+probCI/2)

  # obtain limits of the intervals, proportion corresponds to the part of the variance, lower and upper than the mean
  limIntInf = limIntSup = matrix(nrow=nEff+1,ncol=nFut)
  limIntInf[1,] = binf
  limIntSup[1,] = bsup
  for(i in 1:nEff){
    limIntInf[i+1,] = limIntInf[i,]+vNormRev[i,]*(meanPred-binf)
    limIntSup[i+1,] = limIntSup[i,]-vNormRev[i,]*(bsup-meanPred)
  }

  # figure
  if(is.null(col)){
    default.col = c("yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1")
    col = default.col[1:(nEff+1)]
  }

  # start figure
  if(is.null(ylim)) ylim = range(YStar,na.rm=T)
  plot(-1,-1,xlim=range(Xfut),ylim=ylim,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1,...)

  # add raw climate change projections
  Xmat = QUALYPSOOUT$Xmat
  for(iS in 1:nrow(YStar)){
    lines(Xmat[iS,],YStar[iS,],lwd=0.5,col="gray")
  }

  # add intervals
  for(i in 1:(nEff+1)){
    polygon(c(Xfut,rev(Xfut)),c(limIntInf[i,],rev(limIntSup[i,])),col=col[i],lty=1)
  }

  # mean climate change response
  lines(Xfut,meanPred,col="white",lwd=1)

  # add horizontal lines
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)

  # legend
  if(addLegend){
    namesEff = QUALYPSOOUT$names
    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.'))
  }
}
