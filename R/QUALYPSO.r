###===============================###===============================###
### Guillaume Evin
### 03/11/2026, Grenoble
###  INRAE
### guillaume.evin@inrae.fr
###
### These functions provide an assessment of single-member and incomplete
###  ensembles of climate projections.
### It provides unbiased estimates of climate change responses of all
### simulation chains and of all uncertainty variables. It additionally propagates
### uncertainty due to missing information in the estimates.
###
### REFERENCES
###
### First development of the Bayesian approach:
### - Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
### Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
### Journal of Climate. J. Climate, 32, 2423–2440. <doi:10.1175/JCLI-D-18-0606.1>.
###
### The notations of the package corresponding to the following article:
### - Evin, G., Hingray, et al. "Uncertainty sources in a large ensemble of hydrological projections: Regional
###  Climate Models and Internal Variability matter, Hydrol. Earth Syst. Sci., 30, 1023–1051, 
### https://doi.org/10.5194/hess-30-1023-2026, 2026. 
###===============================###===============================###


#==============================================================================
#' Provide Qstar matrix containing Helmert contrasts
#'
#' Provide matrix containing Helmert contrasts (see Eq. A7 in Evin et al., 2019).
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
get_Qstar_mat = function(p){
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
#' Provide the matrix Q derived from Qstar
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
get_Qmat = function(p){
  # get Qstar
  Qstar = get_Qstar_mat(p)

  # get Q: Eq. A6 in Evin et al. (2019)
  Q = Qstar %*% MASS::ginv(expm::sqrtm(t(Qstar) %*% Qstar))

  return(Q)
}


#==============================================================================
#' Fit the climate response for each simulation chain 
#'
#' Fit the climate response for each simulation chain of an ensemble of \code{nS} projections. 
#' Each simulation chain is a time series of \code{nY} time steps (e.g. number of years).
#'
#' @param Y matrix of simulation chains: \code{nS} x \code{nY}
#' @param typeClimateResponse type of climate response fitted to the projections.
#' @param parClimateResponse parameter of the model applied for the extraction of the climate response.
#' @param Xmat matrix of predictors corresponding to the projections, e.g. time or global temperature.
#' @param Xfut values of the predictor over which the ANOVA will be applied.
#' @param typeChangeVariable type of change variable: "abs" (absolute, value by default) or "rel" (relative)
#'
#' @return list with the following fields for each simulation chain:
#' \itemize{
#'   \item \strong{CR_tfut}: \code{nS x nF}, climate responses
#'   \item \strong{CCR_tfut}: \code{nS x nF}, climate change responses
#'   \item \strong{D_t}: \code{nS x nY}, deviation from the climate change response
#'   \item \strong{CHANGEVAR_t}: \code{nS x nY}, change variable
#'   \item \strong{IV_tfut}: \code{nF}, internal variability
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
fit_ClimateResponse = function(Y, typeClimateResponse, parClimateResponse, Xmat, Xfut, typeChangeVariable){

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
  CCR_tfut = CR_tfut = matrix(nrow=nS,ncol=nFut)
  D_t = CHANGEVAR_t = matrix(nrow=nS,ncol=nY)

  for(iS in 1:nS){
    # projection for this simulation chain
    Ys = Y[iS,]
    Xs = Xmat[iS,]

    # fit a smooth signal (smooth cubic splines)
    zz = !is.na(Ys)

    # data.frame with Xs and Ys
    dfpred = data.frame(x=Xs[zz],y=Ys[zz])

    if(typeClimateResponse == "smooth.spline"){
      # call smooth.spline
      smooth.spline.out<- stats::smooth.spline(x=dfpred$x,y=dfpred$y,df=parClimateResponse)

      # fitted responses at the points of the fit (for D_t), at Xfut and at Xref
      CR_t = predict(smooth.spline.out, Xs)$y
      CR_tfut_S = predict(smooth.spline.out, Xfut)$y
      CR_tref = predict(smooth.spline.out, Xref)$y
    }else if(typeClimateResponse == "tweedie"){
      if(any(dfpred$y < 0)){
        stop("the tweedie distributions are designed for non-negative, right-skewed 
        data with a mass at zero. They do not natively support negative input values 
        for the response variable, as they are typically used for continuous positive data.
        Please use other options for typeClimateResponse.")
      }
      # Generalized Linear Model (GLM) with Tweedie distribution
      glm_out <- stats::glm(y ~ x, 
        family = statmod::tweedie(var.power = parClimateResponse, link.power = 0), 
        data = dfpred)
      
      # fitted responses at the points of the fit (for D_t), at Xfut and at Xref
      CR_t = predict(glm_out, newdata = data.frame(x = Xs), type = "response")
      CR_tfut_S = predict(glm_out, newdata = data.frame(x = Xfut), type = "response")
      CR_tref = predict(glm_out, newdata = data.frame(x = Xref), type = "response")      
    }else if(typeClimateResponse == "poly"){
      # polynomial regression (OLS)
      lm_out = stats::lm(y ~ stats::poly(x, parClimateResponse), data = dfpred)
      
      # fitted responses at the points of the fit (for D_t), at Xfut and at Xref
      CR_t = predict(lm_out, newdata = data.frame(x = Xs))
      CR_tfut_S = predict(lm_out, newdata = data.frame(x = Xfut))
      CR_tref = predict(lm_out, newdata = data.frame(x = Xref))      
    }else if(typeClimateResponse == "lowess"){
      # fit a lowess model with a 0.3 degree of smoothing, see:
      # Scherrer, Simon C., Cees de Valk, Michael Begert, Stefanie Gubler, Sven Kotlarski, et
      # Mischa Croci-Maspoli. « Estimating trends and the current climate mean in a changing
      # climate ». Climate Services 33 (2024): 100428. https://doi.org/10.1016/j.cliser.2023.100428.
      lowess_out <- lowess(x=dfpred$x,y=dfpred$y, f=parClimateResponse, iter=0, delta=0)

      # Create an interpolation function
      unique_df <- unique(data.frame(x=lowess_out$x, y=lowess_out$y))
      interpolation_fn <- stats::approxfun(unique_df$x, unique_df$y)

      # fitted responses at the points of the fit (for D_t), at Xfut and at Xref
      CR_t = interpolation_fn(Xs)
      CR_tfut_S = interpolation_fn(Xfut)
      CR_tref = interpolation_fn(Xref)      
    }

    # store climate response for this simulation chain
    CR_tfut[iS,] = CR_tfut_S

    # Climate change response: CCR_tfut, and internal variability expressed as a change: D_t
    if(typeChangeVariable=='abs'){
      # Eq. 5
      CCR_tfut[iS,] = CR_tfut_S-CR_tref
      D_t[iS,] = Ys-CR_t
      CHANGEVAR_t[iS,] = Ys-CR_tref
    }else if(typeChangeVariable=='rel'){
      # Eq. 6
      CCR_tfut[iS,] = CR_tfut_S/CR_tref-1
      D_t[iS,] = (Ys-CR_t)/CR_tref
      CHANGEVAR_t[iS,] = (Ys-CR_tref)/CR_tref
    }else{
      stop("fit_ClimateResponse: argument type.change.var must be equal to 'abs' (absolute changes) or 'rel' (relative changes)")
    }
  }

  # Heteroscedasticity: assume gaussian deviations from the climate responses
  # also assumes that climate responses have removed the main trends (mu=0)
  df_fit = data.frame(x = as.vector(Xmat), y = as.vector(D_t))
  gamlss_out = gamlss(y ~ x, sigma.formula = ~x, family = gamlss.dist::NO(), data = df_fit,mu.fix = TRUE, mu.start = 0,
  control = gamlss.control(trace = FALSE))

  # predict this model for Xfut
  sigma_pred = predict(object = gamlss_out, what = "sigma", type = "response", newdata = data.frame(x=Xfut), data = df_fit)
  IV_tfut = sigma_pred^2


  # return objects
  return(list(CCR_tfut=CCR_tfut,D_t=D_t,CHANGEVAR_t=CHANGEVAR_t,CR_tfut=CR_tfut,IV_tfut=IV_tfut))
}


#==============================================================================
#' Perform ANOVA using a Bayesian approach for a future time i
#'
#' Partition sources of uncertainty in climate change responses for one lead time or one grid point.
#'
#' @param CCR_tfut_i vector of \code{nS} climate change response for one lead time: \code{nS x 1}
#' @param nMCMC number of MCMC simulation required
#' @param listScenario list containing specifications, provided by \code{\link{process_scen}}
#'
#' @return  list with the following fields:
#' \itemize{
#'   \item \strong{mu}: vector of length \code{nMCMC}, mean climate change response
#'   \item \strong{sigma2}: vector of length \code{nMCMC}, variance of the residual terms
#'   \item \strong{effect}: list with \code{nEffect} elements, where each element corresponds to a different factor
#'   Each element is a matrix \code{nMCMC} x \code{nMainFactor}, and \code{nMainFactor} is the number of main effects (e.g. number of GCMs, RCMs, etc.)
#' }
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @author Guillaume Evin
Bayesian_ANOVA_i = function(CCR_tfut_i, nMCMC, listScenario){
  #============= retrieve objects related to the scenarios =====
  listEff=listScenario$listEff
  scen=listScenario$scen
  scenComp=listScenario$scenComp
  nFactor=listScenario$nFactor
  nEffect=listScenario$nEffect
  nComp=listScenario$nComp
  isMissing=listScenario$isMissing
  nMissing=listScenario$nMissing
  iMatchScen=listScenario$iMatchScen
  indexEffInCompScen=listScenario$indexEffInCompScen
  Qmat=listScenario$Qmat


  #=============  initialize effects =================
  effect.POST = list()

  if(any(is.na(CCR_tfut_i))){
    # when i corresponds to the reference period, CCR_tfut_i = 0 and no ANOVA can be performed
    # return 0 instead of NA for representations
    mu.POST = sigma2.POST = rep(0,nMCMC)
    for(i.eff in 1:nFactor) effect.POST[[i.eff]] = matrix(0,nrow=nMCMC,ncol=nEffect[i.eff])
    return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST))
  }

  if(all(CCR_tfut_i==0)){
    # when i corresponds to the reference period, CCR_tfut_i = 0 and no ANOVA can be performed
    # return 0 instead of NA for representations
    mu.POST = sigma2.POST = rep(0,nMCMC)
    for(i.eff in 1:nFactor) effect.POST[[i.eff]] = matrix(0,nrow=nMCMC,ncol=nEffect[i.eff])
    return(list(mu=mu.POST,sigma2=sigma2.POST,effect=effect.POST))
  }else{
    # phi to complete
    phi2comp = CCR_tfut_i[iMatchScen]
    # initialise matrix and arrays
    mu.POST = sigma2.POST = vector(length=nMCMC)
    for(i.eff in 1:nFactor) effect.POST[[i.eff]] = matrix(nrow=nMCMC,ncol=nEffect[i.eff])
  }



  ###################################################################
  # Bayesian sampling is done with Gibbs algorithm. Details about the choice of the prior distributions,
  # the hyperparameters, and full conditional posterior distributions are given in Appendix in Evin et al. (2019)

  #============= hyper-parameters: see Appendix A, section g. =============
  mu.mu = mean(CCR_tfut_i,na.rm=T)
  largevar = 16*var(as.vector(CCR_tfut_i),na.rm=T)
  sig2.mu = largevar
  lam.sig = 0.5
  X.diff = CCR_tfut_i-mean(CCR_tfut_i,na.rm=T)
  for(i.eff in 1:nFactor){
    eff.hat = aggregate(x = X.diff, by = list(scen[,i.eff]), FUN = "mean")
    X.diff = X.diff - eff.hat$x[match(scen[,i.eff],eff.hat$Group.1)]
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
  mat.eff = matrix(nrow=nComp,ncol=nFactor)
  for(i.eff in 1:nFactor){
    n = nEffect[i.eff]
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

    for(i.eff in 1:nFactor){
      # number of effects for each factor (i.e. number of GCMs, RCMs)
      n = nEffect[i.eff]
      # marginal sums
      mar.diff = Rfast::group(X.shift,indexEffInCompScen[,i.eff],method="sum")
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
#' Process scen
#'
#' Process factors describing each scenario (i.e each simulation chain)
#'
#' @param scen data.frame \code{nS} x \code{nFactor} with the \code{nFactor} factors (e.g. type of GCM) for each of the \code{nS} x \code{nS} scenarios
#'
#' @return list of preprocessed objects (\code{listEff, scen, scenComp, nFactor, nEffect, nComp, isMissing, nMissing, iMatchScen,
#' indexEffInCompScen, Qmat})
#'
#' @author Guillaume Evin
process_scen = function(scen){
  # number of scenarios
  nS = nrow(scen)

  # list of effects
  nFactor = ncol(scen)
  listEff = list()
  for(i in 1:nFactor) listEff[[i]] = unique(scen[,i])
  nEffect = unlist(lapply(listEff,length))

  # possible combinations of main effects
  scenComp = expand.grid(listEff)
  nComp = nrow(scenComp)


  #########  Missing scenarios   #########
  vScenComp <- apply(scenComp, 1, paste, collapse='.')
  vscen <- apply(scen, 1, paste, collapse='.')
  isMissing = !vScenComp%in%vscen
  nMissing = sum(isMissing)


  #########   vector of projections to complete with data augmentation   #########
  iMatchScen = match(vScenComp,vscen)


  #########   matrix of effect index: for each main effect, index of 'scenComp' (combinations of scenarios) related to listEff   #########
  indexEffInCompScen = matrix(nrow=nComp,ncol=nFactor)
  for(i.eff in 1:nFactor){
    indexEffInCompScen[,i.eff] = match(scenComp[,i.eff],listEff[[i.eff]])
  }


  #########   transformation matrices Q  #########
  Qmat = list()
  for(i.eff in 1:nFactor){
    Qmat[[i.eff]] = get_Qmat(nEffect[i.eff])
  }


  return(list(listEff=listEff,
              scen=scen, scenComp=scenComp,
              nFactor=nFactor, nEffect=nEffect, nComp=nComp,
              isMissing=isMissing, nMissing=nMissing,
              iMatchScen=iMatchScen,
              indexEffInCompScen=indexEffInCompScen,
              Qmat=Qmat))
}


#==============================================================================
#' Check validity of options provided in listOption
#'
#' Check if input options provided in \code{\link{QUALYPSO}} are valid and assign default values if missing.
#'
#' @param listOption list of options
#'
#' @return List containing the complete set of options.
#'
#' @author Guillaume Evin
check_option = function(listOption){
  if(is.null(listOption)){
    listOption = list()
  }

  # typeClimateResponse
  if('typeClimateResponse' %in% names(listOption)){
    typeClimateResponse = listOption[['typeClimateResponse']]
    if(!(typeClimateResponse%in%c('smooth.spline','poly','lowess','tweedie'))){
      stop("typeClimateResponse must be equal to: 'smooth.spline', 'poly', 'lowess', or 'tweedie'")
    }
  }else{
    listOption[['typeClimateResponse']] = "smooth.spline"
  }

  # parClimateResponse
  if('parClimateResponse' %in% names(listOption)){
    parClimateResponse = listOption[['parClimateResponse']]
    if(!(is.numeric(parClimateResponse) | length(parClimateResponse)>1)){
      stop("parClimateResponse must be one numerical value")
    }
  }else{
    if(listOption[['typeClimateResponse']] == "smooth.spline"){
      listOption[['parClimateResponse']] = 4
    }else if(listOption[['typeClimateResponse']] == "poly"){
      listOption[['parClimateResponse']] = 1
    }else if(listOption[['typeClimateResponse']] == "lowess"){
      listOption[['parClimateResponse']] = 0.4
    }else if(listOption[['typeClimateResponse']] == "tweedie"){
      listOption[['parClimateResponse']] = 1.5
    }
  }

  # ANOVAmethod
  if('ANOVAmethod' %in% names(listOption)){
    ANOVAmethod = listOption[['ANOVAmethod']]
    if(!(ANOVAmethod%in%c('Bayesian','lm'))) stop("ANOVAmethod must be equal to 'Bayesian' or 'lm'")
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
  listOption[['version']] = 'v3.2'

  return(listOption)
}


#==============================================================================
#' Perform ANOVA using a Bayesian approach
#' 
#' Partition uncertainty in climate responses using an ANOVA inferred with a Bayesian approach.
#'
#' @param CCR_tfut matrix of climate change responses: \code{nS} x \code{nFut}.
#' @param scen data.frame \code{nS} x \code{nFactor} for the \code{nS} scenarios and the \code{nFactor} factors (e.g. type of GCM)
#' @param listOption list of options (see \code{\link{QUALYPSO}})
#' @param namesEff names of the main effects
#'
#' @return list with the following fields:
#'
#' \itemize{
#'   \item \strong{M_tfut}: List of estimates for the grand mean:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut} of posterior means
#'      \item \strong{SD}: vector of length \code{nFut} of posterior standard dev.
#'      \item \strong{CI}: matrix \code{nFut} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'      \item \strong{QUANT}: matrix \code{nFut} x \code{nQ} of quantiles related
#'      to the probabilities \code{quantilePosterior} given in \code{listOption}
#'   }
#'   \item \strong{RV_tfut}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut} of posterior means
#'      \item \strong{SD}: vector of length \code{nFut} of posterior standard dev.
#'      \item \strong{CI}: matrix \code{nFut} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'      \item \strong{QUANT}: matrix \code{nFut} x \code{nQ} of quantiles related
#'      to the probabilities \code{quantilePosterior} given in \code{listOption}
#'   }
#'   \item \strong{E_tfut}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item \strong{MEAN}: matrix \code{nFut} x \code{nEffect} of posterior means
#'      \item \strong{SD}: matrix \code{nFut} x \code{nEffect} of posterior standard dev.
#'      \item \strong{CI}: array \code{nFut} x 2 x \code{nEffect} of credible
#'      intervals of probability \code{probCI} given in \code{listOption}.
#'      \item \strong{QUANT}: array \code{nFut} x \code{nQ} x \code{nEffect} of
#'      quantiles related to the probabilities \code{quantilePosterior} given in
#'      \code{listOption}
#'   }
#'   \item \strong{VARE_tfut}: variability related to the main effects (i.e.
#'   variability between the different RCMs, GCMs,..). Matrix \code{nFut} x
#'   \code{nEffect}
#'   \item \strong{CONTRIB_VARE_tfut}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{nFut} x \code{nEffect}
#'   \item \strong{listOption}: list of options used to obtained these results
#'    (obtained from \code{\link{check_option}})
#'   \item \strong{listScenario}: list of scenario factors
#'    (obtained from \code{\link{process_scen}})
#' }
#'
#' @export
#'
#' @author Guillaume Evin
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
Bayesian_ANOVA = function(CCR_tfut,scen,listOption=NULL,namesEff){
  #########  process input #########
  # number of future times
  n = dim(CCR_tfut)[2]

  # number of MCMC samples
  nKeep = listOption$nKeep
  vec.keep = (listOption$nBurn+1):listOption$nMCMC

  # Process scenarios data.frame to get different objects
  listScenario = process_scen(scen = scen)
  nFactor = listScenario$nFactor
  nEffect = listScenario$nEffect

  #########  Apply ANOVA for each time step: parallel computation over time steps #########
  # apply parallel computation only if more than one cluster has been indicated
  Anova.POST = list()
  for(i in 1:n){
    Anova.POST[[i]] = Bayesian_ANOVA_i(CCR_tfut_i=CCR_tfut[,i], nMCMC=listOption$nMCMC,listScenario = listScenario)
  }


  #-------------------------------------------------------------------------
  # POSTPROCESS: EXTRACT SOME QUANTITIES
  #-------------------------------------------------------------------------
  qCI = c((1-listOption$probCI)/2,0.5+listOption$probCI/2)
  qPost = listOption$quantilePosterior


  #============================
  # M_tfut
  #============================
  M_tfut = list()

  # extract posterior
  mu.POST = matrix(nrow=n,ncol=nKeep)
  for(i in 1:n) mu.POST[i,]=Anova.POST[[i]]$mu[vec.keep]

  # mean
  M_tfut$MEAN = apply(mu.POST,1,mean)

  # sd
  M_tfut$SD = apply(mu.POST,1,sd)

  # CI
  M_tfut$CI = t(apply(mu.POST,1,quantile,probs=qCI))

  # QUANT
  M_tfut$QUANT = t(apply(mu.POST,1,quantile,probs=qPost))


  #============================
  # RV_tfut
  #============================
  # variance of residual effects, With the Bayesian estimation, this variance
  # is part of the inference (sigma^2) and we take the mean of the posterior as the point estimate
  RV_tfut = list()

  # extract posterior
  sigma2.POST = matrix(nrow=n,ncol=nKeep)
  for(i in 1:n) sigma2.POST[i,]=Anova.POST[[i]]$sigma2[vec.keep]

  # mean
  RV_tfut$MEAN = apply(sigma2.POST,1,mean)

  # sd
  RV_tfut$SD = apply(sigma2.POST,1,sd)

  # CI
  RV_tfut$CI = t(apply(sigma2.POST,1,quantile,probs=qCI))

  # QUANT
  RV_tfut$QUANT = t(apply(sigma2.POST,1,quantile,probs=qPost))


  #==============================================================
  # E_tfut + VARE_tfut + CONTRIB_VARE_tfut
  #==============================================================
  # main effects
  E_tfut = list()

  # variance of the main effects
  VARE_tfut = matrix(nrow=n,ncol=nFactor)

  # contribution to the variance of individual effects
  CONTRIB_VARE_tfut = list()

  # main effects
  for(i.eff in 1:nFactor){
    eff = namesEff[i.eff]
    #============
    # E_tfut
    #============
    # trim posterior distributions of the main effects
    eff.POST = array(dim=c(n,nKeep,nEffect[i.eff]))
    for(i in 1:n) eff.POST[i,,]=Anova.POST[[i]]$effect[[i.eff]][vec.keep,]

    # retrieve estimates for the main effects
    lHatMain =  list()
    lHatMain$MEAN = apply(eff.POST,c(1,3),mean)
    lHatMain$SD = apply(eff.POST,c(1,3),sd)
    eff.ci = apply(eff.POST,c(1,3),quantile,probs=qCI)
    lHatMain$CI = aperm(eff.ci, c(2,3,1))
    eff.quant = apply(eff.POST,c(1,3),quantile,probs=qPost)
    lHatMain$QUANT = aperm(eff.quant, c(2,3,1))
    E_tfut[[eff]] = lHatMain

    #============
    # VARE_tfut
    #============
    # predictive variance: mean of variances, which correspond to the mean of squares since the mean is 0 by constraint (see Eq 16, 17 and 18)
    VARE_tfut[,i.eff] = apply(eff.POST^2,1,mean)

    #============
    # CONTRIB_VARE_tfut
    #============
    matContrib = matrix(nrow=n,ncol=nEffect[i.eff])
    for(iIndEff in 1:nEffect[i.eff]){
      # individual effect squared
      indEff2 = apply(eff.POST[,,iIndEff]^2,1,mean)
      # variance of this main effect x number of factors
      infEff2Tot = (VARE_tfut[,i.eff]*nEffect[i.eff])
      # contribution of this individual effect, as a percentage
      matContrib[,iIndEff] = indEff2/infEff2Tot
    }
    CONTRIB_VARE_tfut[[eff]] = matContrib
  }


  # return results
  return(list(M_tfut=M_tfut,
              RV_tfut=RV_tfut,
              E_tfut=E_tfut,
              VARE_tfut=VARE_tfut,
              CONTRIB_VARE_tfut=CONTRIB_VARE_tfut,
              listOption=listOption,
              listScenario=listScenario))
}


#==============================================================================
#' Perform ANOVA using a linear model
#'
#' Partition uncertainty in climate responses using an ANOVA inferred with a linear model.
#'
#' @param CCR_tfut matrix of climate change responses: \code{nS} x \code{nFut}.
#' @param scen data.frame \code{nS} x \code{nFactor} for each of the \code{nS} scenarios and the \code{nFactor} factors (e.g. type of GCM) 
#' @param listOption list of options (see \code{\link{QUALYPSO}})
#' @param namesEff names of the main effects
#'
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{M_tfut}: List of estimates for the grand mean:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut} of means
#'      \item \strong{SD}: vector of length \code{nFut} of standard dev.
#'      \item \strong{CI}: matrix \code{nFut} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}.
#'   }
#'   \item \strong{RV_tfut}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut}
#'   }
#'   \item \strong{E_tfut}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item \strong{MEAN}: matrix \code{nFut} x \code{nEffect}
#'   }
#'   \item \strong{VARE_tfut}: variability related to the main effects (i.e.
#'   variability between the different RCMs, GCMs,..). Matrix \code{nFut} x
#'   \code{nEffect}
#'   \item \strong{CONTRIB_VARE_tfut}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{nFut} x \code{nEffect}
#'   \item \strong{listOption}: list of options used to obtained these results
#'    (obtained from \code{\link{check_option}})
#'   \item \strong{listScenario}: list of scenario factors
#'    (obtained from \code{\link{process_scen}})
#' }
#'
#' @export
#'
#' @author Guillaume Evin
lm_ANOVA = function(CCR_tfut,scen,listOption=NULL,namesEff){
  #########  process input #########
  # number of lead time
  n = dim(CCR_tfut)[2]

  # Process scenarios data.frame to get different objects
  listScenario = process_scen(scen = scen)
  nFactor = listScenario$nFactor
  nEffect = listScenario$nEffect
  listEff = listScenario$listEff

  #########  Apply lm for each time step #########
  # contrasts
  list.contrasts = list()
  for(j in 1:length(namesEff)){
    list.contrasts[[namesEff[j]]] = contr.sum
  }

  # formula
  formula = paste0("CCR_tfut ~ ",paste0(namesEff,collapse = " + "))

  # build data for the call to the lm function
  lm.data = scen

  lm.out = lm.sum = list()
  for(i in 1:n){
    lm.data$CCR_tfut=CCR_tfut[,i]
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
  Qt <- qt(pCI, lm.out[[1]]$df.residual, lower.tail = TRUE)

  #============================
  # M_tfut
  #============================
  M_tfut = list()
  M_tfut$MEAN = rep(0,n)
  M_tfut$SD = rep(0,n)
  M_tfut$CI = matrix(data=0,nrow = n, ncol = 2)

  # extract estimates for each lead time
  for(i in 1:n){
    lm.sum.i = lm.sum[[i]]
    M_tfut$MEAN[i] = lm.sum.i$coefficients[1,1]
    M_tfut$SD[i] = lm.sum.i$coefficients[1,2]
    M_tfut$CI[i,] = M_tfut$MEAN[i] + M_tfut$SD[i]*Qt
  }


  #============================
  # RV_tfut
  #============================
  # variance of residual effects, With the Bayesian estimation, this variance
  # is part of the inference (sigma^2) and we take the mean of the posterior as the point estimate
  RV_tfut = list()
  RV_tfut$MEAN = rep(0,n)

  # extract estimates for each lead time
  for(i in 1:n){
    lm.sum.i = lm.sum[[i]]
    RV_tfut$MEAN[i] = lm.sum.i$sigma^2
  }


  #============================
  # E_tfut
  #============================
  # main effects
  E_tfut = list()

  # variance of the main effects
  VARE_tfut = matrix(data=0,nrow=n,ncol=nFactor)

  # contribution to the variance of individual effects
  CONTRIB_VARE_tfut = list()

  # main effects
  for(i.eff in 1:nFactor){
    eff = namesEff[i.eff]
    #============
    # E_tfut
    #============
    # retrieve estimates for the main effects
    lHat = list()
    lHat$MEAN = matrix(data = 0, nrow = n, ncol = nEffect[i.eff])
    lHat$SD = matrix(data = 0, nrow = n, ncol = nEffect[i.eff])
    lHat$CI = array(data = 0, dim = c(n,nEffect[i.eff],2))
    for(i in 1:n){
      lm.out.i = lm.out[[i]]
      # match initial list of effects (ListEff) with the list of effects from lm
      matchEff = match(listEff[[i.eff]],listEff.LM[[eff]])

      # find coefficients from lm corresponding to this effect (without the reference effect)
      e.raw = lm.out.i$coefficients[lm.out.i$assign==i.eff]
      # compute coefficients for all effects (including the reference)
      Cmat <- ContrSumMat(as.factor(listEff.LM[[eff]]))
      e.unsorted = (Cmat %*% e.raw)[, 1]
      lHat$MEAN[i,] = e.unsorted[matchEff]

      # find standard errors
      vcov.raw = vcov(lm.out.i)
      vcov.sel = vcov.raw[lm.out.i$assign==i.eff,,drop=FALSE]
      vcov.sel = vcov.sel[,lm.out.i$assign==i.eff,drop=FALSE]
      var_bc <- Cmat %*% vcov.sel %*% t(Cmat)
      std.error.unsorted <- sqrt(diag(var_bc))
      lHat$SD[i,] = std.error.unsorted[matchEff]
      lHat$CI[i,,] = lHat$MEAN[i,] + outer(lHat$SD[i,],Qt)
    }
    E_tfut[[eff]] = lHat

    #============
    # VARE_tfut
    #============
    # predictive variance: mean of variances, which correspond to the mean of squares since the mean is 0 by constraint (see Eq 16, 17 and 18)
    VARE_tfut[,i.eff] = apply(E_tfut[[eff]]$MEAN^2,1,mean)

    #============
    # CONTRIB_VARE_tfut
    #============
    matContrib = matrix(data=0,nrow=n,ncol=nEffect[i.eff])
    AllEff2 = E_tfut[[eff]]$MEAN^2
    for(iIndEff in 1:nEffect[i.eff]){
      # individual effect squared
      indEff2 = AllEff2[,iIndEff]
      # variance of this main effect x number of factors
      infEff2Tot = (VARE_tfut[,i.eff]*nEffect[i.eff])
      # contribution of this individual effect, as a percentage
      matContrib[,iIndEff] = indEff2/infEff2Tot
    }
    CONTRIB_VARE_tfut[[eff]] = matContrib
  }


  # return results
  list.output = list(M_tfut=M_tfut,
                     RV_tfut=RV_tfut,
                     E_tfut=E_tfut,
                     VARE_tfut=VARE_tfut,
                     CONTRIB_VARE_tfut=CONTRIB_VARE_tfut,
                     listOption=listOption,
                     listScenario=listScenario)
}

#==============================================================================
# ContrSumMat
ContrSumMat <- function (fctr, sparse = FALSE) {
  if (!is.factor(fctr)) stop("'fctr' is not a factor variable!")
  N <- nlevels(fctr)
  Cmat <- contr.sum(N, sparse = sparse)
  dimnames(Cmat) <- list(levels(fctr), seq_len(N - 1))
  Cmat
}

#==============================================================================
#' Main function
#'
#' Partition uncertainty in climate responses using an ANOVA applied to climate change responses. The main
#' function of the package needs at least two arguments: \code{Y} is an ensemble of climate projections, i.e. a matrix
#' \code{nS} x \code{nY} where \code{nS} projections are provided for \code{nY} years or future time steps; 
#' \code{scen} is a \code{\link[base]{data.frame}} which provides \code{nFactor} factors for each projection.
#'
#' @param Y matrix \code{nS} x \code{nY} of climate projections.
#' @param scen data.frame \code{nS} x \code{nFactor} with the \code{nFactor} factors
#' (e.g. type of GCM) for each of the \code{nS} scenarios. The number of factors
#'  \code{nFactor} corresponds to the number of main effects that will be included in the ANOVA model.
#' @param X (optional) predictors corresponding to the projections, e.g. time or global temperature.
#' It can be a vector if the predictor is the same for all scenarios (e.g. \code{X=2001:2100}) or
#' a matrix of the same size as Y if these predictors are different for the different scenarios. By default,
#' a vector \code{1:nY} is created.
#' @param Xfut (optional) \code{nF} values of the predictor over which the ANOVA will be applied. It must be
#' a vector of values within the range of values of X. By default, it corresponds to X if X is a vector,
#' \code{1:nY} if X is \code{NULL} or a vector of 10 values equally spaced between the minimum and
#' maximum values of X if X is a matrix. The first value defines the reference period or warming level.
#' @param listOption (optional) list of options
#' \itemize{
#'    \item \strong{typeClimateResponse}: type of climate response fitted to the projections. Can be either
#' "smooth.spline" (\code{\link[stats]{smooth.spline}}) where the parameter \code{df} in (1,nY] is the 
#' equivalent number of degrees of freedom (value by default), "poly" (polynomial model fitted with \code{\link[stats]{lm}})
#' where the parameter \code{degree} is an integer indicating the degree of the polynomial (\code{\link[stats]{poly}}), 
#' "lowess" (\code{\link[stats]{lowess}}) where the parameter \code{f} in (0,1) is the proportion of points which influence
#' the smooth at each value, or "tweedie" (\code{\link[statmod]{tweedie}}), a Generalized Linear Model with the Tweedie 
#' distribution where the parameter \code{var.power} is the index of power variance function (positive value). The tweedie 
#' model is the only model that handles zeros in the climate projections, typically with \code{var.power} in (1,2).
#'   \item \strong{parClimateResponse}: parameter of the model applied for the extraction of the climate response (i.e. \code{df}, \code{degree},
#' \code{f}, or \code{var.power}). Default values are \code{df=4}, \code{degree=1}, \code{f=0.4} or \code{var.power=1.5} for the smoothing spline, 
#' polynomial, lowess, or tweedie models, respectively.
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative).
#'   \item \strong{ANOVAmethod}: ANOVA method: "lm" (default) applies a linear model (\code{\link[stats]{lm}}) to estimate the main effects. 
#' "Bayesian" applies the method described in Evin et al. (2020).   
#'   \item \strong{nBurn}: if \code{ANOVAmethod=="Bayesian"}, number of burn-in samples (default: 1000).
#'   If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: if \code{ANOVAmethod=="Bayesian"}, number of kept samples (default: 2000).
#'   If \code{nKeep} is too small, MCMC samples might not represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{probCI}: probability (in [0,1]) for the confidence intervals, \code{probCI = 0.9} by default.
#'   \item \strong{quantilePosterior}: vector of probabilities (in [0,1]) for which
#'   we compute the quantiles from the posterior distributions
#'    \code{quantilePosterior = c(0.005,0.025,0.05,0.1,0.25,0.33,0.5,0.66,0.75,0.9,0.95,0.975,0.995)} by default.
#'   \item \strong{ClimateResponse}: NULL by default. If it is provided, it must correspond to the outputs
#'   of \code{\link{fit_ClimateResponse}}, i.e. a list with \code{CR_tfut} [nS x nF], \code{CCR_tfut} [nS x nF],
#'  \code{D_t} [nS x nY], \code{CHANGEVAR_t} [nS x nY], and \code{IV_tfut} [nF].
#' }
#'
#' @return  List providing the results for each of the \code{nFut} values of \code{Xfut}:
#' \itemize{
#'   \item \strong{CR_tfut}: Fitted climate responses [nS x nF],
#'   \item \strong{CCR_tfut}: Climate change responses [nS x nF],
#'   \item \strong{D_t}: Deviation from the climate change responses [nS x nY],
#'   \item \strong{CHANGEVAR_t}: Change variable from the projections [nS x nY],
#'   \item \strong{M_tfut}: List of estimates for the grand mean:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut} of means.
#'      \item \strong{SD}: vector of length \code{nFut} of standard dev.
#'      \item \strong{CI}: matrix \code{nFut} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption}
#'      \item \strong{QUANT}: matrix \code{nFut} x \code{nQ} of quantiles of
#'      probability \code{quantilePosterior} given in \code{listOption} if
#'      \code{ANOVAmethod=="Bayesian"}.
#'   }
#'   \item \strong{E_tfut}: List of estimates for the main effects. For each
#'   main effect (GCM, RCM,..), each element of the list contains a list with:
#'   \itemize{
#'      \item \strong{MEAN}: matrix \code{nFut} x \code{nEffect}
#'      \item \strong{SD}: matrix \code{nFut} x \code{nEffect} of standard dev.
#'      \item \strong{CI}: array \code{nFut} x 2 x \code{nEffect} of credible
#'      intervals of probability \code{probCI} given in \code{listOption}
#'      \item \strong{QUANT}: array \code{nFut} x \code{nQ} x \code{nEffect} of
#'      quantiles of probability \code{quantilePosterior} given in
#'      \code{listOption} if \code{ANOVAmethod=="Bayesian"}.
#'   }
#'   \item \strong{VARE_tfut}: Matrix \code{nFut} x \code{nEffect} giving the variability
#'  related to the main effects (i.e. variability between the different RCMs, GCMs,..).
#'   \item \strong{CONTRIB_VARE_tfut}: Contribution of each individual effect
#'   to its component (percentage), e.g. what is the contribution of GCM1 to the
#'    variability related to GCMs. For each main effect (GCM, RCM,..), each
#'    element of the list contains a matrix \code{nFut} x \code{nEffect}
#'   \item \strong{RV_tfut}: List of estimates for the variance of the
#'   residual errors:
#'   \itemize{
#'      \item \strong{MEAN}: vector of length \code{nFut}.
#'      \item \strong{SD}: vector of length \code{nFut} of standard dev.
#'      if \code{ANOVAmethod=="Bayesian"}.
#'      \item \strong{CI}: matrix \code{nFut} x 2 of credible intervals of
#'      probability \code{probCI} given in \code{listOption} if
#'      \code{ANOVAmethod=="Bayesian"}.
#'      \item \strong{QUANT}: matrix \code{nFut} x \code{nQ} of quantiles of
#'      probability \code{quantilePosterior} given in \code{listOption} if
#'      \code{ANOVAmethod=="Bayesian"}.
#'   }
#'   \item \strong{IV_tfut}: Internal variability
#'   \item \strong{CCRV_tfut}: total uncertainty CCRV(t) i.e. the sum of variability related to the main effects
#'   and residual variability     
#'   \item \strong{FV_CCRV_tfut}: Fraction of total uncertainty CCRV(t) (in variance) explained by the different sources of uncertainty
#'   \item \strong{CCRV_IV_tfut}: total uncertainty CCRV(t) + IV(t), i.e. the sum of variability related to the main effects
#'   residual variability and internal variability     
#'   \item \strong{FV_CCRV_IV_tfut}: Fraction of total uncertainty CCRV(t) + IV(t) (in variance) explained by the different sources of uncertainty
#'   \item \strong{RESERR}: differences between the climate change responses and the additive anova formula (grand mean + main effects)
#'   \item \strong{Xmat}: matrix of predictors
#'   \item \strong{Xfut}: future predictor values
#'   \item \strong{namesEff}: names of the main effects
#'   \item \strong{Y}: matrix of available combinations given as inputs
#'   \item \strong{listOption}: list of options used to obtained these results
#'   (output of \code{\link{check_option}})
#'   \item \strong{listScenario}: list of scenario factors
#'   (output of \code{\link{process_scen}})
#' }
#'
#' @examples
#'
#' #____________________________________________________________
#' # CALLING QUALYPSO
#' #____________________________________________________________
#' # call main QUALYPSO function: two arguments are mandatory:
#' # - Y: Climate projections for nS scenarios and nY time steps. Y is a matrix nS x nY, we
#' # run QUALYPSO nY times, for each time step.
#' # - scen: matrix or data.frame of available combinations nS x nFactor. The number of
#' # factors nFactor corresponds to the number of main effects that will be included in the
#' # ANOVA model. In the following example, we have nFactor=2 main effects corresponding to the GCMs
#' # and RCMs.
#'
#' # Many options can be specified in the argument "listOption".
#' 
#' #____________________________________________________________
#' # Example 1: Synthetic scenarios
#' #____________________________________________________________
#' # create nS=3 fictive climate scenarios with 2 GCMs and 2 RCMs, for a period of nY=20 years
#' n=20
#' t=0:n/n
#'
#' # GCM effects (sums to 0 for each t)
#' effGCM1 = t*2
#' effGCM2 = t*1
#' effGCM3 = t*-3
#' 
#' # RCM effects (sums to 0 for each t)
#' effRCM1 = t*1
#' effRCM2 = t*-1
#' 
#' # These climate scenarios are a sum of effects and a random gaussian noise
#' scenGCM1RCM1 = effGCM1 + effRCM1 + rnorm(n=n+1,sd=0.5)
#' scenGCM1RCM2 = effGCM1 + effRCM2 + rnorm(n=n+1,sd=0.5)
#' scenGCM2RCM1 = effGCM2 + effRCM1 + rnorm(n=n+1,sd=0.5)
#' scenGCM2RCM2 = effGCM2 + effRCM2 + rnorm(n=n+1,sd=0.5)
#' scenGCM3RCM1 = effGCM3 + effRCM1 + rnorm(n=n+1,sd=0.5)
#' Y_synth = rbind(scenGCM1RCM1,scenGCM1RCM2,scenGCM2RCM1,scenGCM2RCM2,scenGCM3RCM1)
#' 
#' # Here, scen indicates that the first scenario is obtained with the combination of the
#' # GCM "GCM1" and RCM "RCM1", the second scenario is obtained with the combination of
#' # the GCM "GCM1" and RCM "RCM2" and the third scenario is obtained with the combination
#' # of the GCM "GCM2" and RCM "RCM1".
#' scen_synth = data.frame(GCM=c('GCM1','GCM1','GCM2','GCM2','GCM3'),
#' RCM=c('RCM1','RCM2','RCM1','RCM2','RCM1'))
#'
#'
#' # run QUALYPSO
#' QUALYPSO_synth = QUALYPSO(Y=Y_synth, scen=scen_synth, X=2000:2020,
#' listOption = list(typeClimateResponse="poly"))
#'
#' # plot main GCM effects
#' plotQUALYPSO_effect(QUALYPSO_synth,nameEff="GCM",xlab="Years")
#' lines(2000:2020, effGCM1,lty=2,lwd=2,col="black")
#' lines(2000:2020, effGCM2,lty=2,lwd=2,col="red")
#' lines(2000:2020, effGCM3,lty=2,lwd=2,col="green")
#'
#' # plot main RCM effects
#' plotQUALYPSO_effect(QUALYPSO_synth,nameEff="RCM",xlab="Years")
#' lines(2000:2020, effRCM1,lty=2,lwd=2,col="black")
#' lines(2000:2020, effRCM2,lty=2,lwd=2,col="red")
#'
#' # plot fraction of total variance for the differences sources of uncertainty
#' plotQUALYPSO_FV_CCRV_IV(QUALYPSO_synth,xlab="Years")
#'
#' #____________________________________________________________
#' # Example 2: climate projections of mean winter (DJF) temperature
#' # with time as the predictor
#' #____________________________________________________________
#'
#' # list of options
#' listOption = list(typeChangeVariable='abs')
#'
#' # call QUALYPSO
#' QUALYPSO_time = QUALYPSO(Y=Y_DJFTas,scen=scen_DJFTas,X=1971:2099,
#'                          Xfut=1990:2099,listOption=listOption)
#' 
#' # diagnostic of the climate responses
#' plotQUALYPSO_CR(QUALYPSO_time)
#'
#' # diagnostic of the internal variability
#' plotQUALYPSO_IV(QUALYPSO_time)
#' 
#' # grand mean effect
#' plotQUALYPSO_M(QUALYPSO_time,xlab="Years")
#'
#' # main GCM effects
#' plotQUALYPSO_effect(QUALYPSO_time,nameEff="GCM",xlab="Years")
#'
#' # main RCM effects
#' plotQUALYPSO_effect(QUALYPSO_time,nameEff="RCM",xlab="Years")
#'
#' # variance decomposition
#' plotQUALYPSO_FV_CCRV_IV(QUALYPSO_time,xlab="Years")
#'
#' #____________________________________________________________
#' # Example 3: climate projections of mean winter (DJF) temperature
#' # with global temperature as the predictor
#' #____________________________________________________________
#'
#' # list of options
#' listOption = list(typeChangeVariable='abs')
#'
#' # call QUALYPSO
#' QUALYPSO_wl = QUALYPSO(Y=Y_DJFTas,scen=scen_DJFTas,X=X_DJFTas_WL,
#'                               Xfut=seq(from=0.5,to=3,by=0.1),listOption=listOption)
#'
#' # diagnostic of the climate responses
#' plotQUALYPSO_CR(QUALYPSO_wl)
#' 
#' # diagnostic of the internal variability
#' plotQUALYPSO_IV(QUALYPSO_wl)
#' 
#' # grand mean effect
#' plotQUALYPSO_M(QUALYPSO_wl,xlab="Global warming (Celsius)")
#'
#' # main GCM effects
#' plotQUALYPSO_effect(QUALYPSO_wl,nameEff="GCM",xlab="Global warming (Celsius)")
#'
#' # main RCM effects
#' plotQUALYPSO_effect(QUALYPSO_wl,nameEff="RCM",xlab="Global warming (Celsius)")
#'
#' # variance decomposition
#' plotQUALYPSO_FV_CCRV_IV(QUALYPSO_wl,xlab="Global warming (Celsius)")
#'
#' @references Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie (2020)
#' Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation.
#' Journal of Climate. <doi:10.1175/JCLI-D-18-0606.1>.
#'
#' @export
#'
#' @author Guillaume Evin
QUALYPSO = function(Y,scen,X=NULL,Xfut=NULL,listOption=NULL){
  ######### Check inputs and assign default values ##########

  # Check list of options
  listOption = check_option(listOption)

  # Check  dimensions
  check_arguments(Y, scen, X)

  # Y is a matrix: Scenario x Time
  nS = nrow(Y)
  nY = ncol(Y)

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
    # repeat the vector to obtain a matrix
    Xmat = matrix(rep(X,nS),byrow=T,nrow=nS,ncol=nY)
  }else{
    Xmat = X
  }


  # Xfut are the values for X used to compute absolute or relative changes and indicating
  #  the future lead time / warming levels. if Xfut is provided, we check that it does
  # lead to extrapolation outside the range covered by X
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


  ##############################################
  # check presence of NAs in climate projections
  # check is some simulation chains are entirely missing
  hasAllNa = apply(Y,1,function(x) all(is.na(x)))

  if(any(hasAllNa)){
    warning(paste0('Error in QUALYPSO: some projections have only NAs in Y: ',paste(which(hasAllNa),collapse = ',')))
    return(NULL)
  }

  ##############################################
  # extract climate response
  if(is.null(listOption$ClimateResponse)){
    ClimateResponse = fit_ClimateResponse(Y,
                                      typeClimateResponse=listOption$typeClimateResponse,
                                      parClimateResponse=listOption$parClimateResponse, 
                                      Xmat=Xmat, Xfut=Xfut,
                                      typeChangeVariable=listOption$typeChangeVariable)
  }else{
    ClimateResponse = listOption$ClimateResponse
  }

  # extract quantities from these fits
  CCR_tfut = ClimateResponse$CCR_tfut
  D_t = ClimateResponse$D_t
  CHANGEVAR_t = ClimateResponse$CHANGEVAR_t
  CR_tfut = ClimateResponse$CR_tfut

  # internal variability
  IV_tfut = ClimateResponse$IV_tfut

  # CCR_tfut for the ANOVA
  CCR_tfut_ANOVA = CCR_tfut

  # names of the main effect
  if(is.null(colnames(scen))){
    namesEff = paste0("Eff",1:ncol(scen))
  }else{
    namesEff = colnames(scen)
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
  formula = paste0("CCR_tfut ~ ",paste0(namesEff,collapse = " + "))
  #lm
  lm.data = scen
  lm.data$CCR_tfut=CCR_tfut_ANOVA[,ncol(CCR_tfut_ANOVA)]
  lm.out = lm(formula, lm.data,contrasts=list.contrasts)
  if(any(is.na(lm.out$coefficients))){
    stop("singular fit encountered: the effects cannot be estimated (ill-posed problem)")
  }

  ##################
  # ANOVA on CCR_tfut
  if(listOption$ANOVAmethod=="Bayesian"){
    anova = Bayesian_ANOVA(CCR_tfut = CCR_tfut_ANOVA, scen = scen, listOption = listOption, namesEff = namesEff)
  }else{
    anova = lm_ANOVA(CCR_tfut = CCR_tfut_ANOVA, scen = scen, listOption = listOption, namesEff = namesEff)
  }



  #====================================================================================================
  # further computation using point estimates

  # first retrieve point estimates and some quantities
  effhat = anova$E_tfut
  muHat = anova$M_tfut$MEAN
  listEff = anova$listScenario$listEff
  nFactor = length(effhat)
  nP = dim(CCR_tfut_ANOVA)[2]

  # retrieve residual errors: differences between the climate change responses
  # and the additive effects (grand mean + main effect)
  mat.eff = matrix(nrow=nP,ncol=nFactor)
  RESERR = matrix(nrow=nP,ncol=nS)
  for(iS in 1:nS){
    for(iE in 1:nFactor){
      indE = which(scen[iS,iE]==listEff[[iE]])
      mat.eff[,iE] = effhat[[namesEff[iE]]]$MEAN[,indE]
    }
    RESERR[,iS] = CCR_tfut_ANOVA[iS,1:nP] - muHat - Rfast::rowsums(mat.eff)
  }

  # variance decomposition if CCRV only
  Vbind = cbind(anova$VARE_tfut,anova$RV_tfut$MEAN)

  # Total variability
  CCRV_tfut = Rfast::rowsums(Vbind)

  # Decomposition of the total uncertainty
  FV_CCRV_tfut = Vbind/replicate(n = ncol(Vbind), CCRV_tfut)
  colnames(FV_CCRV_tfut) = c(namesEff,"RV_tfut")

  # variance decomposition of CCRV + IV
  Vbind = cbind(anova$VARE_tfut,anova$RV_tfut$MEAN,IV_tfut)

  # Total variability
  CCRV_IV_tfut = Rfast::rowsums(Vbind)

  # Decomposition of the total uncertainty
  FV_CCRV_IV_tfut = Vbind/replicate(n = ncol(Vbind), CCRV_IV_tfut)
  colnames(FV_CCRV_IV_tfut) = c(namesEff,"RV_tfut","IV_tfut")


  #############################################
  # return results
  return(list(CCR_tfut=CCR_tfut,
    D_t=D_t,
    CHANGEVAR_t=CHANGEVAR_t,
    CR_tfut=CR_tfut,
    M_tfut=anova$M_tfut,
    E_tfut=anova$E_tfut,
    VARE_tfut=anova$VARE_tfut,
    CONTRIB_VARE_tfut=anova$CONTRIB_VARE_tfut,
    RV_tfut=anova$RV_tfut,
    IV_tfut=IV_tfut,
    CCRV_tfut=CCRV_tfut,
    FV_CCRV_tfut=FV_CCRV_tfut,
    CCRV_IV_tfut=CCRV_IV_tfut,
    FV_CCRV_IV_tfut=FV_CCRV_IV_tfut,
    RESERR=RESERR,
    Xmat=Xmat,Xfut=Xfut,
    namesEff=namesEff,
    Y=Y,listOption=anova$listOption,
    listScenario=anova$listScenario))
}


#==============================================================================
#' Check arguments provided to QUALYPSO
#'
#' Check type and dimensions of the arguments provided to QUALYPSO main function
#'
#' @param Y matrix \code{nS} x \code{nY}
#' @param scen data.frame \code{nS} x \code{nFactor} with the \code{nFactor} factors
#' for each of the \code{nS} scenarios. The number of factors \code{nFactor} corresponds to the 
#' number of main effects that will be included in the ANOVA model.
#' @param X (optional) predictors corresponding to the projections, e.g. time or global temperature.
#' It can be a vector if the predictor is the same for all scenarios (e.g. \code{X=2001:2100}) or
#' a matrix of the same size as Y if these predictors are different for the different scenarios. By default,
#' a vector \code{1:nY} is created.
#'
#' @author Guillaume Evin
check_arguments = function(Y, scen, X){
  # Y
  dimY = dim(Y)
  if(length(dimY)!=2){
    stop(paste('Y must be a matrix nS x nY. Y:\n', Y))
  }

  # scen
  if(!is.data.frame(scen)){
    stop(paste('scen must be a data.frame and not a ',class(scen)))
  }

  # X
  if(!is.null(X)){
    if(is.vector(X)){
      if(length(X)!=dimY[2]){
        stop(paste("if X is a vector, it must have the same length as the number of columns of Y. X:\n", 
        paste(X, collapse = ", ")))
      }
    }else if(length(dim(X))!=2){
      stop(paste("If X is not a vector, it must be a matrix nS x nY: dim(X):\n", paste(dim(X), collapse = ", ")))
    }else if(!all(dim(X)==dimY)){
      stop(paste("If X is not a vector, it must be a matrix nS x nY: dim(X):\n", paste(dim(X), collapse = ", ")))
    }
  }
}

#==============================================================================
#' Plot deviations from the climate change response and an interval representing
#' the magnitude of the estimated internal variability
#'
#' Plot D(t) as a function of X as well as +/- 1.645*sqrt(IV_tfut),
#' where IV_tfut is a function of Xfut.
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
plotQUALYPSO_IV = function(QUALYPSOOUT,lim=NULL,xlab="X",ylab="D(t)",...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut
  Xmat = QUALYPSOOUT$Xmat
  
  # retrieve D_t
  D_t = QUALYPSOOUT$D_t

  # number of scenarios
  nS = nrow(D_t)

  ylim = max(abs(range(D_t)))
  plot(-1, -1, xlim = range(c(Xmat,Xfut)), ylim = c(-ylim,ylim),
         xlab = xlab, ylab = ylab, ...)
  
  for(iS in 1:nS){
    D_t_s = D_t[iS,]
    Xs = Xmat[iS,]
    points(Xs, D_t_s,type="p",pch=20)
  }

  
    # add lines of raw projection and climate projection
    lines(Xfut, -1.644854*sqrt(QUALYPSOOUT$IV_tfut), lwd = 2, col="red")
    lines(Xfut, 1.644854*sqrt(QUALYPSOOUT$IV_tfut), lwd = 2, col="red")
  

    # add legend
    legend("bottomright",legend = c("+/- 1.645*sqrt(IV_tfut)", expression(paste(eta,"*"))),
           lty=c(1,NA), pch = c(NA,20), col=c("red","black"), bty="n")
}



#==============================================================================
#' Plot the climate response for one simulation chain
#'
#' Plot the climate response CR(t) for one simulation chain and the corresponding raw 
#' projection. The choice of the simulation chain can be given or randomly sampled.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param iS index of the projection (integer). If NULL, one scenario is chosen randomly
#' @param lim y-axis limits (default is NULL)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSO_CR = function(QUALYPSOOUT,iS=NULL,lim=NULL,xlab="X",ylab="Y",...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  CR_tfut = QUALYPSOOUT$CR_tfut

  # list of scenarios
  scen = QUALYPSOOUT$listScenario$scen

  # Xmat and Y arguments
  Xmat = QUALYPSOOUT$Xmat
  Y = QUALYPSOOUT$Y

  # number of scenarios
  nS = nrow(Y)

  # choose iS
  if(is.null(iS)){
    iS = sample(1:nS, 1)
  }

  # plot climate response for one projection
  Ys = Y[iS,]
  Xs = Xmat[iS,]
  CR_tfuts = CR_tfut[iS,]

  plot(-1, -1, xlim = range(c(Xs,Xfut)), ylim = range(c(Ys,CR_tfuts)),
        main=paste0(scen[iS,],collapse = " / "),
        xlab = xlab, ylab = ylab, ...)

  # add lines of raw projection and climate projection
  lines(Xs, Ys, lwd = 1)
  lines(Xfut, CR_tfuts, lwd = 3)

  # add legend
  legend("topleft",legend = c("Raw projection", "Climate response"),
          lty=1,lwd=c(1,3),bty="n")
}


#==============================================================================
#' Plot climate change responses for all simulation chains.
#'
#' Plot the climate change responses CCR(t) for all simulation chains.
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
plotQUALYPSO_CCR = function(QUALYPSOOUT,lim=NULL,xlab="",ylab="CCR(t)",...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  CCR_tfut = QUALYPSOOUT$CCR_tfut


  # initiate plot
  if(is.null(lim)) lim = range(CCR_tfut)
  plot(-100,-100,xlim=range(Xfut),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  for(i in 1:nrow(CCR_tfut)){
    lines(Xfut,CCR_tfut[i,],lwd=3,col=i)
  }
}

#==============================================================================
#' Plot the estimated grand mean
#'
#' Plot the estimated grand mean and the corresponding estimation uncertainty.
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
plotQUALYPSO_M = function(QUALYPSOOUT,lim=NULL,col='black',xlab="",
                                 ylab="Grand mean",addLegend=T,...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # retrieve mean
  meanPred = QUALYPSOOUT$M_tfut$MEAN

  # retrieve limits
  binf = QUALYPSOOUT$M_tfut$CI[,1]
  bsup = QUALYPSOOUT$M_tfut$CI[,2]

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
#' Plot estimated ANOVA effects for one factor
#'
#' Plot estimated ANOVA effects for one factor. By default, we plot 
#' the estimation uncertainty corresponding to a probability 0.95.
#'
#' @param QUALYPSOOUT output from \code{\link{QUALYPSO}}
#' @param nameEff name of the main effect to be plotted in \code{QUALYPSOOUT$namesEff}
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
plotQUALYPSO_effect = function(QUALYPSOOUT,nameEff,lim=NULL,
                              col=1:20,xlab="",ylab="Effect",addLegend=TRUE,
                              ...){
  # vector of predictors
  Xfut = QUALYPSOOUT$Xfut

  # index of this effect
  iEff = which(QUALYPSOOUT$namesEff==nameEff)
  if(length(iEff)==0) stop("wrong value for nameEff")

  # retrieve effects
  EffHat = QUALYPSOOUT$E_tfut[[nameEff]]
  nFactor = dim(EffHat$MEAN)[2]

  # retrieve mean
  meanPred = EffHat$MEAN

  # initiate plot
  if(is.null(lim)){
    lim = range(EffHat$CI,na.rm=TRUE)
  }
  plot(-100,-100,xlim=range(Xfut),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)

  for(i in 1:nFactor){
    # colors polygon
    colPoly = adjustcolor(col[i],alpha.f=0.2)

    # add confidence interval
    polygon(c(Xfut,rev(Xfut)),c(EffHat$CI[,i,1],rev(EffHat$CI[,i,2])),col=colPoly,lty=0)
    
    # add median
    lines(Xfut,meanPred[,i],lwd=3,col=col[i])
  }

  # legend
  if(addLegend){
    pctCI = round(QUALYPSOOUT$listOption$probCI*100)
    legend('topleft',bty='n',fill=c(NA,'black'),lwd=c(2,NA),lty=c(1,NA),
             border=c(NA,'black'),col=c('black',NA),
             legend=c('Median',paste0(pctCI,'%CI')))

    legend('bottomleft',bty='n',lwd=2,lty=1,col=col,
           legend=QUALYPSOOUT$listScenario$listEff[[iEff]])
  }
}



#==============================================================================
#' plotQUALYPSO_FV_CCRV_IV
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
plotQUALYPSO_FV_CCRV_IV = function(QUALYPSOOUT,vecEff=NULL,
                                                  col=c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1"),
                                                  xlab="",ylab="% Total Variance",addLegend=TRUE,...){
  # future predictor values
  Xfut = QUALYPSOOUT$Xfut
  nFut = length(Xfut)
  nFactor = QUALYPSOOUT$listScenario$nFactor

  # number of main effects
  if(is.null(vecEff)){
    vecEff = 1:nFactor
  }

  # Variance decomposition
  VARDECOMP = QUALYPSOOUT$FV_CCRV_IV_tfut
  VARDECOMP[,1:nFactor] = VARDECOMP[,vecEff]

  # figure
  col = col[1:(nFactor+2)]
  cum=rep(0,nFut)
  plot(-1,-1,xlim=range(Xfut),ylim=c(0,1),xaxs="i",yaxs="i",las=1,xlab=xlab,ylab=ylab,...)
  for(i in 1:(nFactor+2)){
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