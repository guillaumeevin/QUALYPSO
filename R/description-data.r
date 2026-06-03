#' Annual warming levels simulated by different CMIP5 GCMs
#'
#' Annual warming levels at the planetary scales simulated by different CMIP5
#' GCMs for the period 1971-2099. Warming levels are obtained with respect to
#' the year 1860 (common starting year of the CMIP5 simulations). These warming
#' levels have been obtained with the following steps:
#' \enumerate{
#' \item Annual tas averages simulated by different CMIP5 have first been smoothed
#' using a smoothing spline. Let us denote these smoothed values by
#' tas_GCM(y) for a year y.
#' \item Large discrepancies can be observed for tas_GCM_smooth(y) even in the past
#' due to large first-order biases in the GCM simulations. In order to obtain a
#' common reference, we also consider observed tas estimates at the global scale.
#' HadCRUT5 (Morice et al., 2021, 10.1029/2019JD032361) provides anomalies with
#' respect to the period 1961-1990. An estimate of absolute average temperature
#' for this period is 14°C (Jones et al., 1999, 10.1029/1999RG900002). Smoothed
#' estimates of absolute tas averages are obtained using a smoothing spline and
#' is denoted by tas_obs(y).
#' \item Warming levels are obtained as anomalies with respect to the period 1860
#' and considering a reference year, here 1990, where the warming levels WL are
#' in agreement:
#' WL(y) = tas_GCM(y)-tas_GCM(1990)+tas_obs(1990)-tas_obs(1860)
#' }
#'
#' @name X_DJFTas_WL
#' @docType data
#' @usage data(X_DJFTas_WL)
#' @format matrix 20 scenarios x 129 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"X_DJFTas_WL"

#' data.frame indicating which GCMs and RCMs have been used for the 20 climate projections
#' of mean winter temperature over CEU 
#'
#' scen_DJFTas gives the GCM and RCM which have been used for the 20 climate
#' projections (obtained with the RCP8.5)
#'
#' @name scen_DJFTas
#' @docType data
#' @usage data(scen_DJFTas)
#' @format data.frame with 20 rows and two columns: GCM and RCM
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"scen_DJFTas"

#' Mean winter temperature over CEU with 20 GCM/RCM combinations for 1971-2099
#'
#' climate projections of mean winter (DJF) temperature over the SREX region CEU
#' simulated by 20 combinations of CMIP5 GCMs and RCMs for the period 1971-2099
#'
#' @name Y_DJFTas
#' @docType data
#' @usage data(Y_DJFTas)
#' @format matrix 20 scenarios x 129 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Seneviratne, S. I. et al. Changes in Climate Extremes and their
#' Impacts on the Natural Physical Environment, in: Managing the Risks of
#' Extreme Events and Disasters to Advance Climate Change Adaptation: Special
#' Report of the Intergovernmental Panel on Climate Change, edited by: Field,
#' C., Barros, V., Stocker, T., and Dahe, Q., Cambridge University Press,
#' Cambridge, 109-230, https://doi.org/10.1017/CBO9781139177245.006, 2012
"Y_DJFTas"

#' Annual warming levels simulated by different CMIP5 GCMs corresponding to the
#' 18 projections Y_SWE
#'
#' Annual warming levels at the planetary scales simulated by different CMIP5
#' GCMs for the period 1951-2099. Warming levels are obtained with respect to
#' the year 1860 (common starting year of the CMIP5 simulations). See Evin et al. (2025)
#' for further details.
#'
#' @name X_SWE_WL
#' @docType data
#' @usage data(X_SWE_WL)
#' @format matrix 18 scenarios x 149 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Evin, G., E. Le Roux, E. Kamir, and S. Morin. « Estimating changes in 
#' extreme snow load in Europe as a function of global warming levels ». Cold Regions
#'  Science and Technology 231 (2025): 104424. https://doi.org/10.1016/j.coldregions.2025.104424.
"X_SWE_WL"

#' Annual maxima of snow water equivalent for Loire-Atlantique, France, a NUTS-3 region located
#' at a low mean elevation (0 m) and has a suboceanic climate with quite mild and rainy winters.
#' SWE maxima are provided for 18 projections obtained with 9 different combinations of GCMs and RCMs
#' and two emission scenarios (RCP4.5 and RCP8.5).
#'
#' @name Y_SWE
#' @docType data
#' @usage data(Y_SWE)
#' @format matrix 18 scenarios x 149 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Evin, G., E. Le Roux, E. Kamir, and S. Morin. « Estimating changes in 
#' extreme snow load in Europe as a function of global warming levels ». Cold Regions
#'  Science and Technology 231 (2025): 104424. https://doi.org/10.1016/j.coldregions.2025.104424.
"Y_SWE"

#' data.frame indicating which GCM, RCM and RCP scenarios have been used to produce the 18
#' projections of SWE maxima.
#' 
#' SWE maxima are provided for 18 projections obtained with 9 different combinations of GCMs and RCMs
#' and two emission scenarios (RCP4.5 and RCP8.5).
#'
#' @name scen_SWE
#' @docType data
#' @usage data(scen_SWE)
#' @format matrix 18 scenarios x 149 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Evin, G., E. Le Roux, E. Kamir, and S. Morin. « Estimating changes in 
#' extreme snow load in Europe as a function of global warming levels ». Cold Regions
#'  Science and Technology 231 (2025): 104424. https://doi.org/10.1016/j.coldregions.2025.104424.
"scen_SWE"