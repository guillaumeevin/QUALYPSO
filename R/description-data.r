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
#' @name X_globaltas
#' @docType data
#' @usage data(X_globaltas)
#' @format matrix 20 scenarios x 129 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"X_globaltas"


#' Vector of of future warming levels
#'
#' Equally spaced vector of of future warming levels
#'
#' @name Xfut_globaltas
#' @docType data
#' @usage data(Xfut_globaltas)
#' @format vector of length 13
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"Xfut_globaltas"


#' Years 1971-2099 repeated for the 20 scenarios
#'
#' @name X_time_mat
#' @docType data
#' @usage data(X_time_mat)
#' @format matrix 20 scenarios x 129 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"X_time_mat"


#' X_time_vec gives the years corr. to Y, i.e. from 1971 to 2099
#'
#' @name X_time_vec
#' @docType data
#' @usage data(X_time_vec)
#' @format vector of length 129
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"X_time_vec"


#' Xfut_time is a vector of 11 years equally spaced from 1999 to 2099
#'
#' @name Xfut_time
#' @docType data
#' @usage data(Xfut_time)
#' @format vectors of length 11
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"Xfut_time"


#' List of GCM and RCM which have been used for the 20 climate projections
#'
#' scenAvail gives the GCM and RCM which have been used for the 20 climate
#' projections (obtained with the RCP8.5)
#'
#' @name scenAvail
#' @docType data
#' @usage data(scenAvail)
#' @format data.frame with 20 rows and two columns: GCM and RCM
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"scenAvail"


#' Mean winter temperature over CEU with 20 GCM/RCM combinations for 1971-2099
#'
#' climate projections of mean winter (DJF) temperature over the SREX region CEU
#' simulated by 20 combinations of CMIP5 GCMs and RCMs for the period 1971-2099
#'
#' @name Y
#' @docType data
#' @usage data(Y)
#' @format matrix 20 scenarios x 129 years
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Seneviratne, S. I. et al. Changes in Climate Extremes and their
#' Impacts on the Natural Physical Environment, in: Managing the Risks of
#' Extreme Events and Disasters to Advance Climate Change Adaptation: Special
#' Report of the Intergovernmental Panel on Climate Change, edited by: Field,
#' C., Barros, V., Stocker, T., and Dahe, Q., Cambridge University Press,
#' Cambridge, 109-230, https://doi.org/10.1017/CBO9781139177245.006, 2012
"Y"
