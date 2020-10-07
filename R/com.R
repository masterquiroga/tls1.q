# PQ-TLS Conway-Maxwell-Poisson Distribution wrappers

# This wrappers are defined in a way to comply with R distribution standards
# in order to be used for distribution fitting via optimisation
rcom <- compoisson::rcom
qcom <-

dcom <- function(x, lambda, nu, z = NULL) {
  if (length(x) == 0)
    return(numeric(0))
  if (is.na(lambda) || is.na(nu))
    return(rep(NA, length(x)))
  if (lambda < 0 || nu < 0 || is.nan(lambda) || is.nan(nu)) 
    return(rep(NaN, length(x)))
  if (is.null(z) || is.na(z) || z < 0)
    z = compoisson::com.compute.z(lambda, nu)
  return(purrr::map_dbl(x, function(v) {
    if (is.na(v))
      return(NA)
    if (is.nan(v))
      return(NaN)
    return(compoisson::dcom(v, lambda, nu, z))
  }))
}

pcom <- function(q, lambda, nu, z = NULL) {
  if (length(q) == 0)
    return(numeric(0))
  if (is.na(lambda) || is.na(nu))
    return(rep(NA, length(q)))
  if (lambda < 0 || nu < 0 || is.nan(lambda) || is.nan(nu)) 
    return(rep(NaN, length(q)))
  if (is.null(z) || is.na(z) || z < 0)
    z = compoisson::com.compute.z(lambda, nu)
  return(purrr::map_dbl(q, function(v) {
    if (is.na(v))
      return(NA)
    if (is.nan(v))
      return(NaN)
    return(min(sum(dcom(0:floor(v), lambda, nu, z)), 1))
  }))
}



dtweedie <- function(x, xi, mu, phi) tweedie::dtweedie(x, xi, mu, phi)
ptweedie <- function(q, xi, mu, phi) tweedie::ptweedie(q, xi, mu, phi)
qtweedie <- function(p, xi, mu, phi) tweedie::dtweedie(p, xi, mu, phi)
rtweedie <- function(n, xi, mu, phi) tweedie::dtweedie(n, xi, mu, phi)