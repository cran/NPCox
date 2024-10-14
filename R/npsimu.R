#' This is some description of this function.
#' @title Simulation data generation of PH model with time-varying coefficients.
#' @description An example of data generation for nonparametric PH model.
#' @details 'npsimu' is designed for PH model with time-varying coefficients, h(t) = h0(t)exp(b(t)'Z), generating the covariates, observed time and censoring indicator.
#' @param n The number of sample size, which can be self-determined.
#' @param cenpara Censoring parameter, which is supposed to be positive, for adjustment of censoring rate.
#' @return a list that contain covariates, observed time and censoring indicator.
#' @examples
#' data = npsimu(200)

npsimu = function(n, cenpara = NULL){
  
  beta1 = function(t){t^2}
  beta2 = function(t){1-t}
  
  cva  = data.frame(Z1 = runif(n,0.01,1), Z2 = runif(n,0.01,2))
  Z    = as.matrix(cva)
  if(is.null(cenpara)){ cenpara = 4 }
  
  cen   = runif(n,0,cenpara)
  Ti    = array();X = array()
  for(i in 1:n){
    tem   = log(runif(1,0,1))
    lamt  = function(t){  exp(beta1(t)*Z[i,1] + Z[i,2]*beta2(t) - 1)  }
    Ft    = function(t){  integrate(lamt,0,t)$value + tem  }
    Ti[i] = ifelse(Ft(0)*Ft(10)<0,unlist(uniroot(Ft,c(0,10),tol = 1e-4)$root),10)
    X[i]  = min(Ti[i],cen[i])
  }
  delta   = as.numeric(Ti < cen)
  cenrate = 1 - sum(delta)/n
  
  data = list(cva = cva, delta = delta, obstime = X)
  return(data)
}