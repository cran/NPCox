#' This is some description of this function.
#' @title Nonparametric and semiparametric Cox regression model.
#' @description Plot of time-varying coefficient estimated through function 'npcox' or 'spcox'.
#' @details This is the plot function for 
#' @param temp Estimation result from function 'npcox' or 'spcox'.
#' @param CIlevel The default confidence level stands at 0.95.
#' @return The plot of nonparametric coefficients function 'npcox' or 'spcox'.
#' @export
#' @examples 
#' data(pbc)
#' pbc = pbc[(pbc$time < 3000) & (pbc$time > 1500), ] 
#' Z   = pbc[,c("age","edema")]
#' colnames(Z) = c("age","edema")
#' del = pbc$status
#' tim = pbc$time
#' res = npcox(cva = Z,delta = del, obstime = tim, bandwidth = 500)
#' op  = par(mfrow = c(1,2))
#' npplot(res)
#' par(op)


npplot = function(temp, CIlevel = 0.95){
  # print("Please decide the partition of a figure with par(...)")
  coef = as.matrix(temp$temporal_coef)
  see  = as.matrix(temp$temporal_coef_SEE)
  data = as.data.frame(temp$data)
  obs  = data$obstime
  h    = temp$bandwidth
  hmin = max(which(obs<min(obs)+h))
  hmax = min(which(obs>max(obs)-h)) - 1
  
  CIquan = qnorm(1-(1-CIlevel)/2)
  
  r    = ncol(coef)
  n    = nrow(coef)
  cc1  = c(hmin:hmax)
  
  for(j in 1:r){
    uppci = coef[cc1,j]+CIquan*see[cc1,j]
    lowci = coef[cc1,j]-CIquan*see[cc1,j]
    if(sum(is.na(see)) == n*r){
      upp   = max(coef[cc1,j]) + 0.2*abs(max(coef[cc1,j]))
      low   = min(coef[cc1,j]) - 0.2*abs(min(coef[cc1,j]))
    }else{
      upp   = max(uppci) + 0.2*max(uppci)
      low   = min(lowci) - 0.2*min(lowci)
    }
    upp = max(upp,  abs(upp)*0.1)
    low = min(low, -abs(low)*0.1)
    plot(obs[cc1], coef[cc1,j],type = "l",main = colnames(coef)[j], xlab = "Time",
         ylab = "Time-varying coefficients", ylim = c(low, upp))
    lines(obs[cc1], uppci)
    lines(obs[cc1], lowci)
    abline(h = 0, col = "red")
  }
}

