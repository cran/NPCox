#' This is some description of this function.
#' @title Nonparametric and semiparametric Cox regression model.
#' @description Plot of time-varying coefficient estimated through function 'npcox' or 'spcox'.
#' @details This is the plot function for 
#' @param temp Estimation result from function 'npcox' or 'spcox'.
#' @param xrange Illustration range for x-axis.
#' @param CIlevel The default confidence level stands at 0.95.
#' @return The plot of nonparametric coefficients function 'npcox' or 'spcox'.
#' @export
#' @examples 
#' data(pbc)
#' pbc = pbc[(pbc$time < 3000) & (pbc$time > 800), ] 
#' Z   = pbc[,c("age","edema")]
#' colnames(Z) = c("age","edema")
#' del = pbc$status
#' tim = pbc$time
#' res = npcox(cva = Z,delta = del, obstime = tim, bandwidth = 500)
#' op  = par(mfrow = c(1,2))
#' npplot(res)
#' par(op)


npplot = function(temp, xrange = NULL, CIlevel = 0.95){
  # print("Please decide the partition of a figure with par(...)")
  temp1 = temp
  ind = sort(unique(c(which(is.na(temp$temporal_coef[,1])),which(is.na(temp$temporal_coef[,1])))))
  if(length(ind)>0){
    temp1$temporal_coef = temp$temporal_coef[-ind,]
    temp1$temporal_coef_SEE = temp$temporal_coef_SEE[-ind,]
    temp1$data = temp$data[-ind,]
    temp  = temp1
  }
  
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
  if(!is.null(xrange)){
    hmin = max(which(obs<xrange[1]))
    hmax = min(which(obs>xrange[2])) - 1
  }
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
    lines(obs[cc1], uppci, col = 'blue')
    lines(obs[cc1], lowci, col = 'blue')
    abline(h = 0, col = "red")
  }
}


