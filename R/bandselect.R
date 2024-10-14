#' This is some description of this function.
#' @title Bandwidth selection function for 'npcox' and 'spcox'.
#' @description The embeded bandwidth selection function in 'npcox' and 'spcox' together with prediction error calculation. 
#' @details 'bandwidth_selection' function can provide the prediction error calculation with given bandwidth, or produce the optimal bandwidth based on an arithmetic progression of bandwidths.
#' @param cva Covariate Z in h(t) = h0(t)exp(b(t)'Z)
#' @param delta Right censoring indicator for the model
#' @param obstime The observed time = min(censoring time, observed failure time)
#' @param bandwidth Bandwidth for kernel function, which can be specified. The default value is FALSE and can be selected through least prediction error over all subjects.
#' @return Prediction error on given bandwidth or list that contains the optimal bandwidth, an arithmetic progression of bandwidths with corresponding prediction error.
#' @export
#' @examples
#' data(pbc)
#' dta  = na.omit(pbc[,c('time', 'status', 'age', "edema")])
#' dta  = dta[dta$time>600 & dta$time<2000,]
#' dta[,'status'] = sign(dta[,'status'])
#' colnames(dta) = c('time', 'status', 'age',  "edema")
#' res = bandselect(cva = dta[,3:4], delta = dta$status, obstime = dta$time, bandwidth = 700)


bandselect = function(cva, delta, obstime, bandwidth = FALSE){
  
  Z    = cva
  del1 = length(unique(delta))
  
  if(del1 > 2) stop("The delta input contains more than two elements, please check", call. = FALSE)
  if(sum(class(Z) != c('matrix')) == 0) stop("Please transform covariate into matrix with specified variable name", call. = FALSE)
  if(sum( nchar(colnames(Z)) == 0 )){   stop("Please transform covariate into matrix with specified variable name", call. = FALSE) }
  covname = colnames(cva)
  
  # PH part: β(t) estimation
  Kh   = function(h,x){ifelse(abs(x/h)<1, 0.75*(1-(x/h)^2),0)}
  
  PE   = function(bhat,Zs,ds,Xs){
    n    = floor(nrow(bhat)*0.9)
    bhat = bhat[1:n,]
    Zs   = Zs[1:n,]
    ds   = ds[1:n]
    Xs   = Xs[1:n]
    # prediction error calculation, refer to Lu tian(2005)
    bz  = array(0, dim = c(n,n))
    for(i in 1:n){
      for(j in 1:n){
        bz[i,j] = sum(bhat[i,]*Zs[j,])
      }
    }
    pie = array(0, n)
    for(l in 1:n){
      pie[l] = -( bz[l,l] - log( sum(exp(bz[l,l:n])) ) )
    }
    return(sum(pie))
  }
  
  esti = function(h,Zs,ds,Xs){
    
    cri  = 0.001
    nit  = 100 # number of iteration
    conv = array(0,n)
    
    bhat = array(0,dim = c(n,r))
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i) # ite = 2
      
      for(ite in 2:nit){
        
        beta = brec[ite-1,]
        
        G = array(n); Gpie = array(n); # i = n
        for(i in n:min_i){
          Gpie[i] = exp(c(beta%*%Zs[i,]))
          if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] }
        }
        Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
        for(i in n:min_i){
          Gdpie[i,] = Gpie[i]*Zs[i,]
          if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }
        }
        ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
        for(i in non0_i){
          ldpie[i,] = ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]);
          if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] };
        }
        Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
        for(j in n:min_i){
          Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
          if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
        }
        ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
        for(i in non0_i){
          lddpie[i,,] = -ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
          if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
        }
        
        temp1 = ifelse(r>1, det(ldd[max_i,,]), ldd[max_i,,])
        
        if(abs(temp1)>1e-8){
          temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
          brec[ite,] = beta - temp
        }else{
          brec[ite,] = beta + 0.01
        }
        
        if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
        if(max(abs(brec[ite,])) > 10){ conv[t] = 1; brec[ite,] = NA; break }
      }
      bhat[t,] = brec[ite,]
    }
    for(i in 1:(hmin-1)){ bhat[i,] = bhat[hmin,] }
    for(i in (hmax+1):n){ bhat[i,] = bhat[hmax,] }
    return(list(bhat = bhat, conv = conv))
    
  }
  
  # Rename of covar and obstime
  X   = obstime
  
  # Some constants
  n   = nrow(Z)
  r   = ncol(Z)
  n1  = length(delta)
  n2  = length(X)
  ind = sum(is.na(Z)) + sum(is.na(delta)) + sum(is.na(X))
  
  # sort order for X
  XZ  = as.matrix(cbind(X,Z,delta))
  sor = XZ[order(XZ[,1]),]
  Xs  = as.numeric(sor[,1])
  Zs  = as.matrix(sor[,(1:r)+1])
  ds  = as.numeric(sor[,r+2])
  
  if(n != n1 || n != n2 || n1 != n2){
    return("Incorrect length of covariates, censoring indicator or observed time")
  }else if(ind > 0){
    return("data contains NA's")
  }else{
    
    if(bandwidth){
      h    = bandwidth
      bhat = esti(h,Zs,ds,Xs)
      if(sum(bhat$conv)){
        for(i in 1:n){
          if(sum(is.na(bhat$bhat[i,]))){
            bhat$bhat[i,] = bhat$bhat[i-1,]
          }
        }
      }
      pe   = PE(bhat$bhat,Zs,ds,Xs)
      
      bandPE = cbind(h, pe)
      colnames(bandPE) = c('Bandwidth', 'Prediction error')
      
      return(bandPE)
    }else{
      
      pb = progress_bar$new(total = 51)
      diff  = Xs
      for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1] }
      diff  = sort(diff, decreasing = T)
      hmi   = as.numeric(quantile(Xs,0.05))
      sep   = (as.numeric(quantile(Xs,0.2)) - hmi)/50
      ha    = seq(hmi, as.numeric(quantile(Xs,0.2)), by = sep)
      PErec = array(0, length(ha))
      for(l in 1:length(ha)){
        pb$tick()
        h    = ha[l]
        bhat = esti(h,Zs,ds,Xs)
        PErec[l] = PE(bhat$bhat,Zs,ds,Xs)
        # print(paste(Sys.time(), '  ', l ))
      }
      
      # Bandwidth selection
      h = ha[which(PErec == min(PErec[which(!is.na(PErec))]))]
      bandPE = cbind(ha, PErec)
      colnames(bandPE) = c('Bandwidth', 'Prediction error')
      # print(paste("Note: The selected bandwidth is ", h, sep = ""))
      lis = list('band-pred' = bandPE, 'optimal bandwidth' = h)
      if(length(h) == 0){
        stop('No available bandwidth can be provided, please specify a bandwidth')
      }else{
        return(lis)
      }
    }
  }
}

