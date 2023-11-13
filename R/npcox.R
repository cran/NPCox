#' This is some description of this function.
#' @title Nonparametric and semiparametric Cox regression model.
#' @description Estimation of proportional hazards (PH) model with time-varying coefficients.
#' @details 'npcox' function is designed for PH model with time-varying coefficients, h(t) = h0(t)exp(b(t)'Z), providing estimation of b(t) and its pointwise standard errors on [bandwidth, max(obstime)-badwidth]. 
#' @param cva Covariate Z in h(t) = h0(t)exp(b(t)'Z)
#' @param delta Right censoring indicator for the model
#' @param obstime The observed time = min(censoring time, observed failure time)
#' @param SE Whether or not the estimation of standard error through resampling method will be done. The default value is FALSE.
#' @param bandwidth Bandwidth for kernel function, which can be specified. The default value is FALSE and can be selected through least prediction error over all subjects.
#' @param resamp Number of resampling for estimation of pointwise standard error. The default value is 100.
#' @return a list that contain the estimation result of temporal coefficients, standard error estimation, selected or predesigned bandwidth, dataset, unconverged time points.
#' @export
#' @examples
#' data(pbc)
#' pbc = pbc[(pbc$time < 3000) & (pbc$time > 800), ] 
#' Z   = pbc[,c("age","edema")]
#' colnames(Z) = c("age","edema")
#' del = pbc$status
#' tim = pbc$time
#' res = npcox(cva = Z,delta = del, obstime = tim, bandwidth = 500)



npcox = function(cva, delta, obstime, SE = FALSE, bandwidth = FALSE, resamp = 100){
  
  Z   = cva
  if(sum(class(Z) != c('matrix')) == 0) stop("Please transform covariate into matrix with specified variable name", call. = FALSE)
  if(sum( nchar(colnames(Z)) == 0 )){   stop("Please transform covariate into matrix with specified variable name", call. = FALSE) }
  covname = colnames(cva)
  
  # PH part: β(t) estimation
  Kh   = function(h,x){ifelse(abs(x/h)<1, 0.75*(1-(x/h)^2),0)}
  
  PE   = function(bhat,Zs,ds,Xs){
    # prediction error calculation, refer to Lutian(2005)
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
  
  band = function(){
    
    if(bandwidth){
      h = bandwidth
    }else{
      diff  = Xs
      for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1] }
      diff  = sort(diff,decreasing = T)
      hmi   = as.numeric(quantile(Xs,0.05))
      sep   = (as.numeric(quantile(Xs,0.2)) - hmi)/50
      ha    = seq(hmi, as.numeric(quantile(Xs,0.2)), by = sep)
      PErec = array(0, length(ha))
      for(l in 1:length(ha)){
        h    = ha[l]
        bhat = esti(h,Zs,ds,Xs)
        PErec[l] = PE(bhat$bhat,Zs,ds,Xs)
        # print(paste(Sys.time(), '  ', l ))
      }
      # Bandwidth selection
      h = ha[which(PErec == min(PErec[which(!is.na(PErec))]))]
      # cbind(ha, PErec)
      print(paste("Note: The selected bandwidth is ", h, sep = ""))
    }
    return(h)
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
      max_i = max(non0_i)
      
      for(ite in 2:nit){
        
        beta = brec[ite-1,]
        
        G = array(n); Gpie = array(n)
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
        
        if(abs(det(ldd[max_i,,]))>1e-8){
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
  
  esti_sd = function(h,Zs,ds,Xs){
    
    cri   = 0.001
    nit   = 100
    
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    bsd  = array(0,dim = c(n,r))
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      
      betaM = array(0,dim = c(M,r))
      conv2 = array(0,M)
      for(k in 1:M){
        Gi = rexp(n)
        for(ite in 2:nit){
          
          beta = brec[ite-1,]
          
          G = array(n); Gpie = array(n)
          for(i in n:min_i){
            Gpie[i] = exp(c(beta%*%Zs[i,]))
            if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] } # View(cbind(G,Gpie))
          }
          Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
          for(i in n:min_i){
            Gdpie[i,] = Gpie[i]*Zs[i,]
            if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }; # View(cbind(Gd,Gdpie))
          }
          ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
          for(i in non0_i){
            ldpie[i,] =  Gi[i]*ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]); # View(cbind(ldpie,ld))
            if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] }; # View(cbind(ldpie,ds,Ker))
          }
          Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
          for(j in n:min_i){
            Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
            if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
          }
          ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
          for(i in non0_i){
            lddpie[i,,] = -Gi[i]*ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
            if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
          }
          
          if(abs(det(ldd[max_i,,]))>1e-8){
            temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
            brec[ite,] = beta - temp
          }else{
            brec[ite,] = beta + 0.01
          }
          
          if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
          if(max(abs(brec[ite,]))>10){ conv2[k] = 1; break }
        }
        betaM[k,] = beta
      }
      if(sum(conv2) > 0 & sum(conv2) < (M-1) ){
        bsd[t,] = apply(betaM[-which(conv2==1),], 2, sd)
      }else if(sum(conv2) >= (M-1)){
        bsd[t,] = NA
      }else{
        bsd[t,] = apply(betaM, 2, sd)
      }
    }
    
    for(i in 1:(hmin-1)){ bsd[i,] = bsd[hmin,] }
    for(i in (hmax+1):n){ bsd[i,] = bsd[hmax,] }
    return(bsd)
  }
  
  # Rename of covar and obstime
  X   = obstime
  M   = resamp # if(!resamp){M = 100}else{}
  
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
    h    = band()
    if(length(h) == 0){
      stop('No available bandwidth can be provided, please specify a bandwidth')
    }else{
      tem1 = esti(h,Zs,ds,Xs)
      bhat = tem1$bhat
      conv = tem1$conv
      
      if(SE){
        #print("Commencing estimation of standard error")
        bsd = esti_sd(h,Zs,ds,Xs)
      }else{
        bsd  = array(NA, dim = c(nrow(bhat), ncol(bhat)) )
      }
      
      ind     = conv
      convind = which(ind == 1)
      if(length(convind) > 0){
        bhat[convind,] = NA
        bsd[convind,]  = NA
      }
      data = cbind(Zs, Xs, ds)
      colnames(data) = c(covname,'obstime','delta')
      colnames(bhat) = covname
      colnames(bsd)  = paste(covname,"_SEE", sep = "")
      
      return(list(temporal_coef = bhat, temporal_coef_SEE = bsd, bandwidth = h, data = data,
                  inconverged = convind,
                  "Time points (not converged)" = paste('obstime = ', Xs[convind], sep = "")))
    }
  }
}


