#' This is some description of this function.
#' @title Nonparametric and semiparametric Cox regression model.
#' @description Estimation of proportional hazards (PH) model with time-varying coefficients and constant coefficients.
#' @details 'spcox' is designed for PH model with both time-varying and constant coefficients, h(t) = h0(t)exp(b(t)'Z1 + c*Z2), providing estimation of b(t), c and their standard errors.
#' @param cva_cons Covariate Z1 with constant coefficeint c in h(t) = h0(t)exp(c'Z1 + b(t)'Z2)
#' @param cva_time Covariate Z2 with time-varying coefficeint b(t) in h(t) = h0(t)exp(c'Z1 + b(t)'Z2)
#' @param delta Right censoring indicator for the model
#' @param obstime The observed time = min(censoring time, observed failure time)
#' @param SE Whether or not the estimation of standard error through resampling method will be done. The default value is FALSE.
#' @param bandwidth Bandwidth for kernel function, which can be specified. The default value is FALSE and can be selected through least prediction error over all subjects.
#' @param resamp Number of resampling for estimation of pointwise standard error. The default value is 100.
#' @return a list that contain the estimation result of both temporal (hat{b}(t)) and constant (hat{c}) coefficients, standard error estimation, selected/predesigned bandwidth, dataset, unconverged time points.
#' @export
#' @examples 
#' data(pbc)
#' pbc = pbc[(pbc$time < 3000) & (pbc$time > 1500), ] 
#' Z1  = as.matrix(pbc[,5])
#' Z2  = as.matrix(pbc[,c('albumin')])
#' colnames(Z1) = c('age')
#' colnames(Z2) = c('albumin')
#' del = pbc$status
#' tim = pbc$time
#' res1 = spcox(cva_cons = Z1, cva_time = Z2, delta = del, obstime = tim, bandwidth = 500)

spcox = function(cva_cons, cva_time, delta, obstime, SE = FALSE, bandwidth = FALSE, resamp = 100){
  
  # PH part: β(t) estimation
  # Bandwidth assignment or selection function
  Kh  = function(h,x){ifelse(abs(x/h)<1, 0.75*(1-(x/h)^2),0)}
  fm2 = function(x){x = array(x,dim = c(length(x),1)); return(x%*%t(x))}
  hmx = function(){ min(which(Xs>max(obstime)-h)) }# floor(.9*n)min(which(Xs>max(Xs)-h)) - 1
  
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
      # print("--------------------------------------------------------------------")
      # print(paste("Note: The bandwidth is predesigned as ", bandwidth, sep = ""))
      h = bandwidth
    }else{
      # Prediction error for bandwidth selection
      # print("--------------------------------------------------------------------")
      # print("Note: No predesigned bandwidth, commencing bandwidth selection procedure")
      
      diff  = Xs
      for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1] }
      diff  = sort(diff,decreasing = T)
      hmi   = diff[5]
      sep   = (Xs[n]/3 - Xs[1])/30
      ha    = seq(max(hmi,Xs[1]), Xs[n]/3, by = sep)
      PErec = array(0, length(ha))
      for(l in 1:length(ha)){
        h    = ha[l]
        #print(paste("Calculating PE for h = ",h, sep = ""))
        bhat = esti(h,Zs,ds,Xs)
        PErec[l] = PE(bhat,Zs,ds,Xs)
      }
      # Bandwidth selection
      h = ha[which(PErec == min(PErec))]
      # print(paste("Note: The selected bandwidth is ", h, sep = ""))
    }
    return(h)
  }
  
  cons = function(bhat,h,Zs,ds,Xs){
    bhatcheck = sum(is.na(bhat))
    if(bhatcheck > 0){
      ind1 = which(is.na(bhat[,1]))
      for(l in ind1){
        bhat[l,] = bhat[l-1,]
      }
    }
    
    # esti of the partly constant effect
    S0 = array(0, n)
    S1 = array(0, dim = c(n,r))
    S2 = array(0, dim = c(n,r,r))
    bz = array(0, dim = c(n,n))
    for(j in 1:n){
      for(i in 1:n){
        bz[j,i] = sum(bhat[j,]*Zs[i,])
      }
    }
    for(j in 1:n){
      S0[j] = sum(exp(bz[j,j:n]))
      if(j == n){
        S1[j,]  = exp(bz[j,j:n])*Zs[n,]
        S2[j,,] = exp(bz[j,j:n])*fm2(Zs[n,])
      }else{
        for(i in j:n){
          S1[j,]  = S1[j,]  + exp(bz[j,i])*Zs[i,]
          S2[j,,] = S2[j,,] + exp(bz[j,i])*fm2(Zs[i,])
        }
      }
    }
    Vbt = array(0, dim = c(n,r,r))
    for(j in 1:n){ Vbt[j,,] = S2[j,,]/S0[j] - fm2(S1[j,]/S0[j] ) }
    
    hmin = max(which(Xs<min(Xs)+h))
    hmax = hmx()
    Ibt  = array(0, dim = c(n,r,r))
    Iinv = array(0, dim = c(n,r,r))
    for(j in hmin:hmax){
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[j]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      for(i in min_i:max_i){
        Ibt[j,,] = Ibt[j,,] + Vbt[i,,]*Ker[i]*ds[i]/n
      }
      Iinv[j,,] = solve(Ibt[j,,])
    }
    Jt    = array(0, dim = c(n,r1,r1))
    for(i in hmin:hmax){ Jt[i,,] = solve(Iinv[i,1:r1,1:r1]) }# v(c(Jt))
    diff  = Xs
    for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1]}
    Jtint = array(0,dim = c(r1,r1))
    for(i in hmin:hmax){ Jtint = Jtint + Jt[i,,]*diff[i] }
    wop   = array(0, dim = c(n,r1,r1))
    Jtint_inv = solve(Jtint)
    for(i in hmin:hmax){ wop[i,,] = Jtint_inv%*%Jt[i,,] }
    beta1 = rep(0,r1)
    for(i in hmin:hmax){ beta1 = beta1 + wop[i,,]%*%bhat[i,1:r1]*diff[i]}
    beta1_sd = sqrt(diag(Jtint_inv)) # /n
    
    return(data.frame(constant_coef = beta1, constant_coef_SEE = beta1_sd))
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
  # covname = c(colnames(cva_cons), colnames(cva_time))
  Z1  = cva_cons
  Z2  = cva_time
  if(sum(class(Z1) != c('matrix')) == 0) stop("Please transform cva_cons into matrix with specified variable name", call. = FALSE)
  if(sum(class(Z2) != c('matrix')) == 0) stop("Please transform cva_time into matrix with specified variable name", call. = FALSE)
  if(sum( nchar(colnames(Z1)) == 0 )){   stop("Please transform cva_cons into matrix with specified variable name", call. = FALSE)  }
  if(sum( nchar(colnames(Z2)) == 0 )){   stop("Please transform cva_time into matrix with specified variable name", call. = FALSE)  }
  X   = obstime
  M   = resamp
  
  # Some constants
  Z   = cbind(Z1,Z2)
  covname = colnames(Z)
  
  n   = nrow(Z)
  r1  = ncol(Z1)
  r2  = ncol(Z2)
  r   = r1 + r2
  n1  = length(delta)
  n2  = length(X)
  ind = sum(is.na(Z)) + sum(is.na(delta)) + sum(is.na(X))
  
  # sort order for X
  XZ  = cbind(X,Z,delta)
  sor = XZ[order(XZ[,1]),]
  Xs  = sor[,1]
  Zs  = sor[,(1:r)+1]
  ds  = sor[,r+2]
  
  if(n != n1 || n != n2 || n1 != n2){
    return("Incorrect length of covariates, censoring indicator or observed time")
  }else if(ind > 0){
    return("data contains NA's")
  }else{
    h = band()
    #print("Commencing estimation of temporal coefficient")
    tem1 = esti(h,Zs,ds,Xs)
    bhat = tem1$bhat
    convind = tem1$conv
    
    if(SE){
      #print("Commencing estimation of standard error")
      bsd = esti_sd(h,Zs,ds,Xs)
    }else{
      bsd = array(NA, dim = c(nrow(bhat), ncol(bhat)) )
    }
    
    #print("Commencing estimation of constant effect")
    conres  = cons(bhat,h,Zs,ds,Xs)
    if(length(convind) > 0){
      bhat[convind,] = NA
      bsd[convind,]  = NA
    }
    data = cbind(Zs, Xs, ds)
    colnames(data) = c(covname,'obstime','delta')
    colnames(bhat) = covname
    colnames(bsd)  = paste(covname,"_SEE", sep = "")
    
    return(list(temporal_coef = bhat[,(1:r2)+r1], temporal_coef_SEE = bsd[,(1:r2)+r1],
                constant_coef = conres$constant_coef, constant_coef_SEE = conres$constant_coef_SEE,
                bandwidth = h, data = data,
                "Time points (not converged)" = paste('obstime = ', Xs[convind], sep = "")))
  }
}

