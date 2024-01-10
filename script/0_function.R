pacman::p_load(tidyverse, dplyr, ggplot2, ggthemr, rio, here, magrittr, haven, labelled, sjlabelled,
               sp, raster, rgeos, purrr, tmap, spdep, rgdal, ncdf4)

pacman::p_load(splines, mgcv, nlme, caTools, dlnm, data.table)

extract_labels <- function(dat){
df_labels <- lapply(dat, attr, "label") %>% unlist() %>% as.data.frame()
df_labels <- data.frame(var = row.names(df_labels), label = df_labels$.)

return(df_labels)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

equal_knots <- function(dat, n){
  rangesr <-range(dat,na.rm=T)
  rangest <-range(dat,na.rm=T)

  # BASIS FOR rainfall
  boundr <- rangesr
  varknots <- boundr[1] + diff(boundr)/n*(1:(n-1))
  return(varknots)
}

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-quasipoisson MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

fqbic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  n <- model$y %>% length() # number of observations
  D <- length(coef(model)) # number of model parameters
  qbic <- -2*loglik + D*log(n)
  return(qbic)
}


perform_dlnm <- function(dat_comb,
                         var_knots_rain = quantile(dat_comb$rain)[2:4], # default ns
                         var_knots_temp = quantile(dat_comb$temp)[2:4], # default ns
                         var_spline_rain = "ns",
                         var_spline_temp = "ns",
                         var_degree_rain,
                         var_degree_temp,
                         lag_knots = c(2,4,6,8,10), # assuming constant lag effects in these intervals
                         time_df = 4*5,
                         final = F){# to allow for two peaks

  range_rain <-range(dat_comb$rain,na.rm=T)
  range_temp <-range(dat_comb$temp,na.rm=T)
  lag <- c(0,12) # 0 week to 12 weeks

  # BASIS FOR rainfall
  bound_rain <- range_rain
  #hist(dat_comb$rain)
  if(var_spline_rain %in% c("ns",  "thr")){
    argvar_rain <- list(type= var_spline_rain,
                        knots= var_knots_rain,
                        cen=0
                        )
  }else if(var_spline_rain == "bs"){
    argvar_rain <- list(type = var_spline_rain,
                        knots = var_knots_rain,
                        cen=0,
                        degree = var_degree_rain # cubic spline (=3)
                        )
  }else if(var_spline_rain == "hthr"){
    argvar_rain <- list(type= "thr",
                        side = "h",
                        knots= var_knots_rain,
                        cen=0
    )
  }

  arglag_rain <- list(type="ns",
                      knots = lag_knots,
                      df=6)

  ##########*
  # BASIS FOR TEMPERATURE:
  #hist(dat_comb$temp)
  bound_temp <- range_temp

  if(var_spline_temp %in% c("ns", "thr")){
  argvar_temp <- list(type= var_spline_temp,
                      knots= var_knots_temp,
                      cen=0
                      )
  }else if(var_spline_temp == "bs"){
    argvar_temp <- list(type = var_spline_temp,
                        knots = var_knots_temp,
                        cen=0,
                        degree = var_degree_temp) # cubic spline (3)

  }else if(var_spline_temp == "hthr"){
    argvar_temp <- list(type= "thr",
                        side = "h",
                        knots= var_knots_temp,
                        cen=0
    )
  }

  arglag_temp <- list(type= "ns",
                      knots = lag_knots
                      ) # you don't need to put the knots here. (knots include the complexity)

  suppressWarnings({
    cb_temp <- crossbasis(dat_comb$temp, lag= lag, argvar = argvar_temp, arglag = arglag_temp)
    cb_rain <- crossbasis(dat_comb$rain, lag = lag, argvar = argvar_rain, arglag = arglag_rain)
  })

  model <- glm(n ~ cb_temp +  cb_rain +
                   ns(time, df=time_df), # time as natural spline (df captures the season)
                 family=quasipoisson(),
                 dat_comb)

  print(summary(model))

  # Check if teh RR esimate is reasonable
  cp_rain <- crosspred(cb_rain,
                       model,
                       from=bound_rain[1],
                       to=bound_rain[2], by=.1,
                       cen=median(dat_comb$rain))
  cp_temp <- crosspred(cb_temp,
                       model,
                       from=bound_temp[1],
                       to=bound_temp[2], by=.1,
                       cen=median(dat_comb$temp))

  if(!final){
  return(list(
  model = model,
  aic = fqaic(model),
  bic = fqbic(model),
  valid = if(max(cp_rain$allRRfit)==Inf){"RR estimate is not valid"}
  ))
  }else{
    return(list(
      model = model,
      cb_rain = cb_rain,
      cb_temp = cb_temp,
      aic = fqaic(model),
      bic = fqbic(model)
    ))
  }
}


attrdl <- function(x,basis, name, cases,model=NULL,coef=NULL,vcov=NULL,type="af", cen, range = NULL,
                   dir="back",tot=TRUE, sim=FALSE,nsim=5000) {
  ################################################################################
  #

  .getcoef <- getFromNamespace("getcoef", "dlnm")
  .getvcov <- getFromNamespace("getvcov", "dlnm")
  .getlink <- getFromNamespace("getlink", "dlnm")
  .seqlag <- getFromNamespace("seqlag", "dlnm")
  .mkXpred <- getFromNamespace("mkXpred", "dlnm")

  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen

  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    xlag <- if(dir=="back") tsModel::Lag(x,seq(lag[1],lag[2])) else
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(x)!=diff(lag)+1) stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # cases: TRANFORM IN MEAN OF FUTURE CASES IF dir="forw"
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    cases <- rowMeans(cases)
  } else {
    if(dir=="forw") cases <- rowMeans(tsModel:::Lag(cases,-seq(lag[1],lag[2])))
  }
  #
  ################################################################################
  #
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED
  if(!is.null(model)) {
    cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if(ncol(basis)==1L) cond <- name
    model.class <- class(model)
    coef <- .getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- .getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- .getlink(model,model.class)
    if(!model.link %in% c("log","logit"))
      stop("'model' must have a log or logit link function")
  }
  #  }
  #
  # CHECK IF REDUCED ESTIMATES ARE PROVIDED
  red <- length(coef)!=ncol(basis)
  #
  ################################################################################
  #
  # COMPUTE CROSS-BASIS (OR ONEBASIS IF REDUCED)
  basisnew <- if(!red)
    do.call(crossbasis,list(x=xlag,lag=lag,argvar=attr(basis,"argvar"),
                            arglag=attr(basis,"arglag"))) else
                              do.call(onebasis,c(list(x=x),attr(basis,"argvar")))
  #
  # CHECK DIMENSIONS
  if(length(coef)!=ncol(basisnew))
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if(any(dim(vcov)!=c(length(coef),length(coef))))
    stop("arguments 'coef' and 'vcov' do no match")
  if(red&&dir=="back") stop("only dir='forw' allowed for reduced estimates")
  #
  ################################################################################
  #
  # COMPUTE AF AND AN
  af <- 1-exp(-rowSums(as.matrix(basisnew%*%coef)))
  an <- af*cases
  #
  # TOTAL
  if(tot) {
    isna <- is.na(an)
    an <- sum(an,na.rm=T)
    af <- an/sum(cases[!isna])
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    ansim <- apply(coefsim,2, function(coefi)
      sum((1-exp(-drop(basisnew%*%coefi)))*cases,na.rm=T))
    afsim <- ansim/sum(cases[!isna])
  }
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af
  }
  #
  return(res)
}

#
attrdl <- function(x,basis,cases,model=NULL,coef=NULL,vcov=NULL,model.link=NULL,
                   type="af",dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000) {
  ################################################################################
  #
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- attr(basis,"argvar")$cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    xlag <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2])) else
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(x)!=diff(lag)+1) stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # cases: TRANFORM IN MEAN OF FUTURE CASES IF dir="forw"
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    cases <- rowMeans(cases)
  } else {
    if(dir=="forw") cases <- rowMeans(tsModel:::Lag(cases,-seq(lag[1],lag[2])))
  }
  #
  ################################################################################
  #
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED
  if(!is.null(model)) {
    name <- deparse(substitute(basis))
    cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if(ncol(basis)==1L) cond <- name
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
  }
  #
  # CHECK IF REDUCED ESTIMATES ARE PROVIDED
  red <- length(coef)!=ncol(basis)
  #
  ################################################################################
  #
  # COMPUTE CROSS-BASIS (OR ONEBASIS IF REDUCED)
  basisnew <- if(!red)
    do.call(crossbasis,list(x=xlag,lag=lag,argvar=attr(basis,"argvar"),
                            arglag=attr(basis,"arglag"))) else
                              do.call(onebasis,c(list(x=x),attr(basis,"argvar")))
  #
  # CHECK DIMENSIONS
  if(length(coef)!=ncol(basisnew))
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if(any(dim(vcov)!=c(length(coef),length(coef))))
    stop("arguments 'coef' and 'vcov' do no match")
  if(red&&dir=="back") stop("only dir='forw' allowed for reduced estimates")
  #
  ################################################################################
  #
  # COMPUTE AF AND AN
  af <- 1-exp(-rowSums(as.matrix(basisnew%*%coef)))
  an <- af*cases
  #
  # TOTAL
  if(tot) {
    isna <- is.na(an)
    an <- sum(an,na.rm=T)
    af <- an/sum(cases[!isna])
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    ansim <- apply(coefsim,2, function(coefi)
      sum((1-exp(-drop(basisnew%*%coefi)))*cases,na.rm=T))
    afsim <- ansim/sum(cases[!isna])
  }
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af
  }
  #
  return(res)
}

#



draw_lag_plot <- function(res,
                          dat,
                          cen_rain,
                          cen_temp){
  model <- res$model
  cb_rain <- res$cb_rain
  cb_temp <- res$cb_temp
  bound_rain <- range(dat$rain, na.rm=T)
  bound_temp <- range(dat$temp, na.rm=T)

  # Association between the rainfall and malaria incidence
  cp_rain <- crosspred(cb_rain,
                       model,
                       from=bound_rain[1],
                       to=bound_rain[2], by=.1,
                       cen=cen_rain)
  cr_rain <- crossreduce(cb_rain,
                         model,
                         from=bound_rain[1],
                         to=bound_rain[2], by=.1,
                         cen=cen_temp)

  ## Overall 3D contour plot

  # Lag-specific
  lag_12_rain <- crossreduce(cb_rain,model,type="lag",value=12, from=bound_rain[1],
                             to=bound_rain[2],bylag=0.1,cen=cen_rain)
  lag_8_rain <- crossreduce(cb_rain,model,type="lag",value=8, from=bound_rain[1],
                            to=bound_rain[2],bylag=0.1,cen=cen_rain)
  lag_6_rain <- crossreduce(cb_rain,model,type="lag",value=6, from=bound_rain[1],
                            to=bound_rain[2],bylag=0.1,cen=cen_rain)
  lag_4_rain <- crossreduce(cb_rain,model,type="lag",value=4, from=bound_rain[1],
                            to=bound_rain[2],bylag=0.2,cen=cen_rain)
  lag_2_rain <- crossreduce(cb_rain,model,type="lag",value=2, from=bound_rain[1],
                            to=bound_rain[2],bylag=0.2,cen=cen_rain)
  lag_0_rain <- crossreduce(cb_rain,model,type="lag",value=0, from=bound_rain[1],
                            to=bound_rain[2],bylag=0.2,cen=cen_rain)


  par(mfrow=c(2,3))
  mar=c(0,0,0,0)

  plot(lag_0_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_0_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 0 Weeks",cex=1)


  plot(lag_2_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_2_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 2 Weeks",cex=1)


  plot(lag_4_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_4_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 4 Weeks",cex=1)


  plot(lag_6_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_6_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 6 Weeks",cex=1)


  plot(lag_8_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_8_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 8 Weeks",cex=1)

  plot(lag_12_rain,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(.4,4),xlab=expression("Rainfall (mm)"))
  lines(lag_12_rain,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 12 Weeks",cex=1)


  ## Association between the temperature and the malaria incidence
  cp_temp <- crosspred(cb_temp,
                       model,
                       from=bound_temp[1],
                       to=bound_temp[2], by=.1,
                       cen=cen_temp)

  cr_temp <- crossreduce(cb_temp,
                         model,
                         from=bound_temp[1],
                         to=bound_temp[2], by=.1,
                         cen=cen_temp)

  # Lag-specific
  lag_12_temp <- crossreduce(cb_temp,model,type="lag",value=12, from=bound_temp[1],
                             to=bound_temp[2],bylag=0.1,cen=cen_temp)
  lag_8_temp <- crossreduce(cb_temp,model,type="lag",value=8, from=bound_temp[1],
                            to=bound_temp[2],bylag=0.1,cen=cen_temp)
  lag_6_temp <- crossreduce(cb_temp,model,type="lag",value=6, from=bound_temp[1],
                            to=bound_temp[2],bylag=0.1,cen=cen_temp)
  lag_4_temp <- crossreduce(cb_temp,model,type="lag",value=4, from=bound_temp[1],
                            to=bound_temp[2],bylag=0.2,cen=cen_temp)
  lag_2_temp <- crossreduce(cb_temp,model,type="lag",value=2, from=bound_temp[1],
                            to=bound_temp[2],bylag=0.2,cen=cen_temp)
  lag_0_temp <- crossreduce(cb_temp,model,type="lag",value=0, from=bound_temp[1],
                            to=bound_temp[2],bylag=0.2,cen=cen_temp)


  par(mfrow=c(2,3))
  mar=c(0,0,0,0)

  plot(lag_0_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_0_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 0 Weeks",cex=1)


  plot(lag_2_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_2_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 2 Weeks",cex=1)


  plot(lag_4_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_4_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 4 Weeks",cex=1)


  plot(lag_6_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_6_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 6 Weeks",cex=1)


  plot(lag_8_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_8_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 8 Weeks",cex=1)

  plot(lag_12_temp,type="n",ci="n",
       ylab="Relative Risk",
       ylim=c(0,4.5),xlab=expression("Temperature ("*~degree*C*")"))
  lines(lag_12_temp,col="red",lty=1,lwd=2,ci="area",
        ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
  mtext(text="Lag = 12 Weeks",cex=1)

}

draw_lag_plot2 <- function(model,
                           dat,
                           rain = F,
                           tmin= F,
                           tmax= F,
                           temp = F,
                           cb_rain=NA,
                           cb_tmin=NA,
                           cb_tmax=NA,
                           cb_temp = NA,
                           cen_rain=NA,
                           cen_tmin=NA,
                           cen_tmax=NA,
                           cen_temp = NA){

  if(tmin){
    bound_tmin = range(dat$tmin)
    cp_tmin <- crosspred(cb_tmin,
                         model,
                         from=bound_tmin[1],
                         to=bound_tmin[2], by=.1,
                         cen=cen_tmin)
    cr_tmin <- crossreduce(cb_tmin,
                           model,
                           from=bound_tmin[1],
                           to=bound_tmin[2], by=.1,
                           cen=cen_tmin)

    # Tmin
    lag_12_tmin <- crossreduce(cb_tmin,model,type="lag",value=12, from=bound_tmin[1],
                               to=bound_tmin[2],bylag=0.1,cen=cen_tmin)
    lag_8_tmin <- crossreduce(cb_tmin,model,type="lag",value=8, from=bound_tmin[1],
                              to=bound_tmin[2],bylag=0.1,cen=cen_tmin)
    lag_6_tmin <- crossreduce(cb_tmin,model,type="lag",value=6, from=bound_tmin[1],
                              to=bound_tmin[2],bylag=0.1,cen=cen_tmin)
    lag_4_tmin <- crossreduce(cb_tmin,model,type="lag",value=4, from=bound_tmin[1],
                              to=bound_tmin[2],bylag=0.2,cen=cen_tmin)
    lag_2_tmin <- crossreduce(cb_tmin,model,type="lag",value=2, from=bound_tmin[1],
                              to=bound_tmin[2],bylag=0.2,cen=cen_tmin)
    lag_0_tmin <- crossreduce(cb_tmin,model,type="lag",value=0, from=bound_tmin[1],
                              to=bound_tmin[2],bylag=0.2,cen=cen_tmin)


    par(mfrow=c(2,3))
    mar=c(0,0,0,0)

    plot(lag_0_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_0_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 0 Weeks",cex=1)


    plot(lag_2_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_2_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 2 Weeks",cex=1)


    plot(lag_4_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_4_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 4 Weeks",cex=1)


    plot(lag_6_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_6_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 6 Weeks",cex=1)


    plot(lag_8_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_8_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 8 Weeks",cex=1)

    plot(lag_12_tmin,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Minimum Temperature ("*~degree*C*")"))
    lines(lag_12_tmin,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 12 Weeks",cex=1)
  }
  if(tmax){
    ## Tmax
    bound_tmax = range(dat$tmax)

    cp_tmax <- crosspred(cb_tmax,
                         model,
                         from=bound_tmax[1],
                         to=bound_tmax[2], by=.1,
                         cen=cen_tmax)
    cr_tmax <- crossreduce(cb_tmax,
                           model,
                           from=bound_tmax[1],
                           to=bound_tmax[2], by=.1,
                           cen=cen_tmax)

    lag_12_tmax <- crossreduce(cb_tmax,model,type="lag",value=12, from=bound_tmax[1],
                               to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
    lag_8_tmax <- crossreduce(cb_tmax,model,type="lag",value=8, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
    lag_6_tmax <- crossreduce(cb_tmax,model,type="lag",value=6, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
    lag_4_tmax <- crossreduce(cb_tmax,model,type="lag",value=4, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.2,cen=cen_tmax)
    lag_2_tmax <- crossreduce(cb_tmax,model,type="lag",value=2, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.2,cen=cen_tmax)
    lag_0_tmax <- crossreduce(cb_tmax,model,type="lag",value=0, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.2,cen=cen_tmax)


    par(mfrow=c(2,3))
    mar=c(0,0,0,0)

    plot(lag_0_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_0_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 0 Weeks",cex=1)


    plot(lag_2_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_2_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 2 Weeks",cex=1)


    plot(lag_4_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_4_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 4 Weeks",cex=1)


    plot(lag_6_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_6_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 6 Weeks",cex=1)


    plot(lag_8_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_8_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 8 Weeks",cex=1)

    plot(lag_12_tmax,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
    lines(lag_12_tmax,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 12 Weeks",cex=1)
  }
  if(temp){
    ## temp
    bound_temp = range(dat$temp)

    cp_temp <- crosspred(cb_temp,
                         model,
                         from=bound_temp[1],
                         to=bound_temp[2], by=.1,
                         cen=cen_temp)
    cr_temp <- crossreduce(cb_temp,
                           model,
                           from=bound_temp[1],
                           to=bound_temp[2], by=.1,
                           cen=cen_temp)

    lag_12_temp <- crossreduce(cb_temp,model,type="lag",value=12, from=bound_temp[1],
                               to=bound_temp[2],bylag=0.1,cen=cen_temp)
    lag_8_temp <- crossreduce(cb_temp,model,type="lag",value=8, from=bound_temp[1],
                              to=bound_temp[2],bylag=0.1,cen=cen_temp)
    lag_6_temp <- crossreduce(cb_temp,model,type="lag",value=6, from=bound_temp[1],
                              to=bound_temp[2],bylag=0.1,cen=cen_temp)
    lag_4_temp <- crossreduce(cb_temp,model,type="lag",value=4, from=bound_temp[1],
                              to=bound_temp[2],bylag=0.2,cen=cen_temp)
    lag_2_temp <- crossreduce(cb_temp,model,type="lag",value=2, from=bound_temp[1],
                              to=bound_temp[2],bylag=0.2,cen=cen_temp)
    lag_0_temp <- crossreduce(cb_temp,model,type="lag",value=0, from=bound_temp[1],
                              to=bound_temp[2],bylag=0.2,cen=cen_temp)


    par(mfrow=c(2,3))
    mar=c(0,0,0,0)

    plot(lag_0_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_0_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 0 Weeks",cex=1)


    plot(lag_2_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_2_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 2 Weeks",cex=1)


    plot(lag_4_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_4_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 4 Weeks",cex=1)


    plot(lag_6_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_6_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 6 Weeks",cex=1)


    plot(lag_8_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_8_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 8 Weeks",cex=1)

    plot(lag_12_temp,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,3),xlab=expression("Average Temperature ("*~degree*C*")"))
    lines(lag_12_temp,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 12 Weeks",cex=1)
  }
  if(rain){
    bound_rain = range(dat$rain)
    cp_rain <- crosspred(cb_rain,
                         model,
                         from=bound_rain[1],
                         to=bound_rain[2], by=.1,
                         cen=cen_rain)
    cr_rain <- crossreduce(cb_rain,
                           model,
                           from=bound_rain[1],
                           to=bound_rain[2], by=.1,
                           cen=cen_rain)

    #### Lag-specific plot ###############

    # Rain
    lag_12_rain <- crossreduce(cb_rain,model,type="lag",value=12, from=bound_rain[1],
                               to=bound_rain[2],bylag=0.1,cen=cen_rain)
    lag_8_rain <- crossreduce(cb_rain,model,type="lag",value=8, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.1,cen=cen_rain)
    lag_6_rain <- crossreduce(cb_rain,model,type="lag",value=6, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.1,cen=cen_rain)
    lag_4_rain <- crossreduce(cb_rain,model,type="lag",value=4, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.2,cen=cen_rain)
    lag_2_rain <- crossreduce(cb_rain,model,type="lag",value=2, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.2,cen=cen_rain)
    lag_0_rain <- crossreduce(cb_rain,model,type="lag",value=0, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.2,cen=cen_rain)


    par(mfrow=c(2,3))
    mar=c(0,0,0,0)

    plot(lag_0_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_0_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 0 Weeks",cex=1)


    plot(lag_2_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_2_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 2 Weeks",cex=1)


    plot(lag_4_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_4_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 4 Weeks",cex=1)


    plot(lag_6_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_6_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 6 Weeks",cex=1)


    plot(lag_8_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_8_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 8 Weeks",cex=1)

    plot(lag_12_rain,type="n",ci="n",
         ylab="Relative Risk",
         ylim=c(0,4),xlab=expression("Rainfall (mm)"))
    lines(lag_12_rain,col="red",lty=1,lwd=2,ci="area",
          ci.arg=list(col=rgb(red=0, green=0, blue=1,0.3),col="lightblue"))
    mtext(text="Lag = 12 Weeks",cex=1)
  }


}


