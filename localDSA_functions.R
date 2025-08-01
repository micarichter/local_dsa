## SEIR model ------------------------------------------------------------------
# Kermack-McKendrick SEIR epidemic ordinary differential equations
KMode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # ODEs
    dS <- -beta * S * I
    dE <- beta * S * I - delta * E
    dI <- delta * E - gamma * I
    dR <- gamma * I
    
    # return ODE solution
    list(c(dS, dE, dI, dR))
  })
}

# find approximate solution to the Kermack-McKendrick SEIR ODEs
# diagnostic() and summary() methods for deSolve still work on epidemic
SEIRepidemic <- function(
    beta = 2, delta = 0.5, gamma = 1, rhoE = 0.02, rhoI = 0.05, rhoR = 0, 
    tmin = 0, tmax = 60, tstep = 0.1
) {
  # rho is the prevalence of infection at tmin
  parameters <- c(beta = beta, delta = delta, gamma = gamma, rhoE = rhoE, 
                  rhoI = rhoI, rhoR = rhoR)
  times <- seq(tmin, tmax, by = tstep)
  state <- c(S = 1 - rhoE - rhoI - rhoR, E = rhoE, I = rhoI, R = rhoR)
  KMsolve <- ode(y = state, times = times, func = KMode,
                 parms = parameters)
  
  # collect attributes from deSolve and data.frame into epidemic
  epidemic <- as.data.frame(KMsolve)
  mostattributes(epidemic) <- append(attributes(epidemic), attributes(KMsolve))
  class(epidemic) <- c("data.frame", "deSolve")
  attr(epidemic, "parameters") <- parameters
  attr(epidemic, "tstep") <- tstep
  
  # return epidemic
  epidemic
}

# conveniently retrieve epidemic parameters and R0 on log or original scales
# tstart is the beginning of observation for the human sensor network
SEIRparams <- function(epidemic, tstart, ln = FALSE) {
  params <- attr(epidemic, "parameters")
  R0 <- with(as.list(params), beta / gamma)
  params <- c(params, R0 = R0)
  if (missing(tstart)) {
    tstart <- min(epidemic$time)
  } else {
    rhoEfun <- approxfun(epidemic$time, epidemic$E)
    params["rhoE"] <- rhoEfun(tstart)
    rhoIfun <- approxfun(epidemic$time, epidemic$I)
    params["rhoI"] <- rhoIfun(tstart)
    rhoRfun <- approxfun(epidemic$time, epidemic$R)
    params["rhoR"] <- rhoRfun(tstart)
  }
  if (ln == TRUE) {
    params <- log(params)
    names(params) <- sapply(
      names(params), function(x) paste("ln", x, sep = ""))
  }
  c(params, tstart = tstart)
}

# plot SEIR solution
SEIRplot <- function(
    epidemic, ylim = c(0, 1), curves = "SEIR",
    Scolor = "#23CE6B", Ecolor = "#E3D888", Icolor = "#79ADDC", Rcolor = "#9D5C63", ...
) {
  # epidemic is an object returned by SEIR epidemic
  Stype <- ifelse(grepl("S", curves), "l", "n")
  plot(
    epidemic$time, epidemic$S, type = Stype, col = Scolor, lwd = 2, ylim = ylim,
    xlab = "Time", ylab = "Prevalence", ...)
  Etype <- ifelse(grepl("E", curves), "l", "n")
  lines(epidemic$time, epidemic$E, type = Etype, col = Ecolor, lwd = 2)
  Itype <- ifelse(grepl("I", curves), "l", "n")
  lines(epidemic$time, epidemic$I, type = Itype, col = Icolor, lwd = 2)
  Rtype <- ifelse(grepl("R", curves), "l", "n")
  lines(epidemic$time, epidemic$R, type = Rtype, col = Rcolor, lwd = 2)
  grid()
}


## human sensor network -------------------------------------------------------
# generate human sensor data (observation from tstart to tstop)
HSdat <- function(n, epidemic, tstart, tstop) {
  tmin <- min(epidemic$time)
  tmax <- max(epidemic$time)
  # set tstart and tstop if not specified; check validity if specified
  if (missing(tstart)) tstart <- tmin
  if (missing(tstop)) tstop <- tmax
  # if (!all(c(tstart, tstop) %in% epidemic$time)) {
  #   stop("The start and stop times should both be in epidemic$time.")
  # } else
  if (tstart >= tstop) {
    stop("The observation window (tstart, tstop] must have positive width.")
  } else if (tstart < tmin | tstop > tmax) {
    stop(paste(
      "The observation window (tstart, tstop] must be inside (",
      tmin, ", ", tmax, "]."))
  }
  
  # generate Etimes for complete epidemic
  # Erisk indicates risk of exposure at tstart (?)
  # Estat indicates occurrence of exposure on or before tstart (?)
  # Erisk = 1 and Estat = 0 for right-censored data
  # Erisk = 0 and Estat = 1 for left-censored data
  # Erisk = Estat = 1 for observed exposures
  # Erisk = Estat = 0 should not happen
  cumincE <- 1 - epidemic$S
  invcdf_SE <- approxfun(cumincE, epidemic$time, yleft = tmin, yright = tmax)
  unifsamp <- runif(n)
  Etime_complete <- invcdf_SE(unifsamp)
  Erisk_complete <- (Etime_complete > tmin)
  Estat_complete <- (unifsamp <= max(cumincE))
  
  # E or I at tmin iff Erisk_complete = 0 and Estat_complete = 1
  # Irisk indicates risk of infection at Itime
  # Istat indicates occurrence of infection at Itime
  # Irisk = 1 and Istat = 0 for right-censored data
  # Irisk = 0 and Istat = 1 for left-censored data
  # Irisk = Istat = 1 for observed recoveries
  # Irisk = Istat = 0 for individuals who never get infected
  params <- SEIRparams(epidemic)
  tmin_EIR <- (Erisk_complete == 0) & (Estat_complete == 1)
  rhoE <- params["rhoE"]
  rhoI <- params["rhoI"]
  rhoR <- params["rhoR"]
  tmin_IR <- tmin_EIR & rbinom(n, 1, (rhoI + rhoR) / (rhoE + rhoI + rhoR))
  Irisk_complete <- Estat_complete & !tmin_IR
  
  # generate latent periods, Itimes, and Istat for complete epidemic
  latpd <- ifelse(Irisk_complete, rexp(n, rate = params["delta"]), 0)
  Itime_uncens <- Etime_complete + latpd
  Itime_complete <- pmin(Itime_uncens, tmax)
  Istat_complete <- Estat_complete & (Itime_uncens <= tmax)
  
  # Rrisk indicates risk of recovery at Rtime
  # Rstat indicates occurrence of recovery at Rtime
  # Rrisk = 1 and Rstat = 0 for right-censored data
  # Rrisk = 0 and Rstat = 1 for left-censored data
  # Rrisk = Rstat = 1 for observed recoveries
  # Rrisk = Rstat = 0 for individuals who never get infected
  params <- SEIRparams(epidemic)
  # tmin_IR <- (Irisk_complete == 0) & (Istat_complete == 1)
  tmin_R <- tmin_IR & rbinom(n, 1, rhoR / (rhoI + rhoR))
  Rrisk_complete <- Istat_complete & !tmin_R
  
  # generate infectious periods, Rtimes, and Rstat for complete epidemic
  infpd <- ifelse(Rrisk_complete, rexp(n, rate = params["gamma"]), 0)
  Rtime_uncens <- Itime_complete + infpd
  Rtime_complete <- pmin(Rtime_uncens, tmax)
  Rstat_complete <- Istat_complete & (Rtime_uncens <= tmax)
  
  # restrict epidemic data to (tstart, tstop]
  Erisk <- (Etime_complete > tstart)
  Estat <- Estat_complete & (Etime_complete <= tstop)
  Etime <- ifelse(Estat, ifelse(Erisk, Etime_complete, tstart), tstop)
  
  # The same logic applies to Irisk and Istat. 
  Irisk <- (Itime_complete > tstart) & Estat # test
  Istat <- Istat_complete & (Itime_complete <= tstop)
  Itime <- ifelse(Istat, ifelse(Irisk, Itime_complete, tstart), tstop)
  
  # The same logic applies to Rrisk and Rstat.
  Rrisk <- (Rtime_complete > tstart) & Istat # test
  Rstat <- Rstat_complete & (Rtime_complete <= tstop)
  Rtime <- ifelse(Rstat, ifelse(Rrisk, Rtime_complete, tstart), tstop)
  
  # human sensors data
  hsdata <- data.frame(
    # Etime_complete = Etime_complete, Erisk_complete = Erisk_complete,
    # Estat_complete = Estat_complete,
    # Itime_complete = Itime_complete, Irisk_complete = Irisk_complete,
    # Istat_complete = Istat_complete,
    # Rtime_complete = Rtime_complete, Rrisk_complete = Rrisk_complete,
    # Rstat_complete = Rstat_complete,
    Etime = Etime, Erisk = Erisk, Estat = Estat,
    Itime = Itime, Irisk = Irisk, Istat = Istat,
    Rtime = Rtime, Rrisk = Rrisk, Rstat = Rstat
  )
  attr(hsdata, "tstart") <- tstart
  attr(hsdata, "tstop") <- tstop
  hsdata
}

HSsubset <- function(hsdat, tstart, tstop) {
  with(hsdat, {
    # restrict epidemic data to (tstart, tstop]
    Erisk <- (Etime > tstart)
    Estat <- Estat & (Etime <= tstop)
    Etime <- ifelse(Estat, ifelse(Erisk, Etime, tstart), tstop)
    
    # The same logic applies to Irisk and Istat.
    Irisk <- (Itime > tstart) & Estat # test
    Istat <- Istat & (Itime <= tstop)
    Itime <- ifelse(Istat, ifelse(Irisk, Itime, tstart), tstop)
    
    # The same logic applies to Rrisk and Rstat.
    Rrisk <- (Rtime > tstart) & Istat # test
    Rstat <- Rstat & (Rtime <= tstop)
    Rtime <- ifelse(Rstat, ifelse(Rrisk, Rtime, tstart), tstop)
    
    # human sensors data
    hsdata <- data.frame(
      # Etime_complete = Etime_complete, Erisk_complete = Erisk_complete,
      # Estat_complete = Estat_complete,
      # Itime_complete = Itime_complete, Irisk_complete = Irisk_complete,
      # Istat_complete = Istat_complete,
      # Rtime_complete = Rtime_complete, Rrisk_complete = Rrisk_complete,
      # Rstat_complete = Rstat_complete,
      Etime = Etime, Erisk = Erisk, Estat = Estat,
      Itime = Itime, Irisk = Irisk, Istat = Istat,
      Rtime = Rtime, Rrisk = Rrisk, Rstat = Rstat
    )
    attr(hsdata, "tstart") <- tstart
    attr(hsdata, "tstop") <- tstop
    hsdata
  })
}

# estimate S, E, I, and R curves from human sensors data
HSsurv <- function(
    hsdata, stype = 2, ctype = 2, conf.type = "log-log", level = 0.95, ...
) {
  Esurv <- survfit(
    Surv(Etime, Estat) ~ 1, data = hsdata, stype = stype, ctype = ctype,
    conf.type = conf.type, conf.int = level, ...)
  Isurv <- survfit(
    Surv(Itime, Istat) ~ 1, data = hsdata, stype = stype, ctype = ctype,
    conf.type = conf.type, conf.int = level, ...)
  Rsurv <- survfit(
    Surv(Rtime, Rstat) ~ 1, data = hsdata, stype = stype, ctype = ctype,
    conf.type = conf.type, conf.int = level, ...)
  event_times <- sort(unique(c(Esurv$time, Isurv$time, Rsurv$time)))
  Eprev <- summary(Isurv, event_times)$surv - summary(Esurv, event_times, extend = T)$surv
  Iprev <- summary(Rsurv, event_times)$surv - summary(Isurv, event_times)$surv
  
  # return list of estimates and function arguments
  list(
    Esurv = Esurv, Isurv = Isurv, Rsurv = Rsurv,
    Eprev = data.frame(time = event_times, prev = Eprev),
    Iprev = data.frame(time = event_times, prev = Iprev),
    stype = stype, ctype = ctype, conf.type = conf.type, level = level, ...)
}

# plot HS estimates against SEIR curves
HSplot <- function(
    hsestimate, epidemic, curves = "SEIR", ylim = c(0, 1),
    Scolor = "#23CE6B", Ecolor = "#E3D888", Icolor = "#79ADDC", Rcolor = "#9D5C63", ...
) {
  Esurv <- hsestimate$Esurv
  Isurv <- hsestimate$Isurv
  Rsurv <- hsestimate$Rsurv
  Eprev <- hsestimate$Eprev
  Iprev <- hsestimate$Iprev
  
  # plot S with confidence limits
  Esurv_type <- ifelse(grepl("S", curves), "s", "n")
  Stype <- ifelse(grepl("S", curves), "l", "n")
  plot(
    Esurv$time, Esurv$surv, type = Esurv_type, lty = "dashed",
    col = Scolor, xlab = "Time", ylab = "Prevalence", ylim = ylim, ...)
  lines(
    Esurv$time, Esurv$lower, type = Esurv_type, lty = "dotted", col = Scolor)
  lines(
    Esurv$time, Esurv$upper, type = Esurv_type, lty = "dotted", col = Scolor)
  lines(epidemic$time, epidemic$S, type = Stype, col = Scolor)
  
  # plot E
  Eprev_type <- ifelse(grepl("E", curves), "s", "n")
  Etype <- ifelse(grepl("E", curves), "l", "n")
  lines(
    Eprev$time, Eprev$prev, type = Eprev_type, lty = "dashed", col = Ecolor)
  lines(epidemic$time, epidemic$E, type = Etype, col = Ecolor)
  
  # plot I
  Iprev_type <- ifelse(grepl("I", curves), "s", "n")
  Itype <- ifelse(grepl("I", curves), "l", "n")
  lines(
    Iprev$time, Iprev$prev, type = Iprev_type, lty = "dashed", col = Icolor)
  lines(epidemic$time, epidemic$I, type = Itype, col = Icolor)
  
  # plot R and add grid
  Rsurv_type <- ifelse(grepl("R", curves), "s", "n")
  Rtype <- ifelse(grepl("R", curves), "l", "n")
  lines(
    Rsurv$time, 1 - Rsurv$surv, type = Rsurv_type,
    lty = "dashed", col = Rcolor)
  lines(
    Rsurv$time, 1 - Rsurv$upper, type = Rsurv_type,
    lty = "dotted", col = Rcolor)
  lines(
    Rsurv$time, 1 - Rsurv$lower, type = Rsurv_type,
    lty = "dotted", col = Rcolor)
  lines(epidemic$time, epidemic$R, type = Rtype, col = Rcolor)
  
  # add grid
  grid()
}


## maximum likelihood estimation of SEIR parameters ----------------------------
# DSA log likelihood
nloglikDSA <- function(pvec, data, tstep) {
  # SEIR parameters
  pvec <- as.numeric(pvec)
  lnbeta <- pvec[1]
  lndelta <- pvec[2]
  lngamma <- pvec[3]
  xrhoE <- exp(pvec[4])
  xrhoI <- exp(pvec[5])
  xrhoR <- exp(pvec[6])
  
  # human sensors start and stop times
  tstart <- attr(data, "tstart")
  tstop <- attr(data, "tstop")
  if (missing(tstep)) tstep <- (tstop - tstart) / 200
  
  # solve ODE
  times <- seq(tstart, tstop, tstep)
  rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
  rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
  rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
  state <- c(S = 1 - rhoE - rhoI - rhoR, E = rhoE, I = rhoI, R = rhoR)
  params <- c(beta = exp(lnbeta), delta = exp(lndelta), gamma = exp(lngamma), 
              rhoE, rhoI, rhoR)
  odesolve <- ode(y = state, times = times, func = KMode, parms = params)
  logSfun <- approxfun(
    odesolve[, "time"], log(odesolve[, "S"]),
    yleft = log(1 - rhoE - rhoI - rhoR), yright = log(odesolve[nrow(odesolve), "S"]))
  # logEfun <- approxfun(
  #   odesolve[, "time"], log(odesolve[, "E"]), rule = 2,
  #   yleft = log(rhoE))
  logIfun <- approxfun(
    odesolve[, "time"], log(odesolve[, "I"]), 
    yleft = log(rhoI), yright = log(odesolve[nrow(odesolve), "I"]))
  lnrhoE <- log(rhoE)
  lnrhoI <- log(rhoI)
  lnrhoR <- log(rhoR)
  
  with(data, {
    loglikS <- sum(logSfun(Etime[Erisk == 1])) 
    exposed <- (Erisk == 1) & (Estat == 1)
    loglikE <- (lnbeta * sum(exposed) + sum(logIfun(Etime[exposed])))
    loglikE <- loglikE + lnrhoE * sum((Erisk == 0) & (Irisk == 1))
    infected <- (Irisk == 1) & (Istat == 1)
    loglikI <- lndelta * sum(infected) - exp(lndelta) * sum(Itime - Etime) 
    loglikI <- loglikI + lnrhoI * sum((Irisk == 0) & (Rrisk == 1))
    recovered <- (Rrisk == 1) & (Rstat == 1)
    loglikR <- lngamma * sum(recovered) - exp(lngamma) * sum(Rtime - Itime)
    loglikR <- loglikR + lnrhoR * sum((Rrisk == 0) & (Rstat == 1))
    
    # return negative log likelihood
    -(loglikS + loglikE + loglikI + loglikR)
  })
}

# DSA maximum likelihood estimates
DSAmle <- function(data, init, tstep, level = 0.95, ...) {
  if (missing(init)) {
    init <- c(0, 0, 0, 0, 0, 0)
  }
  names(init) <- c("lnbeta", "lndelta", "lngamma", "lnxrhoE", "lnxrhoI", "lnxrhoR")
  mle <- optim(
    init, nloglikDSA, data = data, tstep = tstep, hessian = TRUE,  ...)
  
  # point estimate, covariance matrix, and log likelihood value
  point <- data.frame(point = mle$par)
  cov <- solve(mle$hessian)
  ll <- -mle$value
  
  # interval estimates
  alpha <- 1 - level
  stderr <- sqrt(diag(cov))
  lower <- point - qnorm(1 - alpha / 2) * stderr
  upper <- point + qnorm(1 - alpha / 2) * stderr
  lower_name <- paste(toString(round(alpha / 2 * 100, 2)), "%", sep = "")
  upper_name <- paste(toString(round((1 - alpha / 2) * 100, 2)), "%", sep = "")
  interval <- data.frame(lower, upper)
  names(interval) <- c(lower_name, upper_name)
  
  # R0 point and interval estimates
  R0coefs <- c(1, 0, -1, 0, 0, 0)
  lnR0 <- R0coefs %*% mle$par
  lnR0_stderr <- sqrt(R0coefs %*% cov %*% R0coefs)
  lnR0_lower <- lnR0 - qnorm(1 - alpha / 2) * lnR0_stderr
  lnR0_upper <- lnR0 + qnorm(1 - alpha / 2) * lnR0_stderr
  lnR0 <- c(point = lnR0, lower = lnR0_lower, upper = lnR0_upper)
  stderr <- c(stderr, lnR0 = lnR0_stderr)
  
  # return list of results
  list(
    point = point, interval = interval, lnR0 = lnR0,
    stderr = stderr, cov = cov, ll = ll)
}

DSAsummary <- function(DSAest, invln = FALSE) {
  tabl <- rbind(cbind(DSAest$point, DSAest$interval), lnR0 = DSAest$lnR0)
  if (invln == TRUE) {
    tabl <- exp(tabl)
    rownames(tabl) <- sapply(
      rownames(tabl), function(x) substring(x, 3, nchar(x)))
  }
  tabl
}

## maximum likelihood model fit and predictions -------------------------------
## model fit and predictions --------------------------------------------------
DSApredict <- function(point, times, samples, level = 0.95) {
  # point estimates are lnbeta, lndelta, lngamma, lnxrhoE, lnxrhoI, and lnxrhoR
  # times is a vector of times at which the SEIR curves should be plotted
  # samples is a set of multivariate normal or MC posterior samples
  # samples and level are ignored unless ci = TRUE
  
  # run ODE at point estimate
  xrhos <- exp(point[4:6])
  rho_total <- 1 + sum(xrhos)
  rhoE <- xrhos[1] / rho_total
  rhoI <- xrhos[2] / rho_total
  rhoR <- xrhos[3] / rho_total
  params <- c(beta = exp(point[1]), delta = exp(point[2]), gamma = exp(point[3]),
              rhoE = rhoE, rhoI = rhoI, rhoR = rhoR)
  
  state <- c(S = 1 - rhoE - rhoI - rhoR, E = rhoE, I = rhoI, R = rhoR)
  KMsolve <- ode(y = state, times = times, func = KMode, parms = params)
  
  # calculate confidence limits if needed and return epidemic
  epidemic <- as.data.frame(KMsolve)
  if (!missing(samples)) {
    bounds <- DSApred_ci(samples, times, level)
    epidemic <- cbind(epidemic, bounds)
    attr(epidemic, "level") <- level
  }
  
  # set attributes and return epidemic
  mostattributes(epidemic) <- append(attributes(epidemic), attributes(KMsolve))
  class(epidemic) <- c("data.frame", "deSolve")
  attr(epidemic, "params") <- params
  attr(epidemic, "times") <- times
  
  # return epidemic
  epidemic
}

DSApred_mlesamp <- function(mle, nsamp = 1000) {
  # sample from MLE multivariate normal approximation
  mu <- mle$point$point
  Sigma <- mle$cov
  mvrnorm(nsamp, mu, Sigma, empirical = TRUE)
}

DSApred_ci <- function(samples, times, level = 0.95) {
  # samples has one row for each sample and four columns for the parms
  # run epidemic for each sample and keep values of S, E, I, and R
  nsamp <- nrow(samples)
  S <- matrix(nrow = nsamp, ncol = length(times))
  colnames(S) <- times
  E <- matrix(nrow = nsamp, ncol = length(times))
  colnames(E) <- times
  I <- matrix(nrow = nsamp, ncol = length(times))
  colnames(I) <- times
  R <- matrix(nrow = nsamp, ncol = length(times))
  colnames(R) <- times
  Rt <- matrix(nrow = nsamp, ncol = length(times))
  colnames(Rt) <- times
  
  for (i in 1:nrow(samples)) {
    parameters <- exp(samples[i, ])
    names(parameters) <- c("beta", "delta", "gamma", "xrhoE", "xrhoI", "xrhoR")
    beta <- parameters["beta"]
    delta <- parameters["delta"]
    gamma <- parameters["gamma"]
    xrhoE <- as.numeric(parameters["xrhoE"])
    xrhoI <- as.numeric(parameters["xrhoI"])
    xrhoR <- as.numeric(parameters["xrhoR"])
    rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
    rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
    rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
    state <- c(S = 1 - rhoE - rhoI - rhoR, E = rhoE, I = rhoI, R = rhoR)
    KMsolve <- ode(
      y = state, times = times, func = KMode, parms = parameters)
    S[i, ] <- KMsolve[, "S"]
    E[i, ] <- KMsolve[, "E"]
    I[i, ] <- KMsolve[, "I"]
    R[i, ] <- KMsolve[, "R"]
    Rt[i, ] <- as.matrix(S[i, ] * beta / gamma)
  }
  
  # find lower and upper quantiles of S, I, and R
  alpha <- 1 - level
  lowerq <- function(x) quantile(x, alpha / 2)
  upperq <- function(x) quantile(x, 1 - alpha / 2)
  median <- function(x) quantile(x, 0.5)
  lowerS <- apply(S, 2, lowerq)
  upperS <- apply(S, 2, upperq)
  lowerE <- apply(E, 2, lowerq)
  upperE <- apply(E, 2, upperq)
  lowerI <- apply(I, 2, lowerq)
  upperI <- apply(I, 2, upperq)
  lowerR <- apply(R, 2, lowerq)
  upperR <- apply(R, 2, upperq)
  lowerRt <- apply(Rt, 2, lowerq)
  upperRt <- apply(Rt, 2, upperq)
  Rt <- apply(Rt, 2, median)
  
  # return confidence limits
  bounds <- data.frame(
    time = times,
    lowerS = lowerS, upperS = upperS,
    lowerE = lowerE, upperE = upperE,
    lowerI = lowerI, upperI = upperI,
    lowerR = lowerR, upperR = upperR,
    lowerRt = lowerRt, upperRt = upperRt,
    Rt = Rt,
    row.names = NULL)
  attr(bounds, "level") <- level
  
  # return bounds
  bounds
}


DSAplot <- function(
    prediction, epidemic, curves = "SEIR", ci, xlim, ylim = c(0, 1),
    Scolor = "#23CE6B", Ecolor = "#E3D888", Icolor = "#79ADDC", Rcolor = "#9D5C63", ...
) {
  # set ci and xlim if missing
  if (missing(ci)) {
    # check whether prediction has a "level" attribute
    ci <- ifelse(is.null(attr(prediction, "level")), FALSE, TRUE)
  }
  if (missing(xlim)) {
    xlim <- c(min(prediction$time), max(prediction$time))
  }
  
  # plot SEIR curves
  SEIRplot(
    epidemic, ylim, curves = curves, xlim = xlim,
    Scolor = Scolor, Ecolor = Ecolor, Icolor = Icolor, Rcolor = Rcolor, ...)
  abline(v = attr(prediction, "tstart"), col = "darkgray")
  abline(v = attr(prediction, "tstop"), col = "darkgray")
  
  # plot DSA estimated S
  Stype <- ifelse(grepl("S", curves), "l", "n")
  lines(
    prediction$time, prediction$S, type = Stype, lty = "dashed", col = Scolor)
  
  # plot DSA estimated E
  Etype <- ifelse(grepl("E", curves), "l", "n")
  lines(
    prediction$time, prediction$E, type = Stype, lty = "dashed", col = Ecolor)
  
  # plot DSA estimated I
  Itype <- ifelse(grepl("I", curves), "l", "n")
  lines(
    prediction$time, prediction$I, type = Itype, lty = "dashed", col = Icolor)
  
  # plot DSA estimated R
  Rtype <- ifelse(grepl("R", curves), "l", "n")
  lines(
    prediction$time, prediction$R, type = Rtype, lty = "dashed", col = Rcolor)
  
  if (ci == TRUE) {
    # plot DSA confidence intervals for S
    lines(
      prediction$time, prediction$lowerS, type = Stype,
      lty = "dotted", col = Scolor)
    lines(
      prediction$time, prediction$upperS, type = Stype,
      lty = "dotted", col = Scolor)
    
    # plot DSA confidence intervals for E
    lines(
      prediction$time, prediction$lowerE, type = Etype,
      lty = "dotted", col = Ecolor)
    lines(
      prediction$time, prediction$upperE, type = Etype,
      lty = "dotted", col = Ecolor)
    
    # plot DSA confidence intervals for I
    lines(
      prediction$time, prediction$lowerI, type = Itype,
      lty = "dotted", col = Icolor)
    lines(
      prediction$time, prediction$upperI, type = Itype,
      lty = "dotted", col = Icolor)
    
    # plot DSA confidence intervals for R
    lines(
      prediction$time, prediction$lowerR, type = Rtype,
      lty = "dotted", col = Rcolor)
    lines(
      prediction$time, prediction$upperR, type = Rtype,
      lty = "dotted", col = Rcolor)
  }
}

# for EpiCast data, no parameters so omit ode comparison
EPIplot <- function(
    hsestimate, curves = "SEIR", ylim = c(0, 1),
    Scolor = "#23CE6B", Ecolor = "#E3D888", Icolor = "#79ADDC", Rcolor = "#9D5C63", ...
) {
  Esurv <- hsestimate$Esurv
  Isurv <- hsestimate$Isurv
  Rsurv <- hsestimate$Rsurv
  Eprev <- hsestimate$Eprev
  Iprev <- hsestimate$Iprev
  
  # plot S with confidence limits
  Esurv_type <- ifelse(grepl("S", curves), "s", "n")
  Stype <- ifelse(grepl("S", curves), "l", "n")
  plot(
    Esurv$time, Esurv$surv, type = Esurv_type, lty = "dashed",
    col = Scolor, xlab = "Time", ylab = "Prevalence", ylim = ylim, ...)
  lines(
    Esurv$time, Esurv$lower, type = Esurv_type, lty = "dotted", col = Scolor)
  lines(
    Esurv$time, Esurv$upper, type = Esurv_type, lty = "dotted", col = Scolor)
  
  # plot E
  Eprev_type <- ifelse(grepl("E", curves), "s", "n")
  Etype <- ifelse(grepl("E", curves), "l", "n")
  lines(
    Eprev$time, Eprev$prev, type = Eprev_type, lty = "dashed", col = Ecolor)
  
  # plot I
  Iprev_type <- ifelse(grepl("I", curves), "s", "n")
  Itype <- ifelse(grepl("I", curves), "l", "n")
  lines(
    Iprev$time, Iprev$prev, type = Iprev_type, lty = "dashed", col = Icolor)
  
  # plot R and add grid
  Rsurv_type <- ifelse(grepl("R", curves), "s", "n")
  Rtype <- ifelse(grepl("R", curves), "l", "n")
  lines(
    Rsurv$time, 1 - Rsurv$surv, type = Rsurv_type,
    lty = "dashed", col = Rcolor)
  lines(
    Rsurv$time, 1 - Rsurv$upper, type = Rsurv_type,
    lty = "dotted", col = Rcolor)
  lines(
    Rsurv$time, 1 - Rsurv$lower, type = Rsurv_type,
    lty = "dotted", col = Rcolor)
  
  # add grid
  grid()
}

# histogram function
stanhist <- function(
    param, true) {
  hist(
    draws$param, prob = TRUE, breaks = "FD", border = "gray", main = NA, 
    xlab = expression(paste("param")),
    ylab = "Density")
  lines(density(draws$param), lty = "dashed")
  grid()
  abline(v = true)
}

# add stat and risk indicators for epicast data
EPIdat <- function(dat, tstart, tstop) {
  tmin <- min(dat$eTime)
  tmax <- max(dat$eTime)
  
  Etime_complete <- dat$eTime
  Erisk_complete <- (Etime_complete > tmin)
  Estat_complete <- dat$estat
  
  Itime_complete <- dat$iTime
  Istat_complete <- dat$istat
  Irisk_complete <- (Estat_complete & Itime_complete > tmin)
  
  Rtime_complete <- dat$rTime
  Rstat_complete <- dat$rstat
  Rrisk_complete <- (Istat_complete & Rtime_complete > tmin)
  
  # restrict epidemic data to (tstart, tstop]
  Erisk <- (Etime_complete > tstart)
  Estat <- Estat_complete & (Etime_complete <= tstop)
  Etime <- ifelse(Estat, ifelse(Erisk, Etime_complete, tstart), tstop)
  
  # The same logic applies to Irisk and Istat. 
  Irisk <- (Itime_complete > tstart) & Estat # test
  Istat <- Istat_complete & (Itime_complete <= tstop)
  Itime <- ifelse(Istat, ifelse(Irisk, Itime_complete, tstart), tstop)
  
  # The same logic applies to Rrisk and Rstat.
  Rrisk <- (Rtime_complete > tstart) & Istat # test
  Rstat <- Rstat_complete & (Rtime_complete <= tstop)
  Rtime <- ifelse(Rstat, ifelse(Rrisk, Rtime_complete, tstart), tstop)
  
  # human sensors data
  hsdata <- data.frame(
    Etime = Etime, Erisk = Erisk, Estat = Estat,
    Itime = Itime, Irisk = Irisk, Istat = Istat,
    Rtime = Rtime, Rrisk = Rrisk, Rstat = Rstat
  )
  attr(hsdata, "tstart") <- tstart
  attr(hsdata, "tstop") <- tstop
  hsdata
}

# fit data and estimate parameters over a sliding window
DSA_window <- function(epidemic) {
  # make parameter plots
  width <- 5
  tstarts <- seq(0, 37, by = 1)
  param_names <- c("lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  
  est <- matrix(nrow = length(tstarts), ncol = 7)
  lwr <- matrix(nrow = length(tstarts), ncol = 7)
  upr <- matrix(nrow = length(tstarts), ncol = 7)
  stderr <- matrix(nrow = length(tstarts), ncol = 7)
  colnames(est) <- param_names
  colnames(lwr) <- param_names
  colnames(upr) <- param_names
  colnames(stderr) <- param_names
  
  for (i in 1:length(tstarts)) {
    try({
      tstart <- tstarts[i]
      hsdat_i <- HSdat(5000, epidemic, tstart = tstart, tstop = tstart + width)
      DSAest_i <- DSAmle(hsdat_i, method = "L-BFGS-B")
      est[i, 1:6] <- DSAest_i$point[, 1]
      est[i, 7] <- DSAest_i$lnR0[1]
      lwr[i, 1:6] <- DSAest_i$interval[, 1] 
      lwr[i, 7] <- DSAest_i$lnR0[2]
      upr[i, 1:6] <- DSAest_i$interval[, 2]
      upr[i, 7] <- DSAest_i$lnR0[3]
      stderr[i, ] <- DSAest_i$stderr
    })
  }
  
  est <- data.frame(est, tstart = tstarts)
  stderr <- data.frame(stderr, tstart = tstarts)
  lwr <- data.frame(lwr, tstart = tstarts)
  upr <- data.frame(upr, tstart = tstarts)
  
  est_long <- est %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "est")
  
  lwr_long <- lwr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "lwr")
  
  upr_long <- upr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "upr")
  
  plot_data <- est_long %>%
    left_join(lwr_long, by = c("tstart", "parameter")) %>%
    left_join(upr_long, by = c("tstart", "parameter"))
  
  p1 <- plot_data %>%
    filter(parameter == "lnbeta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.5), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lnbeta est", 
      title = "log(beta) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p2 <- plot_data %>%
    filter(parameter == "lndelta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.3), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lndelta est", 
      title = "log(delta) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p3 <- plot_data %>%
    filter(parameter == "lngamma") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.2), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lngamma est", 
      title = "log(gamma) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p4 <- plot_data %>%
    filter(parameter == "lnR0") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(2.5), linetype = "solid", color = "black") +
    coord_cartesian(ylim = c(0, 2)) +
    labs(
      x = "tstart", 
      y = "lnR0 est", 
      title = "log(R0) estimates"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  # plot std errors
  s1 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnbeta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(beta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s2 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lndelta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(delta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s3 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lngamma), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(gamma))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s4 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnR0), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(R[0]))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  list(plot_grid(p1, p2, p3, p4) ,plot_grid(s1, s2, s3, s4))
}

# fit epicast data and estimate parameters over a sliding window
EPI_window <- function(dat) {
  width <- 14
  tstarts <- seq(50, 186, by = 1)
  param_names <- c("lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  
  est <- matrix(nrow = length(tstarts), ncol = 7)
  lwr <- matrix(nrow = length(tstarts), ncol = 7)
  upr <- matrix(nrow = length(tstarts), ncol = 7)
  stderr <- matrix(nrow = length(tstarts), ncol = 7)
  colnames(est) <- param_names
  colnames(lwr) <- param_names
  colnames(upr) <- param_names
  colnames(stderr) <- param_names
  
  for (i in 1:length(tstarts)) {
    try({
      tstart <- tstarts[i]
      hsdat_i <- EPIdat(dat[sample(nrow(dat), 5000), ], tstart = tstart, tstop = tstart + width)
      DSAest_i <- DSAmle(hsdat_i, method = "L-BFGS-B")
      est[i, 1:6] <- DSAest_i$point[, 1]
      est[i, 7] <- DSAest_i$lnR0[1]
      lwr[i, 1:6] <- DSAest_i$interval[, 1] 
      lwr[i, 7] <- DSAest_i$lnR0[2]
      upr[i, 1:6] <- DSAest_i$interval[, 2]
      upr[i, 7] <- DSAest_i$lnR0[3]
      stderr[i, ] <- DSAest_i$stderr
    })
  }
  
  est <- data.frame(est, tstart = tstarts)
  stderr <- data.frame(stderr, tstart = tstarts)
  lwr <- data.frame(lwr, tstart = tstarts)
  upr <- data.frame(upr, tstart = tstarts)

  est_long <- est %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "est")
  
  lwr_long <- lwr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "lwr")
  
  upr_long <- upr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "upr")
  
  plot_data <- est_long %>%
    left_join(lwr_long, by = c("tstart", "parameter")) %>%
    left_join(upr_long, by = c("tstart", "parameter"))
  
  p1 <- plot_data %>%
    filter(parameter == "lnbeta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    labs(
      x = "", 
      y = "", 
      title = expression(ln(beta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p2 <- plot_data %>%
    filter(parameter == "lndelta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    labs(
      x = "", 
      y = "", 
      title = expression(ln(delta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p3 <- plot_data %>%
    filter(parameter == "lngamma") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    labs(
      x = "", 
      y = "", 
      title = expression(ln(gamma))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    ) 
  
  p4 <- plot_data %>%
    filter(parameter == "lnR0") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    labs(
      x = "", 
      y = "", 
      title = expression(ln(R[0]))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  e1 <- ggdraw(plot_grid(p1, p2, p3, p4)) +
          draw_label("tstart", x = 0.5, y = 0, vjust = -1, angle = 0, 
                      size = 14) +
                      draw_label("Point estimate", x = 0, 
                      y = 0.5, vjust = 1.5, angle = 90, size = 14)
  
  # plot std errors
  s1 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnbeta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.2)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(beta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s2 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lndelta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.2)) +
    labs(
      x = "", 
      y = "",
      title = expression(ln(delta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s3 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lngamma), color = "blue") +
    coord_cartesian(ylim = c(0, 0.2)) +
    labs(
      x = "", 
      y = "",
      title = expression(ln(gamma))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s4 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnR0), color = "blue") +
    coord_cartesian(ylim = c(0, 0.2)) +
    labs(
      x = "", 
      y = "",
      title = expression(ln(R[0]))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  e2 <- ggdraw(plot_grid(s1, s2, s3, s4)) + 
    draw_label("tstart", x = 0.5, y = 0, vjust = -1, angle = 0, 
               size = 14) +
    draw_label("Standard errors", x = 0, 
               y = 0.5, vjust = 1.5, angle = 90, size = 14)
  
  list(e1, e2)
}

# estimate parameters across a sliding window using synthetic data, and plot
DSA_window <- function(epidemic) {
  # make parameter plots
  width <- 5
  tstarts <- seq(0, 37, by = 1)
  est <- matrix(nrow = length(tstarts), ncol = 7)
  lwr <- matrix(nrow = length(tstarts), ncol = 7)
  upr <- matrix(nrow = length(tstarts), ncol = 7)
  stderr <- matrix(nrow = length(tstarts), ncol = 7)
  colnames(est) <- c(
    "lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  colnames(lwr) <- c(
    "lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  colnames(upr) <- c(
    "lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  colnames(stderr) <- c(
    "lnbeta", "lndelta", "lngamma", "lnrhoE", "lnrhoI", "lnrhoR", "lnR0")
  for (i in 1:length(tstarts)) {
    try({
      tstart <- tstarts[i]
      hsdat_i <- HSdat(5000, epidemic, tstart = tstart, tstop = tstart + width)
      DSAest_i <- DSAmle(hsdat_i, method = "L-BFGS-B")
      est[i, 1:6] <- DSAest_i$point[, 1]
      est[i, 7] <- DSAest_i$lnR0[1]
      lwr[i, 1:6] <- DSAest_i$interval[, 1] 
      lwr[i, 7] <- DSAest_i$lnR0[2]
      upr[i, 1:6] <- DSAest_i$interval[, 2]
      upr[i, 7] <- DSAest_i$lnR0[3]
      stderr[i, ] <- DSAest_i$stderr
    })
  }
  
  est <- data.frame(est)
  est$tstart <- tstarts
  stderr <- data.frame(stderr)
  stderr$tstart <- tstarts
  lwr <- data.frame(lwr)
  lwr$tstart <- tstarts
  upr <- data.frame(upr)
  upr$tstart <- tstarts
  
  est_long <- est %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "est")
  
  lwr_long <- lwr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "lwr")
  
  upr_long <- upr %>%
    pivot_longer(cols = starts_with("ln"), names_to = "parameter", values_to = "upr")
  
  plot_data <- est_long %>%
    left_join(lwr_long, by = c("tstart", "parameter")) %>%
    left_join(upr_long, by = c("tstart", "parameter"))
  
  p1 <- plot_data %>%
    filter(parameter == "lnbeta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.5), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lnbeta est", 
      title = "log(beta) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p2 <- plot_data %>%
    filter(parameter == "lndelta") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.3), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lndelta est", 
      title = "log(delta) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p3 <- plot_data %>%
    filter(parameter == "lngamma") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(0.2), linetype = "solid", color = "black") +
    labs(
      x = "tstart", 
      y = "lngamma est", 
      title = "log(gamma) estimates"
    ) +
    coord_cartesian(ylim = c(-2, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  p4 <- plot_data %>%
    filter(parameter == "lnR0") %>%
    ggplot() +
    geom_ribbon(aes(x = tstart, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(x = tstart, y = est), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = log(2.5), linetype = "solid", color = "black") +
    coord_cartesian(ylim = c(0, 2)) +
    labs(
      x = "tstart", 
      y = "lnR0 est", 
      title = "log(R0) estimates"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  # plot std errors
  s1 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnbeta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(beta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s2 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lndelta), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(delta))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s3 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lngamma), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(gamma))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  s4 <- ggplot(data = stderr) +
    geom_line(aes(x = tstart, y = lnR0), color = "blue") +
    coord_cartesian(ylim = c(0, 0.15)) +
    labs(
      x = "",
      y = "",
      title = expression(ln(R[0]))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  list(plot_grid(p1, p2, p3, p4) ,plot_grid(s1, s2, s3, s4))
}

# end <- 60
# width <- 5
# iter <- 500
# tstarts <- c(5, 15, 25, 35)
# ntstarts <- length(tstarts)
# preds_S <- array(dim = c(6, ntstarts, iter))
# preds_I <- array(dim = c(6, ntstarts, iter))
# bias_S <- array(dim = c(6, ntstarts, iter))
# 
# ptimes <- data.frame(c(15, 25, 30, 40, 50, 60), c(25, 25, 30, 40, 50, 60),
#                      c(35, 35, 35, 40, 50, 60), c(45, 45, 45, 45, 50, 60))
# colnames(ptimes) <- c("five", "fifteen", "twentyfive", "thirtyfive")

# function to make predictions (note: executing this at 500 iterations takes days)
window_predict2 <- function(i) {
  hsdat_i <- HSdat(5000, epidemic, tstart = 0, tstop = 100)
  
  res_S <- array(dim = c(6, ntstarts))
  res_I <- array(dim = c(6, ntstarts))
  lwrS <- array(dim = c(6, ntstarts))
  uprS <- array(dim = c(6, ntstarts))
  lwrI <- array(dim = c(6, ntstarts))
  uprI <- array(dim = c(6, ntstarts))
  
  Sval <- array(dim = c(6, ntstarts))
  Ival <- array(dim = c(6, ntstarts))
  
  empS <- array(dim = c(6, ntstarts))
  empI <- array(dim = c(6, ntstarts))
  bias <- array(dim = c(6, ntstarts))
  
  for (j in 1:ntstarts) {
    tstart <- tstarts[j]
    pred_times <- ptimes[, j]
    
    try({
      # generate human sensors data
      hsdat <- HSsubset(hsdat_i, tstart = tstart, tstop = tstart + width)
      
      # get MLEs
      DSAest <- DSAmle(hsdat, method = "L-BFGS-B")
      
      # make predictions
      mlesamp <- DSApred_mlesamp(DSAest)
      times_pred <- seq(tstart, end, by = 0.1)
      DSApred_mle <- DSApredict(DSAest$point$point, times_pred, mlesamp)
      
      # extract bounds and estimates at each time
      lwrS[, j] <- DSApred_mle$lowerS[match(pred_times, DSApred_mle$time)]
      uprS[, j] <- DSApred_mle$upperS[match(pred_times, DSApred_mle$time)]
      Sval[, j] <- DSApred_mle$S[match(pred_times, DSApred_mle$time)]
      
      lwrI[, j] <- DSApred_mle$lowerI[match(pred_times, DSApred_mle$time)]
      uprI[, j] <- DSApred_mle$upperS[match(pred_times, DSApred_mle$time)]
      Ival[, j] <- DSApred_mle$I[match(pred_times, DSApred_mle$time)]
    })
    
    # get empirical data at time points
    empS[, j] <- epidemic$S[match(pred_times, epidemic$time)]
    empI[, j] <- epidemic$S[match(pred_times, epidemic$time)]
    
    # test if the empirical is within the bounds  
    res_S[, j] <- empS[, j] >= lwrS[, j] & empS[, j] <= uprS[, j]
    res_I[, j] <- empI[, j] >= lwrI[, j] & empI[, j] <= uprI[, j]
    
    # calculate bias
    bias[, j] <- abs(empS[, j] - Sval[, j])
  }
  return(list(res_S = res_S, res_I = res_I, bias = bias))
}

# parallelize computation using parLapply
# n.cores <- detectCores() 
# 
# system.time({
#   clust <- makeCluster(n.cores) 
#   clusterEvalQ(clust, source("localDSA_functions.R"))
#   clusterExport(clust, list("epidemic", "tstarts", "ntstarts", "width", "preds_S", 
#                             "ptimes","preds_I", "end", "iter", "bias_S"))
#   result <- pblapply(1:iter, window_predict2)
#   stopCluster(clust)
# })
#
# put the results in an array
# for (i in 1:iter) {
#   preds_S[,,i] <- result[[i]]$res_S
#   preds_I[,,i] <- result[[i]]$res_I
#   bias_S[,,i] <- result[[i]]$bias
#   
# }
# 
# rownames(preds_S) <- c("time1", "time2", "time3", "time4", "time5", "time6")
# rownames(preds_I) <- c("time1", "time2", "time3", "time4", "time5", "time6")
# 
# colnames(preds_S) <- c("tstart5", "tstart15", "tstart25", "tstart35")
# colnames(preds_I) <- c("tstart5", "tstart15", "tstart25", "tstart35")
# 
# S_true <- apply(preds_S, c(1,2), mean)
# I_true <- apply(preds_I, c(1,2), mean)
# 
# S_true <- as.data.frame(as.table(S_true))
# I_true <- as.data.frame(as.table(I_true))
# colnames(S_true) <- colnames(I_true) <- c("Prediction time", "tstart", 
#                                           "Coverage Probability")
# write.csv(S_true, "s_covs.csv")
# write.csv(I_true, "i_covs.csv")