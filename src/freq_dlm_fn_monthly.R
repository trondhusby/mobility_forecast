## title: functions for mobility forecast
## author: trond
## date: 24.04.2018

rmse <- function(error) {
    sqrt(mean(error^2))
}

mape <- function(y_true, y_pred) {
    mean(abs((y_true - y_pred) / y_true))
}

freq_mod <- function(par, type){
    if (type == 'local_level') {
        dlm <- dlmModPoly(1) + dlmModSeas(12)
        diag(W(dlm))[1:2] <- exp(par[1:2])
        V(dlm) <- exp(par[3])
    } else if (type == 'local_trend') {
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModSeas(12)
        diag(W(dlm))[1:3] <- exp(par[1:3])
        V(dlm) <- exp(par[4])
    } else if (type == 'lt_fourier') {
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModTrig(s = 12)
        diag(W(dlm))[1:3] <- exp(par[1:3])
        V(dlm) <- exp(par[4])
    } else if (type == 'lt_smooth') {
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModSeas(12)
        diag(W(dlm))[2:3] <- exp(par[1:2])
        V(dlm) <- exp(par[3])         
    } else if (type == 'lt_arma') {
       dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModARMA(ar=ARtransPars(par[5:8]), sigma2 = exp(par[9])) +
            dlmModSeas(12)
        diag(W(dlm))[c(1:2, 5)] <- exp(par[2:4])
       V(dlm) <- exp(par[1])
    } else if (type == 'lt_arma_fourier') {
       dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModARMA(ar=ARtransPars(par[5:6]), sigma2 = exp(par[7])) +
            dlmModTrig(s=12)
       diag(W(dlm))[c(1:2, 5)] <- exp(par[2:4])
       V(dlm) <- exp(par[1])
    }  else if (type == 'lt_arma_gibbs') {
       dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModARMA(ar=par[5:6], sigma2 = par[7]) +
            dlmModSeas(12)
        diag(W(dlm))[c(1:2, 5)] <- par[2:4]
        V(dlm) <- par[1]  
    } else if (type == 'lt_arma_2') {
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModSeas(12)  +
            dlmModARMA(ar=par[5:6], sigma2 = exp(par[7]))
        V(dlm) <- exp(par[1])
        diag(W(dlm))[1:3] <- exp(par[2:4])
    }
    return(dlm)
}

smooth_plot_dt <- function(mod) {
    freq_smooth <- dlmSmooth(freq_filt[[mod]])
    if (grepl('arma', mod)) {
        tc <- dropFirst(freq_smooth$s[,1] + freq_smooth$s[,3])
        s <- dropFirst(freq_smooth$s[,5])    
    } else {
        tc <- dropFirst(freq_smooth$s[,1])
        s <- dropFirst(freq_smooth$s[,3])    
    }    
    return(melt(data.table(quarter = time(y_m),
                           mod = mod,
                           data = y_m,
                           trend_cycle = tc,
                           seasonal = s), id.vars = c('mod', 'quarter')
                )
           )
}

fore_plot_dt <- function(mod, plot = 'level') {
    freq_fore <- dlmForecast(freq_filt[[mod]], nAhead = 8)
    if (plot == 'level') {
        ## forecast level
        fore <- freq_fore$f
        fsd <- sqrt(unlist(freq_fore$Q))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    } else if (plot == 'trend'){ ## do something about lt_arma_2
        ## forecast trend
        fore <- freq_fore$a[,1]
        fsd <- sapply(freq_fore$R, function(x) sqrt(x[1,1]))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    }
    if (plot == 'level') {
        dt <- data.table(quarter = c(time(y_m), time(freq_fore$f)),
                     mod = mod,
                     data = c(y_m, rep(NA, length(freq_fore$f))),
                     forecast = c(rep(NA, (length(y_m) -1)), y_m[length(y_m)], fore))
    } else {
        dt <- data.table(quarter = c(time(y_m), time(freq_fore$f)),
                     mod = mod,
                     data = c(smooth_plot_dt(mod)[variable == 'trend_cycle', value], rep(NA, length(freq_fore$f))),
                     forecast = c(rep(NA, length(y_m)), fore))
    }    
    ## reshape data table to long and add confidence interval
    dt_long <- melt(dt, id.vars = c('quarter', 'mod'))[quarter > max(time(y_m)) & variable == 'forecast', forecast_pl := as.numeric(pl)][quarter > max(time(y_m)) & variable == 'forecast', forecast_pu := as.numeric(pu)]
    return(dt_long)
}



gibbs_est <- function(y, a.y, b.y, a.theta, b.theta, shape.y, rate.y,
                     shape.theta, rate.theta,                       
                     dV = 1e-7, m0 = c(level0, slope0, rep(0,13)),
                     C0 = diag(x=c(rep(2,2), rep(1e7,13))),
                     n.sample = 1, thin = 0, save.states = FALSE)
{
    mod <- dlmModPoly(2, dV = dV, dW = rep(1,2)) +
        dlmModARMA(ar = rep(0,2), sigma2 = 1) +
        dlmModSeas(12)
    mod$m0 <- m0
    mod$C0 <- C0
    p <- 15 # dim of state space
    r <- 4 # number of unknown variances
    nobs <- NROW(y)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 )
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop("\"thin\" must be a nonnegative integer")
    ## check hyperpriors for precision(s) of 'theta'
    msg4 <- paste("Either \"a.theta\" and \"b.theta\" or \"shape.theta\"",
                  "and \"rate.theta\" must be specified")
    msg5 <- "Unexpected length of \"shape.theta\" and/or \"rate.theta\""
    msg6 <- "Unexpected length of \"a.theta\" and/or \"b.theta\""
    if (is.null(a.theta))
        if (is.null(shape.theta)) stop(msg4)
        else
            if (is.null(rate.theta)) stop(msg4)
            else
            {
                ## check length of shape.theta and rate.theta
                if (!all(c(length(shape.theta), length(rate.theta)) %in% c(1,r)))
                    warning(msg5)
            }
    else
        if (is.null(b.theta)) stop(msg4)
        else
        {
            if (!all(c(length(a.theta), length(b.theta)) %in% c(1,r)))
                warning(msg6)
            shape.theta <- a.theta^2 / b.theta
            rate.theta <- a.theta / b.theta
        }
    shape.theta <- shape.theta + 0.5 * nobs
    theta <- matrix(0, nobs + 1, p)
    if ( save.states )
        gibbsTheta <- array(0, dim = c(nobs + 1, p, n.sample))
    gibbsV <- vector("numeric", n.sample) 
    gibbsPhi <- matrix(0, nrow = n.sample, ncol = 2)
    gibbsVars <- matrix(0, nrow = n.sample, ncol = r)
    AR2support <- function(u)
    {
        ## stationarity region for AR(2) parameters
        (sum(u) < 1) && (diff(u) < 1) && (abs(u[2]) < 1)
    }
    ARfullCond <- function(u)
    {
        ## log full conditional density for AR(2) parameters
        mod$GG[3:4,3] <- u
        -dlmLL(y, mod) + sum(dnorm(u, sd = c(2,1) * 0.33, log=TRUE))
    }
    it.save <- 0
    for (it in 1:mcmc)
    {
        ## generate AR parameters
        mod$GG[3:4,3] <- arms(mod$GG[3:4,3],
                              ARfullCond, AR2support, 1)
        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)       
        ## generate V
        y.center <- y - tcrossprod(theta[-1, , drop = FALSE], 
            mod$FF)
        SSy <- drop(crossprod(y.center))
        mod$V[] <- 1/rgamma(1, shape = shape.theta,
                       rate = rate.theta + SSy/2
	                )
        ## generate W
        theta.center <- theta[-1,-c(4, 6, 7),drop=FALSE] -
            (theta[-(nobs + 1),,drop=FALSE] %*% t(mod$GG))[,-c(4, 6, 7)]
        SStheta <- drop(sapply( 1 : 4, function(i)
            crossprod(theta.center[,i])))
        diag(mod$W)[c(1:3, 5)] <-
            1 / rgamma(4, shape = shape.theta,
                       rate = rate.theta + 0.5 * SStheta)
        ## save current iteration, if appropriate
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsV[it.save] <- diag(mod$V)
            gibbsPhi[it.save,] <- mod$GG[3:4,3]
            gibbsVars[it.save,] <- diag(mod$W)[c(1:3, 5)]
        }
    }
    if ( save.states )
        return(list(phi = gibbsPhi, dW = gibbsVars, dV = gibbsV, theta = gibbsTheta))
    else
        return(list(phi = gibbsPhi, dW = gibbsVars, dV = gibbsV))
}
