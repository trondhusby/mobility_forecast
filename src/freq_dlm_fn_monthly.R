## title: functions for mobility forecast
## author: trond
## date: 24.04.2018

rmse <- function(y_true, y_pred) {
    sqrt(mean((y_pred - y_true)^2))
}

mape <- function(y_true, y_pred) {
    mean(abs((y_true - y_pred) / y_true))
}

freq_mod <- function(par, type, init_level, init_slope) {
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
            dlmModTrig(s = 12, q = 2)
        diag(W(dlm))[1:3] <- exp(par[1:3])
        ##diag(W(dlm))[which(1:13%%2 != 0)] <- exp(par[1:7])
        V(dlm) <- exp(par[4])
    } else if (type == 'lt_smooth') {
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
            dlmModSeas(12)
        diag(W(dlm))[2:3] <- exp(par[1:2])
        V(dlm) <- exp(par[3])         
    } else if (type == 'lt_arma') {
        level0 <- init_level
        slope0 <- init_slope
        dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2), dW = exp(par[2:3])) +
            dlmModSeas(12, dW = c(exp(par[4]), rep(0, 10))) +
           dlmModARMA(ar=ARtransPars(par[6:7]), sigma = exp(par[5]))
       V(dlm) <- exp(par[1])       
    } else if (type == 'lt_arma_fourier') {
       dlm <- dlmModPoly(2, m0 = c(level0, slope0), C0 = 2 * diag(2)) +
           dlmModTrig(s=12, q=2) +
           dlmModARMA(ar=ARtransPars(par[6:7]))
       diag(W(dlm))[c(1:3, 7)] <- exp(par[2:5])
       V(dlm) <- exp(par[1])       
    } else if (type == 'lt_arma2') {
        dlm <- dlmModPoly(2, m0 = c(init_level, init_slope), C0 = 2 * diag(2)) +            
            dlmModTrig(s=12, q = 2) +            
            dlmModARMA(ar=ARtransPars(par[5:6]))            
        V(dlm) <- exp(par[1])
        diag(W(dlm))[c(2:4, 6, 8:ncol(W(dlm)))] <- 0
        diag(W(dlm))[c(1, 5, 7)] <- exp(par[2:4])
    } else if (type == 'lt_arma3') {
        dlm <- dlmModPoly(2, m0 = c(init_level, init_slope), C0 = 2 * diag(2)) +            
            dlmModTrig(s=12, q = 2) +
            dlmModARMA(ar=c(par[5:6], rep(0, 4), par[7], rep(0, 4), par[8]))
        V(dlm) <- exp(par[1])
        diag(W(dlm))[c(2:4, 6, 8:ncol(W(dlm)))] <- 0
        diag(W(dlm))[c(1, 5, 7)] <- exp(par[2:4])
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

fore_plot_dt <- function(freq_fore, plot = 'dlm') {
     ## time period of forecast
    if (grepl('dlm', plot)) {
        tf <- time(freq_fore$f)[1] - 1/12
    } else {
        tf <- time(freq_fore$mean)[1] - 1/12
    }
    tf_year_month <- as.Date(paste(trunc(tf),
                                   round((tf-floor(tf)) * 12)+1,
                                   '01',
                                   sep = '-'),
                             format = '%Y-%m-%d'
                             )
    y <- window(y_m, end = tf)
    if (grepl('dlm', plot)) {
        ## forecast level
        fore <- exp(freq_fore$f)
        fore_t <- time(freq_fore$f)
        fsd <- exp(sqrt(unlist(freq_fore$Q)))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    } else {     
        fore <- exp(freq_fore$mean)
        fore_t <- time(freq_fore$mean)
        pl <- exp(freq_fore$lower)
        pu <- exp(freq_fore$upper)
    }
    dt <- data.table(t = c(time(y), fore_t),
                     mod = paste(year(tf_year_month), months(tf_year_month)),
                     type = plot,
                     data = c(y, rep(NA, length(fore))),
                     forecast = c(rep(NA, length(y)), fore)
                     )
    ## reshape data table to long and add confidence interval
    dt_long <- melt(dt, id.vars = c('t', 'mod', 'type'))[t > max(time(y)) & variable == 'forecast', forecast_pl := as.numeric(pl)][t > max(time(y)) & variable == 'forecast', forecast_pu := as.numeric(pu)][!is.na(value)]
    return(dt_long)
}


roll_fore_fn <- function(begin, end, periods, mod, cores, pars=NULL) {
    mclapply(time(y_m)[begin:end],
                     function(x) {
                         y <- log(window(y_m, end = x))
                         level0 <- y[1]
                         slope0 <- mean(diff(y))
                         if (mod == 'dlm1') {
                             fit <- dlmMLE(log(y_m),
                                           par=c(-18, -9, -12, -6, -0.5, -0.4),
                                           type = 'lt_arma2', build=freq_mod,
                                           init_level = level0, init_slope = slope0,
                                           lower = c(-Inf, rep(-Inf,7)))
                             mod <- freq_mod(fit$par, type = 'lt_arma2',
                                             init_level = level0, init_slope = slope0)
                             filt <- dlmFilter(y, mod)
                             return(list(fit$par,
                                     dlmForecast(filt, nAhead = periods)
                                     )
                                )                    
                         } else if (mod == 'dlm2') {
                             fit <- dlmMLE(log(y_m),
                                           par=c(-18, -9, -12, -6, -0.5, -0.4, 0.1, 0.1),
                                           type = 'lt_arma3', build=freq_mod,
                                           init_level = level0, init_slope = slope0,
                                           lower = c(-Inf, rep(-Inf,7)))
                             mod <- freq_mod(fit$par, type = 'lt_arma3',
                                             init_level = level0, init_slope = slope0)
                             filt <- dlmFilter(y, mod)
                             return(list(fit$par,
                                     dlmForecast(filt, nAhead = periods)
                                     )
                                    )
                           } else if (mod == 'arima') {
                             fit <- Arima(y, order=c(2,1,5), seasonal=c(2,1,1))
                             return(forecast(fit, periods))
                         } else if (mod == 'tbats') {
                             fit <- tbats(y, use.parallel = F)
                             return(forecast(fit, periods))
                         } else if (mod == 'naive') {                             
                             return(snaive(y, n_periods))
                         } else if (mod == 'hw') {                             
                             return(forecast(HoltWinters(y, seasonal = 'multiplicative'), periods))
                         }
                     },
             mc.cores = cores
             )
}

fore_plot_fn <- function(md='dlm1', err.fn = NULL) {
    dt <- merge(rf_dt,
                 rf_fit,
                 all = T,
                 by = c('mod','type'))[as.numeric(t) > 2013 & type == md]
    ggplot(dt,
           aes(t, value)) +
        geom_line(aes(y = actual_value), col = 'red') +
        geom_line(aes(col = variable)) +
        geom_ribbon(aes(ymin = forecast_pl, ymax = forecast_pu), fill = 'grey60', alpha = 0.2) +
        scale_colour_brewer(palette = 'Set1') +
        facet_wrap(~mod, ncol = 4) +
        theme_bw()
}

filt_plot_dt <- function(dat) {
    data.table(f = dat$f,
                        pl = dat$f + qnorm(0.25) * residuals(dat)$sd,
                        pu = dat$f - qnorm(0.25) * residuals(dat)$sd,
           t = as.numeric(time(y_m)),
           y = log(y_m),
           m = rmse(window(log(y_m), start = c(1998, 1)),
                    window(dat$f, start = c(1998, 1))
                    )
           )
}

err_dt <- function(err_fn) {
    dt1 <- dcast(rf_fit, 'mod ~ type', value.var = err_fn)
    dt2 <- rbindlist(
        list(dt1,
             data.frame(rbind(c(paste('Mean', strsplit(err_fn, split = '_')[[1]][2], sep = ' N='),
                                apply(dt1[,-1], 2, mean)
                                ),
                              c(paste('SD', strsplit(err_fn, split = '_')[[1]][2], sep = ' N='),
                                apply(dt1[,-1], 2, sd)
                                )
                              )
                        )
             )
    )[,
      lapply(.SD,
             function(x) as.numeric(as.character(x))),
      by = mod
      ]
    setnames(dt2, 'mod', 'origin')
    return(dt2)
}    


gibbs_est <- function(y, a.y, b.y, a.theta, b.theta, shape.y, rate.y,
                     shape.theta, rate.theta,                       
                     dV = 1e-7, m0 = c(level0, slope0, rep(0,6)),
                     C0 = diag(x=c(rep(2,2), rep(1e7,6))),
                     n.sample = 1, thin = 0,
                     progressBar = interactive(), save.states = FALSE)
{
    mod <- dlmModPoly(2, dV = dV, dW = rep(0,2)) +
        dlmModTrig(12, q = 2, dV = dV, dW = rep(0,4)) +
        dlmModARMA(ar = rep(-0.9,2), sigma2 = 0)
    mod$m0 <- m0
    mod$C0 <- C0
    p <- 8 # dim of state space
    r <- 7 # number of unknown variances
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
        mod$GG[7:8,7] <- u
        -dlmLL(y, mod) + sum(dnorm(u, sd = c(2,1) * 0.33, log=TRUE))
    }
    if ( progressBar )
        pb <- txtProgressBar(0, mcmc, style = 3)
    it.save <- 0
    for (it in 1:mcmc)
    {
        if ( progressBar )
            setTxtProgressBar(pb, it)
        ## generate AR parameters
        mod$GG[7:8,7] <- arms(mod$GG[7:8,7],
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
        theta.center <- theta[-1,-8,drop=FALSE] -
            (theta[-(nobs + 1),,drop=FALSE] %*% t(mod$GG))[,-8]
        SStheta <- drop(sapply( 1:7, function(i)
            crossprod(theta.center[,i])))
        ##theta.center <- theta[-1,which(mod$FF > 0),drop=FALSE] -
        ##    (theta[-(nobs + 1),,drop=FALSE] %*% t(mod$GG))[,which(mod$FF > 0)]
        ##SStheta <- drop(sapply( 1:4, function(i)
        ##    crossprod(theta.center[,i])))
        diag(mod$W)[-8] <-
            1 / rgamma(7, shape = shape.theta,
                       rate = rate.theta + 0.5 * SStheta)
        ## save current iteration, if appropriate
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsV[it.save] <- diag(mod$V)
            gibbsPhi[it.save,] <- mod$GG[7:8,7]
            gibbsVars[it.save,] <- diag(mod$W)[-8]
        }
    }
    if ( save.states )
        return(list(phi = gibbsPhi, dW = gibbsVars, dV = gibbsV, theta = gibbsTheta))
    else
        return(list(phi = gibbsPhi, dW = gibbsVars, dV = gibbsV))
}



#################################################################
###                                                           ###
###       Unequal variances - t-distributed innovations       ###
###                                                           ###
#################################################################

dlmGibbsDIGt <- function(y, mod, A_y, B_y, A_theta = A_y, B_theta = A_y,
                         nuRange = c(1 : 10, seq(20, 100, by = 10)),
                         alpha_y = 1 / length(nuRange),
                         alpha_theta,
                         n.sample = 1,
                         thin = 0, ind, save.states = FALSE,
                         progressBar = interactive())
#################################################################
#################################################################
###                                                           ###
### y      : data (univariate time series)                    ###
### mod    : model skeleton                                   ###    
### A_y    : upper bound for a_y (prior)                      ###
### B_y    : upper bound for b_y (prior)                      ###
### A_theta: upper bounds for the components of a_theta,      ###
###            recycled if needed (prior)                     ###
### B_theta: upper bounds for the components of b_theta,      ###
###            recycled if needed (prior)                     ###
### nuRange: set of possible values for the degree-of-freedom ###
###            parameter                                      ###
### alpha_y: parameter of the Dirichlet distribution of pi_y  ###
###            (prior)                                        ### 
### alpha_theta: vector parameters of the Dirichlet           ###
###            distributions of the pi_theta's, stored as     ###
###            columns of a matrix; recycled if vector of     ###
###            length 1 or length(nuRange) (prior)            ###
### n.sample: number of MCMC to output                        ###
### thin  : number of MCMC iteration to skip between two      ###        
###           consecutived stored draws                       ###
### ind   : vector of integers specifying the position of the ###
###           unknown, nonconstant variances. If not          ###
###           specified, all are considered unknown           ###
### save.states: logical; should the simulated states be      ###
###           included in the output?                         ###
### progressBar: logical                                      ###
###                                                           ###
#################################################################
#################################################################
{
    m <- NCOL(y)
    nobs <- NROW(y)
    mod$JW <- matrix(0, nrow = ncol(mod$FF), ncol = ncol(mod$FF))
    r <- ncol(mod$FF)
    if ( hasArg(ind) ) {
        ind <- sort(unique(as.integer(ind)))
        s <- 1 : r
        perm <- s[c(ind, s[ !(s %in% ind)])]
        FF(mod) <- mod$FF[, perm, drop = FALSE]
        GG(mod) <- mod$GG[perm, perm, drop = FALSE]
        mod$W <- mod$W[perm, perm, drop = FALSE]
        p <- length(ind)
    }
    else {
        perm <- 1 : r
        p <- r
    }
    diag(mod$JW)[ 1 : p ] <- 1 : p
    mod$JV <- matrix(p + 1)
    mod$X <- matrix(0, nrow = nobs, ncol = p + 1)
    K <- length(nuRange)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 ) 
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop("\"thin\" must be a nonnegative integer")
    if (!all(c(length(A_theta), length(B_theta)) %in% c(1,p)))
        warning("Unexpected length of \"A_theta\" and/or \"B_theta\"")
    A_theta <- rep(A_theta, length.out = p)
    B_theta <- rep(B_theta, length.out = p)
    ablim_theta <- cbind(A_theta, B_theta)
    ablim_y <- c(A_y, B_y)
    if ( !(length(alpha_y) %in% c(1,K)) )
        warning("Unexpected length of \"alpha_y\"")
    alpha_y <- rep(alpha_y, length.out = K)
    if ( hasArg(alpha_theta) ) {
        if ( is.matrix(alpha_theta) ) {
            if ( nrow(alpha_theta) != K || ncol(alpha_theta) != p )
                stop("Wrong dimension of \"alpha_theta\"")
        } else {
            if ( !(length(alpha_theta) %in% c(1,K)) )
                warning("Unexpected length of \"alpha_theta\"")
            alpha_theta <- matrix(rep(alpha_theta, length.out = K*p), nr = K)
        }
    }
    else 
        alpha_theta <- matrix(rep(alpha_y, length.out = K*p), nr = K)
    theta <- matrix(0, nobs + 1, nrow(mod$W))
    ## initialize
    omega_y <- rep(1, nobs)
    nu_y <- rep(100, nobs)
    lambda_y <- 1
    ab_y <- 0.5 * ablim_y
    pi_y <- alpha_y / sum(alpha_y)
    ##
    omega_theta <- matrix(1, nrow = nobs, ncol = p)
    nu_theta <- matrix(100, nrow = nobs, ncol = p)
    lambda_theta <- rep(1, p)
    ab_theta <- 0.5 * ablim_theta
    pi_theta <- alpha_theta / rep(colSums(alpha_theta), each = K) 

    ## memory allocation
    if ( save.states ) 
        gibbsTheta <- array(0, dim = c(dim(theta), n.sample))
    gibbsOmega_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsNu_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsLambda_y <- numeric(n.sample)
    gibbsAB_y <- matrix(0, nrow = n.sample, ncol = 2)
    gibbsPi_y <- matrix(0, nrow = n.sample, ncol = length(nuRange))
    ##
    gibbsOmega_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsNu_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsLambda_theta <- matrix(0, nrow = n.sample, ncol = p)
    gibbsAB_theta <- array(0, dim = c(p, 2, n.sample))
    gibbsPi_theta <- array(0, dim = c(K, p, n.sample))

    ## log target and support for use with `arms'
    ldens.ab <- function(x, lambda, ...) {
        rate <- x[1] / x[2]
        dgamma(lambda, rate * x[1], rate, log = TRUE)
    }
    ind.ab <- function(x, ablim, ...) {
        all( x > 1e-6 & x < ablim )
    }

    ## draw from dirichlet distribution, from package gtools
    rdirichlet <- function (n, alpha) 
    {
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
        sm <- x %*% rep(1, l)
        x/as.vector(sm)
    }

    it.save <- 0
    if ( progressBar )
        pb <- txtProgressBar(0, mcmc, style = 3)
    for (it in 1 : mcmc)
    {
        if ( progressBar )
            setTxtProgressBar(pb, it)
        mod$X[] <- 1 / cbind( omega_theta * rep(lambda_theta, each = nobs),
                             omega_y * lambda_y)

        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)

        ## generate omega_y
        Sq.y.center <- (y - tcrossprod(theta[-1,],mod$FF))^2
        omega_y <- rgamma(nobs, 0.5 * (1 + nu_y),
                                           0.5 * (lambda_y * Sq.y.center + nu_y))
            
        ## generate nu_y
        for (i in 1:nobs) {
            probs <- dgamma(omega_y[i], 0.5 * nuRange, 0.5 * nuRange) * pi_y
            nu_y[i] <- sample(nuRange, size = 1, prob = probs)
        }

        ## generate pi_y
        nuTable <- table(factor(nu_y, levels=nuRange))
        pi_y <- rdirichlet(1, nuTable + alpha_y)
        
        ## generate lambda_y
        SSy <- crossprod(Sq.y.center, omega_y)
        u <- ab_y[1] / ab_y[2]
        lambda_y <- rgamma(1, u * ab_y[1] + 0.5 * nobs,
                           u + 0.5 * SSy)
        
        ## generate a_y and b_y
        ab_y <- arms(ab_y, myldens = ldens.ab, indFunc = ind.ab, n = 1,
                     lambda = lambda_y, ablim = ablim_y)
        
        ## same story for the 'theta' parameters
        ## omega_theta_i
        Sq.theta.center <- (theta[-1,1:p] - tcrossprod(theta[-(nobs + 1), ],
                                                       mod$GG)[,1:p])^2
        omega_theta[] <- rgamma(nobs * p, 0.5 * (1 + nu_theta),
                                0.5 * (rep(lambda_theta, each=nobs) * Sq.theta.center +
                                       nu_theta))
        ## nu_theta_i
        for (j in 1 : nobs)
            for (i in 1 : p)
            {
                probs <- dgamma(omega_theta[j,i], 0.5 * nuRange, 0.5 * nuRange) *
                    pi_theta[, i]
                nu_theta[j, i] <- sample(nuRange, size = 1, prob = probs)
            }

        ## pi_theta_i
        for (i in 1 : p)
        {
            nuTable <- table(factor(nu_theta[, i], levels = nuRange))
            pi_theta[, i] <- rdirichlet(1, nuTable + alpha_theta[, i])
        }

        ## lambda_theta_i
        SStheta <- colSums(Sq.theta.center * omega_theta)
        u <- ab_theta[, 1] / ab_theta[, 2]
        lambda_theta <- rgamma(p, u * ab_theta + 0.5 * nobs,
                               u + 0.5 * SStheta)

        ## a_theta_i & b_theta_i
        for (i in 1 : p)
            ab_theta[i, ] <- arms(ab_theta[i, ], myldens = ldens.ab,
                                 indFunc = ind.ab, n = 1,
                                 lambda = lambda_theta[i], ablim = ablim_theta[i, ])
        
        ## save
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsOmega_y[it.save,] <- omega_y
            gibbsNu_y[it.save,] <- nu_y
            gibbsPi_y[it.save,] <- pi_y
            gibbsLambda_y[it.save] <- lambda_y
            gibbsAB_y[it.save,] <- ab_y
            ##
            gibbsOmega_theta[,,it.save] <- omega_theta
            gibbsNu_theta[,,it.save] <- nu_theta
            gibbsPi_theta[,,it.save] <- pi_theta
            gibbsLambda_theta[it.save,] <- lambda_theta
            gibbsAB_theta[,,it.save] <- ab_theta
       }
    }

    if ( progressBar )
        close(pb)
    if ( save.states )
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta,
                    theta = gibbsTheta[, order(perm), , drop = FALSE]))
    else
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta))

#### end dlmGibbsDIGt ####
}


dlmGibbsDIG_upd <- function (y, mod, a.y, b.y, a.theta, b.theta, shape.y, rate.y, 
    shape.theta, rate.theta, n.sample = 1, thin = 0, ind, save.states = TRUE, 
    progressBar = interactive()) 
{
    msg1 <- "Either \"a.y\" and \"b.y\" or \"shape.y\" and \"rate.y\" must be specified"
    msg2 <- "Unexpected length of \"shape.y\" and/or \"rate.y\""
    msg3 <- "Unexpected length of \"a.y\" and/or \"b.y\""
    msg4 <- paste("Either \"a.theta\" and \"b.theta\" or \"shape.theta\"", 
        "and \"rate.theta\" must be specified")
    msg5 <- "Unexpected length of \"shape.theta\" and/or \"rate.theta\""
    msg6 <- "Unexpected length of \"a.theta\" and/or \"b.theta\""
    msg7 <- "\"thin\" must be a nonnegative integer"
    msg8 <- "multivariate observations are not allowed"
    msg9 <- "inadmissible value of \"ind\""
    if (NCOL(y) > 1) 
        stop(msg8)
    r <- ncol(mod$FF)
    if (hasArg(ind)) {
        ind <- unique(as.integer(ind))
        s <- 1:r
        if (!all(ind %in% s)) 
            stop(msg9)
        perm <- s[c(ind, s[!(s %in% ind)])]
        FF(mod) <- mod$FF[, perm, drop = FALSE]
        GG(mod) <- mod$GG[perm, perm, drop = FALSE]
        W(mod) <- mod$W[perm, perm, drop = FALSE]
        p <- length(ind)
    }
    else {
        perm <- ind <- 1:r
        p <- r
    }
    nobs <- NROW(y)
    if (is.numeric(thin) && (thin <- as.integer(thin)) >= 0) {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else stop(msg7)
    if (!hasArg(a.y)) 
        if (!hasArg(shape.y)) 
            stop(msg1)
        else if (!hasArg(rate.y)) 
            stop(msg1)
        else {
            if (!all(c(length(shape.y), length(rate.y)) == 1)) 
                warning(msg2)
        }
    else if (!hasArg(b.y)) 
        stop(msg1)
    else {
        if (!all(c(length(a.y), length(b.y)) == 1)) 
            warning(msg3)
        shape.y <- a.y^2/b.y
        rate.y <- a.y/b.y
    }
    if (!hasArg(a.theta)) 
        if (!hasArg(shape.theta)) 
            stop(msg4)
        else if (!hasArg(rate.theta)) 
            stop(msg4)
        else {
            if (!all(c(length(shape.theta), length(rate.theta)) %in% 
                c(1, p))) 
                warning(msg5)
        }
    else if (!hasArg(b.theta)) 
        stop(msg4)
    else {
        if (!all(c(length(a.theta), length(b.theta)) %in% c(1, 
            p))) 
            warning(msg6)
        shape.theta <- a.theta^2/b.theta
        rate.theta <- a.theta/b.theta
    }
    shape.y <- shape.y + 0.5 * nobs
    shape.theta <- shape.theta + 0.5 * nobs
    shape.theta <- rep(shape.theta, length.out = p)
    rate.theta <- rep(rate.theta, length.out = p)
    theta <- matrix(0, nobs + 1, r)
    if (save.states) 
        gibbsTheta <- array(0, dim = c(nobs + 1, r, n.sample))
    gibbsV <- vector("numeric", n.sample)
    gibbsW <- matrix(0, nrow = n.sample, ncol = p)
    it.save <- 0
    if (progressBar) 
        pb <- txtProgressBar(0, mcmc, style = 3)
    for (it in 1:mcmc) {
        if (progressBar) 
            setTxtProgressBar(pb, it)
        modFilt <- dlmFilter(y, mod, simplify = TRUE)
        theta[] <- dlmBSample(modFilt)
        y.center <- y - tcrossprod(theta[-1, , drop = FALSE], 
            mod$FF)
        SSy <- drop(crossprod(y.center))
        mod$V[] <- 1/rgamma(1, shape = shape.y, rate = rate.y + 
            0.5 * SSy)
        theta.center <- theta[-1, , drop = FALSE] - tcrossprod(theta[-(nobs + 
            1), , drop = FALSE], mod$GG)
        SStheta <- drop(sapply(1:p, function(i) crossprod(theta.center[, 
            i])))
        SStheta <- colSums((theta[-1, 1:p, drop = FALSE] - tcrossprod(theta[-(nobs + 
            1), , drop = FALSE], mod$GG)[, 1:p])^2)
        diag(mod$W)[1:p] <- 1/rgamma(p, shape = shape.theta, 
            rate = rate.theta + 0.5 * SStheta)
        if (!(it%%every)) {
            it.save <- it.save + 1
            if (save.states) 
                gibbsTheta[, , it.save] <- theta
            gibbsV[it.save] <- diag(mod$V)
            gibbsW[it.save, ] <- diag(mod$W)[1:p]
        }
    }
    colnames(gibbsW) <- paste("W", ind, sep = ".")
    if (progressBar) 
        close(pb)
    if (save.states) 
        return(list(dV = gibbsV, dW = gibbsW, theta = gibbsTheta[, 
            order(perm), , drop = FALSE], mod = mod))
    else return(list(dV = gibbsV, dW = gibbsW))
}
