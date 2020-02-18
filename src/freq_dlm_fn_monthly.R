## title: functions for mobility forecast
## author: trond
## date: 24.04.2018


## measures of forecast error
rmse <- function(y_true, y_pred) {
    sqrt(mean((y_pred - y_true)^2))
}

mape <- function(y_true, y_pred) {
    mean(abs((y_true - y_pred) / y_true))
}

## function to extend stationarity restrictions
ARtransPars_custom <- function(num_vec, non_zeros = NULL) {
    num_vec[non_zeros] <- ARtransPars(num_vec[non_zeros])
    return(num_vec)
}

## utility functions for model estimation
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
        diag(W(dlm))[c(1, 3)] <- exp(par[1:2])
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
            dlmModARMA(ar=ARtransPars_custom(c(par[5:6], rep(0, 4),
                                               par[7], rep(0, 4), par[8]),
                                             c(1, 2, 7, 12))
                       )
        V(dlm) <- exp(par[1])
        ## diag(W(dlm))[c(2:4, 6, 8:ncol(W(dlm)))] <- 0
        ## diag(W(dlm))[c(1, 5, 7)] <- exp(par[2:4])
        diag(W(dlm))[c(1, 3:4, 6, 8:ncol(W(dlm)))] <- 0
        diag(W(dlm))[c(2, 5, 7)] <- exp(par[2:4])
        ## diag(W(dlm))[c(3:4, 6, 8:ncol(W(dlm)))] <- 0
        ##  diag(W(dlm))[c(1, 2, 5, 7)] <- exp(par[2:5])
    } else if (type == 'lt_arma3') {
        dlm <- dlmModPoly(2, m0 = c(init_level, init_slope), C0 = 2 * diag(2)) +            
            dlmModTrig(s=12, q = 2) +
            dlmModARMA(ar=ARtransPars_custom(c(par[5:6], rep(0, 4),
                                               par[7], rep(0, 4), par[8]),
                                             c(1, 2, 7, 12))
                       )
        ##  dlmModARMA(ar=c(par[5:6], rep(0, 4), par[7], rep(0, 4), par[8]))
        V(dlm) <- exp(par[1])
        diag(W(dlm))[c(2:4, 6, 8:ncol(W(dlm)))] <- 0
        diag(W(dlm))[c(1, 5, 7)] <- exp(par[2:4])
        ## diag(W(dlm))[c(1, 3:4, 6, 8:ncol(W(dlm)))] <- 0
        ## diag(W(dlm))[c(2, 5, 7)] <- exp(par[2:4])
        ## diag(W(dlm))[c(3:4, 6, 8:ncol(W(dlm)))] <- 0
        ## diag(W(dlm))[c(1, 2, 5, 7)] <- exp(par[2:5])
    } else if (type == 'lt_reg_wn') {
        dlm <- dlmModPoly(2, m0 = c(init_level, init_slope), C0 = 2 * diag(2)) +            
            dlmModTrig(s=12, q = 2) +
            dlmModReg(w_m)
        V(dlm) <- exp(par[1])
        diag(W(dlm))[c(2:4, 6)] <- 0
        diag(W(dlm))[c(1, 5, 7:8)] <- exp(par[c(2:5)])
    }
    return(dlm)
}

## plot forecasted values
fore_plot_dt <- function(freq_fore, plot = 'DLM2') {
     ## time period of forecast
    if (grepl('DLM', plot)) {
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
    if (grepl('DLM', plot)) {
        ## forecast level
        fore <- exp(freq_fore$f)
        fore_t <- time(freq_fore$f)
        fsd <- exp(sqrt(unlist(freq_fore$Q)))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    } else {
        fore <- exp(freq_fore$mean)
        fore_t <- time(freq_fore$mean)
        pl <- exp(freq_fore$lower[,2])
        pu <- exp(freq_fore$upper[,2])
    }
    dt <- data.table(t = c(time(y), fore_t),
                     mod = paste(year(tf_year_month), months(tf_year_month)),
                     type = plot,
                     #data = c(y, rep(NA, length(fore))),
                     forecast = c(rep(NA, length(y)), as.numeric(fore)),
                     forecast_pl = c(rep(NA, length(y)), as.numeric(pl)),
                     forecast_pu = c(rep(NA, length(y)), as.numeric(pu))
                     )    
    return(dt)
}

roll_fore_fn <- function(begin_f, end_f, periods, mod, cores, pars=NULL) {
    mclapply(seq(begin_f, end_f),
           function(x) {
                 ## y <- log(window(y_m, end = x)) ## why the f does this not work?
                 y <- ts(log(y_m[1:x]), start = c(1995, 1), frequency = 12)
                 level0 <- y[1]
                 slope0 <- mean(diff(y))
                 if (mod == 'DLM0') {
                     tryCatch({
                         fit <- dlmMLE(y,
                                       par=rep(-4, 3),
                                       type = 'lt_smooth', build=freq_mod,
                                       init_level = level0, init_slope = slope0,
                                       lower = c(-Inf, rep(-Inf,7))
                                       )
                         mdl <- freq_mod(fit$par, type = 'lt_smooth',
                                         init_level = level0, init_slope = slope0)
                         filt <- dlmFilter(y, mdl)
                         return(list(fit$par,
                                     dlmForecast(filt, nAhead = periods)
                                     )
                                )
                     },
                     error = function(e) message(paste0('Error in dlm0:', e))
                     )
                 } else if (mod == 'DLM1') {
                     tryCatch({
                         fit <- dlmMLE(y,
                                       ##par=c(-18, -9, -12, -6, -0.5, -0.4),
                                       par=c(-18, -9, -12, -6, -0.5, -0.4, 0.1, 0.1),
                                       type = 'lt_arma2', build=freq_mod,
                                       init_level = level0, init_slope = slope0,
                                       lower = c(-Inf, rep(-Inf,7)))
                         mdl <- freq_mod(fit$par, type = 'lt_arma2',
                                         init_level = level0, init_slope = slope0)
                         filt <- dlmFilter(y, mdl)
                         return(list(fit$par,
                                     dlmForecast(filt, nAhead = periods)
                                     )
                                )
                     },
                     error = function(e) message(paste0('Error in dlm1:', e))
                     )
                 } else if (mod == 'DLM2') {
                     tryCatch({
                         fit <- dlmMLE(y,
                                       par=c(-18, -9, -12, -6, -0.5, -0.4, 0.1, 0.1),
                                       type = 'lt_arma3', build=freq_mod,
                                       init_level = level0, init_slope = slope0,
                                       lower = c(-Inf, rep(-Inf,7)))
                         mdl <- freq_mod(fit$par, type = 'lt_arma3',
                                         init_level = level0, init_slope = slope0)
                         filt <- dlmFilter(y, mdl)
                         return(list(fit$par,
                                     dlmForecast(filt, nAhead = periods)
                                     )
                                )
                     },
                     error = function(e) message(paste0('Error in dlm1:', e))
                     )
                 } else if (mod == 'ARIMA') {
                     fit <- Arima(y, order=c(2,1,5), seasonal=c(2,1,1))
                     return(forecast(fit, periods))
                 } else if (mod == 'structts') {
                     fit <- StructTS(y, type = 'BSM')
                     return(forecast(fit, periods))
                 } else if (mod == 'TBATS') {
                     fit <- tbats(y, use.parallel = F)
                     return(forecast(fit, periods))
                 } else if (mod == 'ETS') {
                     fit <- ets(y)
                     return(forecast(fit, periods))
                 } else if (mod == 'Naïve') {
                     return(snaive(y, periods))
                 } else if (mod == 'Holt-Winters') {
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
    data.table(f = exp(dat$f),
               pl = exp(dat$f + qnorm(0.05, sd = residuals(dat)$sd)),
               pu = exp(dat$f + qnorm(0.95, sd = residuals(dat)$sd)),
                        ## pl = exp(dat$f) + exp(qnorm(0.25) * residuals(dat)$sd),
                        ## pu = exp(dat$f) - exp(qnorm(0.25) * residuals(dat)$sd),
           t = as.numeric(time(y_m)),
           y = y_m,
           m = mape(window(y_m, start = c(1998, 1)),
                    window(exp(dat$f), start = c(1998, 1))
                    )
           )
}

## calculate errors
err_dt <- function(err_fn) {
    dt1 <- dcast(rf_fit, 'mod ~ type', value.var = err_fn)
    dt2 <- rbindlist(
        list(dt1,
             data.frame(rbind(c(paste0('N=', strsplit(err_fn, split = '_')[[1]][2]),
                                round(apply(dt1[,-1], 2, mean), 4)
                                ),
                              c('',
                                paste0('(', round(apply(dt1[,-1], 2, sd), 4), ')')
                                )
                              )
                        )
             )
    )
    ## [,
    ##   lapply(.SD,
    ##          function(x) as.numeric(as.character(x))),
    ##   by = mod
    ##   ]
    setnames(dt2, 'mod', 'origin')
    return(dt2)
}    

##  Monte Carlo forecasts
mc_run <- function(orig, new) {    
    ## estimate parameters
    dlm2_mle <- dlmMLE(log(window(y_m, end = orig)),
                   par=init_arma3,
                   type = 'lt_arma3',
                   init_level = log(window(y_m, end = orig))[1],
                   init_slope = mean(diff(log(window(y_m, end = orig)))),
                   build=freq_mod,
                   hessian = T
                   )
    ## set up model
    dlm2_fit <- freq_mod(dlm2_mle$par, type = 'lt_arma3',
                     init_level = log(window(y_m, end = orig))[1],
                     init_slope = mean(diff(log(window(y_m, end = orig))))
                     )
    ## kalman filter
    dlm2_filt <- dlmFilter(log(window(y_m, end = orig)),
                           dlm2_fit)
    ## monte carlo forecasts
    dlm2_mc_fore <- dlmForecast(dlm2_filt,
                                nAhead = length(seq(orig + 1/12, 2018 + 11/12, by = 1/12)),
                                sampleNew = new)
    ## gather results, calculate mean predictions and intervals
    out <- rbindlist(lapply(seq_along(dlm2_mc_fore$newObs),
                               function(x) {
                                   fore <- dlm2_mc_fore$newObs[[x]]
                                   data.table(year = c(2016,
                                                       floor(time(window(y_m, start = 2017, end = orig))),
                                                       floor(time(fore))),
                                              c(sum(window(y_m, start = 2016, end = 2016 + 11/12)),
                                                as.numeric(window(y_m, start = 2017, end = orig)),
                                                as.numeric(exp(fore)))
                                              )[,
                                                .(x, sum(V2)),
                                                by = year
                                                ]
                               }))[,
                                   ':=' (                                       
                                       mean_fore = mean(V2[year > 2016]),
                                       median_fore = median(V2[year > 2016]),
                                       lo_fore = quantile(V2[year > 2016], 0.05),
                                       up_fore = quantile(V2[year > 2016], 0.95)),
                                   by = year
                                   ][year == 2016,
                                     4L:7L := V2
                                     ][,
                                       origin := orig
                                       ]
    return(out)
}
