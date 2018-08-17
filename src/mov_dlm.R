## title: 
## author: trond
## date: 24.04.2018

library(cbsodataR)
library(data.table)
library(ggplot2)
library(dlm)
library(xts)
library(forecast)
library(tseries)
library(gridExtra)
library(mFilter)
library(viridis)
library(forecastHybrid)

## nice links
## https://www.otexts.org/fpp/8/5
## https://sites.ualberta.ca/~sfossati/e509/files/other/dlm_ex.R
## http://hedibert.org/wp-content/uploads/2015/03/EconometriaAvancada-aula7.pdf
## http://docplayer.net/23997208-Rob-j-hyndman-state-space-models-3-arima-and-regarma-models-and-dlm.html
## http://radhakrishna.typepad.com/TimeSeries_Analysis_State_Space_Methods.pdf
## https://www.r-bloggers.com/packages-for-getting-started-with-time-series-analysis-in-r/
## http://blogs2.datall-analyse.nl/2016/02/11/rcode_extended_kalman_filter/
## https://www.uio.no/studier/emner/sv/oekonomi/ECON5101/v14/lecture-6-hg.pdf
## https://lalas.github.io/quantitativeThoughts/r/2014/09/01/dlmTutorial.html
## http://people.duke.edu/~rnau/411home.htm
## http://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/
## http://files.meetup.com/3576292/DLM-PublishedVersion.pdf
## http://r-statistics.co/Time-Series-Forecasting-With-R.html
## http://www.r-programming.org/papers
## http://definetti.uark.edu/dlm/

#cbs_dt[grepl('Voorraad woningen', Title), .(Identifier, Title, ShortTitle, Modified)]

if (!file.exists('../data/cbs_mig_data.RData')) {
    ## get names of tables from CBS api
    cbs_dt <- data.table(cbs_get_table_list(Language="nl"))
    ## find prognose-related tables
    cbs_dt[grepl('bevolkingsontwikkeling', Title), .(Identifier, Title, ShortTitle, Modified)]
    ## dowload data on population (including moves) and gdp
    bev_dt <- data.table(cbs_get_data('37556'))
    bbp_dt <- data.table(cbs_get_data('82262NED'))   
    bev_dt_hf <- data.table(cbs_get_data('37943ned'))
    bbp_dt_hf <- data.table(cbs_get_data('82601NED'))
    woon_dt <- data.table(cbs_get_data('82235NED'))[, lapply(.SD, function(x) as.numeric(as.character(x)))]
    woon_dt_hf <- data.table(cbs_get_data('83906NED'))
    #woon_v_dt_qu <- data.table(cbs_get_data('81955NED'))
    banen_dt <- data.table(cbs_get_data('83599NED'))
    save(bev_dt, bbp_dt, bev_dt_hf, bbp_dt_hf, woon_dt, woon_dt_hf, banen_dt, file = '../data/cbs_mig_data.RData')
} else {
    load('../data/cbs_mig_data.RData')
}

## convert woningvoorraad to quarterly data and apply smoothing
woon_dt[, qu_2 := BeginstandVoorraad_1 + (EindstandVoorraad_8 - BeginstandVoorraad_1)/3]
woon_dt[, qu_3 := BeginstandVoorraad_1 + 2*(EindstandVoorraad_8 - BeginstandVoorraad_1)/3]

woon_dt_long <- melt(woon_dt[,
                             .(Perioden,
                               BeginstandVoorraad_1,
                               qu_2,
                               qu_3,
                               EindstandVoorraad_8)
                             ],
                     id.vars = 'Perioden')

woon_dt_long[,
             variable := gsub('BeginstandVoorraad_1', '1e kwartaal',
                              gsub('qu_2', '2e kwartaal',
                                   gsub('qu_3', '3e kwartaal',
                                        gsub('EindstandVoorraad_8', '4e kwartaal', variable))))
             ][,
               Perioden := paste0(Perioden, ' ', variable)
               ][,
                 variable := NULL
                 ]

## cheeky: add manually first quarter 2018
woon_dt_long <- rbindlist(list(woon_dt_long,
                          data.table(Perioden = '2018 1e kwartaal',
                                     value = 7757737 / 1000 ))
                          )

setkey(woon_dt_long)

## apply smoothing to >1995
woon_dt_long[substr(Perioden, 1, 4) %in% paste0(c(1995:2018)),
             smooth_value := {
                 ind = 1:.N;
                 predict(loess(value ~ ind, span = 0.1))
             }]

## plot results of smoothing
ggplot(woon_dt_long[substr(Perioden, 1, 4) %in% paste0(c(1995:2018))], aes(Perioden, value, group = 1)) +
    geom_line() +
    geom_line(aes(y=smooth_value), col = 'blue') +
    theme_bw()
           

## merge yearly data
freq_dt <- merge(bev_dt[as.numeric(as.character(Perioden)) > 1976, .(Perioden, BinnenGemeentenVerhuisdePersonen_122, TotaalPersonen_123, TotaalBevolking_4)], 
                 bbp_dt[as.numeric(as.character(Perioden)) > 1976, .(Perioden, BrutoBinnenlandsProduct_108)], # gdp 2010 prices
                 by = 'Perioden', all.x = T, all.y = F)[, lapply(.SD, function(x) as.numeric(x)), by = Perioden]
setnames(freq_dt, c('year', 'within_mun_move', 'between_mun_move', 'pop_bp', 'gdp'))

## population is beginning of period, add end of period
freq_dt[, pop_ep := shift(pop_bp, type = 'lead')]

## add mobility numbers for 2017 from quarterly data
freq_dt[year == 2017, ':='(
        within_mun_move = bev_dt_hf[Perioden == 2017, as.numeric(BinnenGemeentenVerhuisdePersonen_10)/1000],
        between_mun_move = bev_dt_hf[Perioden == 2017, as.numeric(TussenGemeentenVerhuisdePersonen_9)/1000],
        pop_ep = bev_dt_hf[Perioden == 2017, as.numeric(BevolkingAanEindVanGekozenPeriode_8)/1000]
        )
        ]

## calculate mobility frequency and transform to long
freq_dt[, pop_ar := 0.5*(pop_bp+pop_ep)][, freq := 1000*(within_mun_move + between_mun_move)/pop_ep]

freq_dt <- melt(freq_dt, id.vars = 'year')[, value := as.numeric(value)][, year := as.numeric(levels(year)[year])]

## plot
ggplot(freq_dt[variable %in% c('freq'), ], aes(year, value, group = 1)) +
    geom_line()+
    facet_wrap(~variable, scales = 'free', ncol = 1) +
    theme_bw()

ggplot(rbind(freq_dt[variable %in% c('freq'), ],
             data.table(year = c(2018:2020), variable = 'freq', value = NA)),
       aes(year, value, group = 1)) +
    geom_line() +
    geom_smooth(method = 'loess') +
    theme_bw() +
    annotate(label = '?', geom = 'text', col = 'red', x= 2019, y = 110, size = 10) +
#    scale_x_discrete(breaks = c(seq(1980 , 2015, by = 5))) +
    ylab('(yearly) moves per 1000 inhabitants')
        
## merge quarterly pop data with gdp data
freq_dt_qu <- merge(bev_dt_hf[,
                              .(Perioden,
                              (as.numeric(TussenGemeentenVerhuisdePersonen_9) +
                               as.numeric(BinnenGemeentenVerhuisdePersonen_10)),
                              as.numeric(BevolkingAanEindVanGekozenPeriode_8)
                              )
                              ],                    
                    bbp_dt_hf[SoortGegevens == 'Prijsniveau 2010',
                           .(Perioden, as.numeric(BrutoBinnenlandsProduct_2))
                           ],
                           by = 'Perioden', all.x = F, all.y = F
                    )[,
                      year := substr(as.character(Perioden), 1, 4)
                      ]


## merge with house data
freq_dt_qu <- merge(freq_dt_qu,
                    woon_dt_hf[,.(Perioden,
                                  PrijsindexBestaandeKoopwoningen_1,
                                  AantalVerkochteWoningen_4,
                                  GemiddeldeVerkoopprijs_7)],
                    by = 'Perioden', all.x = T, all.y = F)[,
                                                           lapply(.SD, function(x) as.numeric(x)),
                                                           by = Perioden]
## attach woonvoorrad
setkey(woon_dt_long, Perioden)
setkey(freq_dt_qu, Perioden)
freq_dt_qu[woon_dt_long, voorraad := i.smooth_value]
                    

setnames(freq_dt_qu, c('quarter', 'moves', 'pop_ep', 'gdp', 'year', 'hp_index', 'n_sales', 'mean_hp', 'stock'))

setkey(freq_dt, year)
setkey(freq_dt_qu, year)
freq_dt_qu[freq_dt[variable == 'pop_ar'], pop_ar_yr := 1000*i.value]

## find end of period population, population at risk and frequenties
freq_dt_qu[nchar(as.character(year)) < 5, pop_bp := shift(pop_ep, type = 'lag')]
freq_dt_qu[nchar(as.character(year)) > 4, pop_bp := shift(pop_ep, type = 'lag')]
freq_dt_qu[, pop_ar := 0.5*(pop_bp + pop_ep)][, freq := 1000*moves/pop_ep]

## frequency of sales
freq_dt_qu[, sales_freq := n_sales / (10 * stock)]

freq_dt_qu <- melt(freq_dt_qu[nchar(as.character(quarter)) > 4, -'year'], id.vars = 'quarter')[, value := as.numeric(value)]

ggplot(freq_dt_qu[variable %in% c('freq', 'mean_hp', 'hp_index', 'sales_freq')], aes(quarter, value, group = 1)) +
    geom_line() +
    facet_wrap(~variable, scales = 'free', ncol = 2) +
    theme_bw() +
    geom_smooth(method = 'loess') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(breaks = freq_dt_qu[grepl('1e kwartaal', as.character(quarter)), unique(quarter)], labels = 1995:2018)


## create time series and check stationarity
y_qu <-  ts(freq_dt_qu[variable == 'freq', value],
            start = c(1995, 1),
            frequency = 4         
         )

x_qu <-  ts(freq_dt_qu[variable == 'hp_index', value],
            start = c(1995, 1),
            frequency = 4         
)

## check whether lags are correlated

## remove seasonality with stl method
test_dt <- data.table(t = time(y_qu),
                      y_s = stl(y_qu, t.window=5, s.window="periodic")$time.series[,1],
                      y_t = stl(y_qu, t.window=5, s.window="periodic")$time.series[,2],
                      x_s = stl(x_qu, t.window=5, s.window="periodic")$time.series[,1],
                      x_t = stl(x_qu, t.window=5, s.window="periodic")$time.series[,2]
                      )

## use hp filter to find cycle
test_dt[,
        c('y_c', 'x_c') := lapply(.SD, function(i) hpfilter(i)$cycle),
        .SDcols = c('y_t', 'x_t')
        ]

## plot
ggplot(melt(test_dt, id.vars = 't')[variable %in% c('y_c', 'x_c')], aes(t, value)) +
    geom_line() +
    facet_wrap(~variable, scales = 'free', ncol = 1) +
    geom_smooth(method = 'loess') +
    theme_bw()

## plot lags
ggplot(melt(test_dt[, .(time(window(y_c, 1995, 2015.5)),
            window(y_c, 1995, 2015.5),
            window(lag(x_c, 10), 1995)
            )], id.vars = 'V1'), aes(V1, value)) +
    geom_line() +
    facet_wrap(~variable, scales = 'free', ncol = 1) +
    geom_smooth(method = 'loess') +
    theme_bw()
   


## check for correlation between lags of x and y
test_cf <- ccf(test_dt[, y_c], test_dt[, x_c])

## create time series and check stationarity
y <-  ts(freq_dt[variable == 'freq', value], 
         c(freq_dt[, min(as.numeric(as.character(year)))])
         ) 

x <-  ts(freq_dt[variable == 'gdp', value], 
         c(freq_dt[, min(as.numeric(as.character(year)))])
)

## dickey fuller test for stationarity       
adf.test(diff(y, 1), alternative = "stationary")

## differentiate: check whether need for yearly diff
diff_fn <- function(par) {
    if (nchar(time(par)[2]) == 4) {
        nd <- ndiffs(par)
        if(nd > 0) {
            par_star <- diff(par, differences=nd)
        }
    } else {
        nd <- ndiffs(diff(par, 4))
        if(nd > 0) {
            par_star <- diff(diff(par,4), differences=nd)
        } else {
            par_star <- diff(par,4)
        }
    }
    print(nd)
    return(par_star)
}

ystar <- diff_fn(y)
xstar <- diff_fn(x)
ystar_qu <- diff_fn(y_qu)
xstar_qu <- diff_fn(x_qu)


## test hodrick prescott filter
test_dt <- data.table(y, time(y))
test_dt[, ':='(
    ystar_hp_trend = hpfilter(y)$trend,
    ystar_hp_cycle = hpfilter(y)$cycle)
    ]

test_dt <- data.table(y_qu, time(y_qu))
test_dt[, ':='(
    ystar_qu_hp_trend = hpfilter(y_qu)$trend,
    ystar_qu_hp_cycle = hpfilter(y_qu)$cycle)
    ]

ggplot(melt(test_dt, id.vars = 'V2'), aes(V2, value)) +
    geom_line() + 
    facet_wrap(~variable, scales = 'free') +
    theme_bw()
    
## plot ACf and PACF
tsdisplay(diff(y_qu, 4))
tsdisplay(diff(diff(y_qu, 4)))

## autofit arima models
y_fit <- auto.arima(window(ystar, 1978, 2016))
y_fit <- auto.arima(y_qu, xreg = x_qu)
y_fit <- auto.arima(diff(diff(y_qu, 4)), xreg = diff(diff(x_qu, 4)))
y_fit <- auto.arima(diff(y_ext, 1))
tsdisplay(residuals(y_fit))

## try different variants
y_fit <- Arima(y_qu, include.mean = F, order=c(0,1,0), seasonal = c(0,0,1))


y_fit <- Arima(window(y, 1977, 2017), include.mean = F, order=c(1,1,0), seasonal = F)

plot(forecast(y_fit, xreg = rep(mean(xstar_qu), 10), h = 10), include = 80)

plot(forecast(y_fit, xreg = rep(mean(x_qu), 40), h=40),include=80)
plot(forecast(y_fit, xreg = (length(y)+1):(length(y)+41), h=40), include = 80)


## (partial) autocorrelation plots
par(mfrow=c(1,2))
Acf(y_qu, main="")  # relationship between y(t) and y(t-k)
Pacf(y_qu, main="") # relationship between y(t) and y(t-k) after removing the effects of other time lags
tsdisplay(residuals(y_fit))
tsdisplay(arima.errors(y_fit), main="ARIMA errors")

Acf(residuals(y_fit))
Box.test(residuals(y_fit), lag=24, fitdf=4, type="Ljung")
plot(forecast(y_fit))

## cross correlation plots
test_cc <- ccf(hpfilter(window(x, 1977, 2016))$cycle, hpfilter(window(y, 1977, 2016))$cycle)

##################################################################
# Forecast of quarterly numbers
##################################################################

qu_mod_mle <- function(par, type = 'seasonal'){
    if (type == 'seasonal') {
        dlm <- dlmModPoly(2) + dlmModSeas(4)
        diag(W(dlm))[2:3] <- exp(par[1:2])
        V(dlm) <- exp(par[3])
    } else if (type == 'trig') {
        dlm <- dlmModPoly(2, dW = c(0, exp(par[1]))) + dlmModTrig(s = 4, dW = exp(par[2]))
        V(dlm) <- exp(par[3])
    } else if (type == 'reg') {
        dlm <- dlmModReg(x_qu) + dlmModSeas(4)
        diag(W(dlm))[1:2] <- exp(par[1:2])
        V(dlm) <- exp(par[3])
    } else if (type == 'arma') {
        dlm <- dlmModPoly(2) + dlmModSeas(4)  +
            dlmModARMA(ar=par[5:6], sigma2 = exp(par[1]))            
        V(dlm) <- exp(par[2])  
        diag(W(dlm))[2:3] <- exp(par[3:4])
    }
    return(dlm)
}

fit <- dlmMLE(y_qu, parm=rep(0, 6), type = 'arma', build=qu_mod_mle, hessian = T)
freq_fit_mle <- qu_mod_mle(fit$par, type = 'arma')
freq_filt <- dlmFilter(y_qu, freq_fit_mle)
#freq_filt <- dlmFilter( c(y_qu, rep(NA, 8)), freq_fit_mle) ## forecasting with rega
freq_fore_plot(freq_fit_mle, plot = 'trend', type = 'arma')


## use mle parameters as prior in bayesian estimation
qu_mod_gibbs <- function(par, type){
    if (type == 'seasonal') {
        dlm <- dlmModPoly(2) + dlmModSeas(4)
        diag(W(dlm))[2:3] <- par[2:3]
        V(dlm) <- par[1]
    } else if (type == 'trig') {
        dlm <- dlmModPoly(2, dW = c(0, par[2])) + dlmModTrig(s = 4, dW = par[3:5])
        V(dlm) <- par[1]
    }
    return(dlm)
}

freq_fit_bayes <- qu_mod_gibbs(gibbs_output[1:5], 'trig')

## gibbs sampling to find parameters
outGibbs <- function(type) {
    if (type == 'seasonal') {
        dlmGibbsDIG(y_qu, dlmModPoly(2) + dlmModSeas(4),
                    a.y = 1, b.y = as.numeric(freq_fit_mle$V), a.theta = 1,
                    b.theta = as.numeric(diag(freq_fit_mle$W))[2:3],
                    n.sample = 1100, ind = c(2:3),
                    save.states = FALSE)
    } else if (type == 'trig') {
        dlmGibbsDIG(y_qu, dlmModPoly(2) + dlmModTrig(4),
                    a.y = 1, b.y = as.numeric(freq_fit_mle$V), a.theta = 1,
                    b.theta = as.numeric(diag(freq_fit_mle$W))[2:5],
                    n.sample = 1100, ind = c(2:5),
                    save.states = FALSE)
    } else if (type == 'arma') {
    }
    return(outGibbs)
}


burn <- 100
gibbs_output <- mcmcMean(cbind(outGibbs$dV[-(1:burn)], outGibbs$dW[-(1:burn), ]))


freq_fit_bayes <- qu_mod_gibbs()
freq_fore_plot(freq_fit_mle, type = 'arma')

freq_filt <- dlmFilter( c(y_qu, rep(NA, 10)), freq_fit_mle)
freq_filt <- dlmFilter(y_qu, freq_fit_mle)


## diagnosis
tsdiag(freq_filt)

## check whether errors are normal distributed
qqnorm(residuals(freq_filt)$res)
qqline(residuals(freq_filt)$res)

## check standard errors of W
aVar <- solve(fit$hessian)
sqrt(diag(aVar))

## check in-sample rmse
rmse <- function(error) {
    sqrt(mean(error^2))
}

rmse(window(freq_filt$f, 1996) - window(y_qu, 1996))

## plot

freq_fore_plot <- function(model, plot = 'level', type) {
    ## filter, smooth and two-year forecast
    freq_filt <- dlmFilter(y_qu,model)
    freq_smooth <- dlmSmooth(freq_filt)
    freq_fore <- dlmForecast(freq_filt, nAhead = 8)
    if (plot == 'level') {
        ## forecast level
        fore <- freq_fore$f
        fsd <- sqrt(unlist(freq_fore$Q))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    } else if (plot == 'trend'){
        ## forecast trend
        fore <- freq_fore$a[,1]
        fsd <- sapply(freq_fore$R, function(x) sqrt(x[1,1]))
        pl <- fore + qnorm(0.05, sd = fsd)
        pu <- fore + qnorm(0.95, sd = fsd)
    }
    ## gather everything in a data table
    if (type == 'arma') {
        freq_mod_dt <- data.table(t = c(time(y_qu), time(freq_fore$f)),
                                  data = c(y_qu, rep(NA, length(freq_fore$f))),
                                  trend = c(dropFirst(freq_smooth$s[,7]), rep(NA, length(freq_fore$f))),
                                  seasonal = c(dropFirst(freq_smooth$s[,3]), rep(NA, length(freq_fore$f))),
                                  filtered = c(freq_filt$f, rep(NA, length(freq_fore$f))),
                                  forecast = c(rep(NA, length(y_qu)), fore)
                                  )
    } else {
        freq_mod_dt <- data.table(t = c(time(y_qu), time(freq_fore$f)),
                                  data = c(y_qu, rep(NA, length(freq_fore$f))),
                                  trend = c(dropFirst(freq_smooth$s[,3]), rep(NA, length(freq_fore$f))),
                                  seasonal = c(dropFirst(freq_smooth$s[,5]), rep(NA, length(freq_fore$f))),
                                  filtered = c(freq_filt$f, rep(NA, length(freq_fore$f))),
                                  forecast = c(rep(NA, length(y_qu)), fore)
                                  )
    }
       ## reshape data table to long
    freq_mod_dt_long <- melt(freq_mod_dt, id.vars = 't')[!is.na(value) & variable == 'forecast', forecast_pl := as.numeric(pl)][!is.na(value) & variable == 'forecast', forecast_pu := as.numeric(pu)]
    ##plot
    return(ggplot(freq_mod_dt_long, aes(t, value)) +
           geom_line(aes(col = variable)) +
           geom_ribbon(data = freq_mod_dt_long[t >2018 & variable == 'forecast'],
                       aes(x = t, ymin = forecast_pl, ymax = forecast_pu), alpha = 0.2) +
           theme_bw())
}


### Gibbs sampler for a Bayesian Analysis of the Output gap

gdpGibbs <- function(y, a.theta, b.theta, shape.theta, rate.theta,
                     dV = 1e-7, m0 = rep(0,4),
                     C0 = diag(x=c(rep(1e7,2), rep(1,2))),
                     n.sample = 1, thin = 0, save.states = FALSE)
{
    mod <- dlmModPoly(2, dV = dV, dW = rep(1,2)) +
        dlmModARMA(ar = rep(0,2), sigma2 = 1)
    mod$m0 <- m0
    mod$C0 <- C0
    p <- 4 # dim of state space
    r <- 3 # number of unknown variances
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
        ## generate W
        theta.center <- theta[-1,-4,drop=FALSE] -
            (theta[-(nobs + 1),,drop=FALSE] %*% t(mod$GG))[,-4]
        SStheta <- drop(sapply( 1 : 3, function(i)
                               crossprod(theta.center[,i])))
        diag(mod$W)[1:3] <-
            1 / rgamma(3, shape = shape.theta,
                       rate = rate.theta + 0.5 * SStheta)
        ## save current iteration, if appropriate
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsPhi[it.save,] <- mod$GG[3:4,3]
            gibbsVars[it.save,] <- diag(mod$W)[1:3]
        }
    }
    if ( save.states )
        return(list(phi = gibbsPhi, vars = gibbsVars, theta = gibbsTheta))
    else
        return(list(phi = gibbsPhi, vars = gibbsVars))
}


## out-of-sample evaluation:
## use 4 different years and 4 different q as end and predict one year ahead

oo_fore <- list()
for (end_yr in seq(2012.25, 2016, by = 0.25)){
    fit <- dlmMLE(window(y_qu, 1995.00, end_yr), parm=rep(0, 6), type = 'arma', build=qu_mod_mle)
    ##freq_fit <- qu_mod_gibbs(gibbs_output[1:3])
    freq_fit <- qu_mod_mle(fit$par, type = 'arma')
    freq_filt <- dlmFilter(window(y_qu, 1995.00, end_yr),freq_fit)
    freq_smooth <- dlmSmooth(freq_filt)
    freq_fore <- dlmForecast(freq_filt, nAhead = 8)$f
    oo_fore[[paste(end_yr)]] <-  data.table(no = end_yr,
                                            t = c(end_yr, as.numeric(time(freq_fore))),
                                            f = c(window(y_qu, end_yr, end_yr), as.numeric(freq_fore)),
                                            rmse = rmse(window(y_qu, end_yr+0.25) - freq_fore)
                                            )
}
oo_fore <- rbindlist(oo_fore)

p1 <- ggplot(data.table(t = as.numeric(time(window(y_qu, 2010))), window(y_qu,
    2010)), aes(t, V2)) + geom_line(col = 'grey') + theme_bw() + geom_line(data
    = oo_fore, aes(t, f), col = 'red') + scale_x_continuous(breaks = seq(1995,
    2018)) + scale_colour_viridis(option = 1, discrete = T) + facet_wrap(~no) +
    ggtitle('Gibbs')

p2 <- ggplot(oo_fore, aes(no, rmse)) +
    geom_point() +
    geom_line() +
    geom_line(aes(no, mean(rmse)), col = 'blue') + 
    theme_bw() +
    ggtitle('Gibbs')

p3 <- ggplot(data.table(t = as.numeric(time(window(y_qu, 2010))), window(y_qu, 2010)), aes(t, V2)) +
    geom_line(col = 'grey') +
    theme_bw() +
    geom_line(data = oo_fore, aes(t, f), col = 'red') +
    scale_x_continuous(breaks = seq(1995, 2018)) +
    scale_colour_viridis(option = 1,  discrete = T) +
    facet_wrap(~no) +
    ggtitle('MLE')

p4 <- ggplot(oo_fore, aes(no, rmse)) +
    geom_point() +
    geom_line() +
    geom_line(aes(no, mean(rmse)), col = 'blue') + 
    theme_bw() +
    ggtitle('MLE')

grid.arrange(p3, p4, ncol = 2)



##################################################################
# Forecast of yearly numbers
##################################################################

## add 2018 to time series
y_ext <-ts(c(y,
             freq_mod_dt[grepl('2018', t),
                         sum(data, na.rm = T) + sum(forecast, na.rm = T)
                         ]
             ),
           start = 1977)

ystar <- diff_fn(y_ext)

## fit an ensemble model
y_fore <- forecast(hybridModel(y_ext, models='aefnt'), 22, PI.combination = 'mean')
accuracy(y_fore)

fore_dt <- data.table(t = seq(1977,2040),
                      data = c(y, rep(NA, 23)),
                      st_forecast = c(rep(NA, 40), y_ext[41:42], rep(NA, 22)),
                      lt_forecast = c(rep(NA, 42), y_fore$mean)
                      )

fore_dt <- melt(fore_dt, id.vars = 't')[!is.na(value) & variable == 'lt_forecast', forecast_pl := y_fore$lower[,2]][!is.na(value) & variable == 'lt_forecast', forecast_pu := y_fore$upper[,2]]

ggplot(fore_dt, aes(t, value)) +
    geom_line(aes(col = variable)) +
    geom_ribbon(data = fore_dt[t >2018 & variable == 'lt_forecast'],
                aes(x = t, ymin = forecast_pl, ymax = forecast_pu), alpha = 0.2) +
    theme_bw()




## if with dlmw


# set parameter restrictions 
parm_rest <- function(parm){
 	return( c(exp(parm[1])/(1+exp(parm[1])),exp(parm[2])) ) 
	}

# set up SS model
arma_fn <- function(parm){
	parm <- parm_rest(parm)
	dlm <- dlmModPoly(2,dV=0.1,dW=c(0.165,0.165)) +
            dlmModARMA(ar=parm[1], ma=NULL, sigma2=parm[2])
        dlm$C0[1,1] <- solve(1-parm[1]^2)*parm[2] 
	return( dlm )
}

## preferred model
arma_fn <- function(p){
    dlm <- dlmModPoly(2, dV=exp(p[1]), dW=exp(p[2:3])) +
        dlmModARMA(ar=p[4], ma = NULL,sigma=exp(p[5]))
    dlm$C0[1,1] <- 10
	return( dlm )
}



# estimate parameters
#fit <- dlmMLE(ystar, parm=c(0,1), build=arma_fn, hessian=T)
fit <- dlmMLE(ystar, parm=c(0,-1,3,0,0), build=arma_fn, hessian=T)
freq_fit <- arma_fn(fit$par)
freq_filt <- dlmFilter(ystar,freq_fit)
freq_smooth <- dlmSmooth(freq_filt)
## two-year forecast
freq_fore <- dlmForecast(freq_filt, nAhead = 23)
freq_fore$f
plot(freq_smooth$s)

tsdiag(freq_filt)

## check whether errors are normal distributed
qqnorm(residuals(freq_filt)$res)
qqline(residuals(freq_filt)$res)

## two-year forecast
freq_fore <- dlmForecast(freq_filt, nAhead = 23)

plot(freq_smooth$s)

## forecast level
fore <- freq_fore$f
fsd <- sqrt(unlist(freq_fore$Q))
pl <- fore + qnorm(0.05, sd = fsd)
pu <- fore + qnorm(0.95, sd = fsd)


## forecast trend
#fore <- freq_fore$a[,1]
#fsd <- sapply(freq_fore$R, function(x) sqrt(x[1,1]))
#pl <- fore + qnorm(0.05, sd = fsd)
#pu <- fore + qnorm(0.95, sd = fsd)

freq_mod_dt <- data.table(t = c(time(ystar), time(freq_fore$f)),
                          data = c(ystar, rep(NA, length(freq_fore$f))),
                          cycle = c(freq_smooth$s[,1], rep(NA, length(freq_fore$f))),
                          trend = c(freq_smooth$s[,3], rep(NA, length(freq_fore$f))),
                          filtered = c(freq_filt$f, rep(NA, length(freq_fore$f))),
                          forecast = c(rep(NA, length(ystar)), fore)
                          )


freq_mod_dt_long <- melt(freq_mod_dt, id.vars = 't')[!is.na(value) & variable == 'forecast', forecast_pl := as.numeric(pl)][!is.na(value) & variable == 'forecast', forecast_pu := as.numeric(pu)]

ggplot(freq_mod_dt_long, aes(t, value)) +
    geom_line(aes(col = variable)) +
    geom_ribbon(data = freq_mod_dt_long[t >2018 & variable == 'forecast'],
                aes(x = t, ymin = forecast_pl, ymax = forecast_pu), alpha = 0.2) +
    theme_bw()

ts(diffinv(as.numeric(fore), xi = y_ext[42]), 2018)



y_fit <- Arima(y_ext, include.mean = F, order=c(0,1,0))
forecast(y_fit, h=24)


########################### old


## with polynomial trend
grid.arrange(
    ggplot(data.table(year = 1978:2016, ystar = window(ystar, 1978, 2016), freq_smooth$s[-1,c(1,3)], c(freq_filt$f[-1], NA)),
           aes(x = year)) +
    geom_line(aes(y = ystar)) + # data
    geom_line(aes(y = V1), col = 'blue') + # cycle
    geom_line(aes(y = V2), col = 'green') + # smooth
    geom_line(aes(y = V4), col = 'red') + # filtered
    theme_bw(),
    ggplot(test_dt, aes(V3, ystar)) +
    geom_line() +
    geom_line(aes(y = ystar_hp_trend), col = 'blue') +
    geom_line(aes(y = ystar_hp_cycle), col = 'green') +
theme_bw()
)

forecast_fit <- dlmForecast(freq_filt, nAhead=5,method='plain')
fsd <- sqrt(unlist(forecast_fit$Q))
pl <- forecast_fit$f + qnorm(0.05, sd = fsd)
pu <- forecast_fit$f + qnorm(0.95, sd = fsd)

ggplot(melt(data.table(seq(2017,(2016+length(forecast_fit$f))),
                       as.numeric(forecast_fit$f),
                       as.numeric(pl),
                       as.numeric(pu)),
            id.vars = 'V1'),
       aes(V1, value)) +
    geom_line(aes(col = variable)) +
    theme_bw()





fn <- function(params){
  model <- dlmModARMA(ar = params[1],
                      ma = NULL,
                      sigma2 = exp(params[2]),
                      dV = 1e-7)
}
fit <- dlmMLE(ystar,rep(0,2),fn)
#ARMA parameters
print(c(fit$par[1:2]))

# set parameter restrictions 
parm_rest <- function(parm){
  return( c(exp(parm[1])/(1+exp(parm[1])),exp(parm[2])) ) 
}

ssm9 <- function(parm){
  parm <- parm_rest(parm)
  dlm <- dlmModARMA(ar=parm[1], ma=NULL, sigma2=parm[2])
  dlm$C0 <- solve(1-parm[1]^2)*parm[2] 
  return(dlm)
}

ssm9 <- function(parm){
  parm <- ARtransPars(parm)
  dlm <- dlmModARMA(ar=parm[1:5], ma=NULL, sigma2=exp(parm[6]))
#  dlm$C0 <- solve(1-parm[1]^2)*parm[2] 
  return(dlm)
}

ARtransPars(c(5,5))

solve(1-0^2)*1

# Case 1: fit model using 230 observations and 
# forecast next 20 observations using dlmForecast

# estimate parameters
fit9 <- dlmMLE(y=ystar,parm=c(rep(0.5,5), 1),build=ssm9,hessian=T)

# filter and smooth
mod9 <- ssm9(fit9$par)
mod9f <- dlmFilter(ystar,mod9)
mod9s <- dlmSmooth(mod9f)

# forecast using SS model
fore9 <- dlmForecast(mod9f,nAhead=24,method='plain')

## forecast first differences
fore_dt <- data.table(c(time(ystar), time(fore9$f)), 
                      c(ystar, fore9$f),
                      c(rep(NA, 38), ystar[39], fore9$f + qnorm(.95) * sqrt(unlist(fore9$Q))),
                      c(rep(NA, 38), ystar[39], fore9$f - qnorm(.95) * sqrt(unlist(fore9$Q))) 
)
setnames(fore_dt, c('year', 'point_estimate', 'lower_95_CI',  'upper_95_CI'))

ggplot(fore_dt, aes(year, point_estimate)) +
  geom_line()+
  geom_ribbon(aes(ymin=lower_95_CI, ymax=upper_95_CI), alpha = 0.2 )

## recover time series as inverse of 1st diff
freq_fc <- data.table(year = 1977:2040,
                     c(y[-40], diffinv(ts(as.vector(fore9$f), start = 2017), xi=y[40]))
)

ggplot(freq_fc, aes(year, V2)) +
  geom_line() 


# filtered values
xtfilt <- ts(mod9f$m[-1],start=1)
# width of 90% confidence interval
sefilt <- dlmSvd2var(mod9f$U.C,mod9f$D.C)[-1]
width <- qnorm(.95) * sqrt(sapply(sefilt,diag))
# put filtered results together
xtfiltered <- cbind(xtfilt,as.vector(xtfilt)+width%o%c(1,-1))

# smoothed values
xtsmooth <- ts(mod9s$s[-1],start=1)
# width of 90% confidence interval
sesmooth <- dlmSvd2var(mod9s$U.S,mod9s$D.S)[-1]
width <- qnorm(.95) * sqrt(sapply(sesmooth,diag))
# put smoothed results together
xtsmoothed <- cbind(xtsmooth,as.vector(xtsmooth)+width%o%c(1,-1))


## smoothed 
LSARsm <- dlmSmooth(ystar,LSARMod)

smoothed <- LSARsm$s[-1, 1] + xstar * LSARsm$s[-1, 2]

ggplot(data.table(time(ystar), ystar, smoothed), aes(V1, ystar)) +
  geom_line() +
  geom_line(aes(y=smoothed), col = 'red')

#################################################################################
## old
#################################################################################

## retrieve Prognose bevolking; geslacht en leeftijd, 2017-2060
if(!file.exists('data/mov_dt.RData')) {
    bev_progn_dt <- data.table(cbs_get_data('83597NED')) # population (forecast)
    bev_dt <- data.table(cbs_get_data('7461bev')) # population (historic)
    mv_dt <- data.table(cbs_get_data('60048ned')) # moves
    bbp_dt <- data.table(cbs_get_data('82601NED')) # gdp
    save(bev_progn_dt, bev_dt, mv_dt, bbp_dt, file = 'data/mov_dt.RData')
} else {
    load(file= 'data/mov_dt.RData')
}

## calculate frequency and add gdp
setkey(mv_dt, Perioden)
setkey(bev_dt, Perioden)
 
mv_dt[bev_dt[Geslacht == 'Totaal mannen en vrouwen' & Leeftijd == 'Totaal', ], Pop_total := i.TotaleBevolking_1 ]
mv_dt[bev_dt[Geslacht == 'Mannen' & Leeftijd == 'Totaal', ], Pop_mannen := i.TotaleBevolking_1 ]
mv_dt[bev_dt[Geslacht == 'Vrouwen' & Leeftijd == 'Totaal', ], Pop_vrouwen := i.TotaleBevolking_1 ]
mv_dt[RegioS == 'Nederland', Verhuismobiliteit_mannen := as.numeric(TotaalMannen_61) + as.numeric(TotaalMannen_101)]
mv_dt[RegioS == 'Nederland', Verhuismobiliteit_vrouwen := as.numeric(TotaalVrouwen_66) + as.numeric(TotaalVrouwen_113)]
 
freq_dt <- mv_dt[RegioS == 'Nederland', .(RegioS, Perioden, Verhuismobiliteit_5,
                                          Verhuismobiliteit_mannen, Verhuismobiliteit_vrouwen, 
                                          Pop_total, Pop_mannen, Pop_vrouwen)]

freq_dt[, Freq_total := as.numeric(1000*as.numeric(Verhuismobiliteit_5) / as.numeric(Pop_total))]
freq_dt[, Freq_mannen := as.numeric(1000*as.numeric(Verhuismobiliteit_mannen) / as.numeric(Pop_mannen))]
freq_dt[, Freq_vrouwen := as.numeric(1000*as.numeric(Verhuismobiliteit_vrouwen) / as.numeric(Pop_vrouwen))]

setkey(bbp_dt, Perioden)
freq_dt[bbp_dt[!grepl('kwartaal', Perioden) & SoortGegevens == 'Werkelijke prijzen',],
        Bbp := as.numeric(i.BrutoBinnenlandsProduct_2)]

freq_dt <- melt(freq_dt, id.vars = c('RegioS', 'Perioden'))[, variable := gsub('5', 'total', variable)][, RegioS := NULL][, value := as.numeric(value)]

for (i in c('total', 'mannen', 'vrouwen')) {
  freq_dt[grepl(i, variable), gender := i ][, variable := gsub(paste0('_',i), '', variable)]
}


## plots
ggplot(freq_dt, aes(Perioden, value, group = gender)) +
  geom_line(aes(col = gender)) +
  facet_wrap(~variable, scales = "free") +
  theme_bw()


y <- datasets::Nile
x <- cbind(c(rep(0, 27), rep(1, length(y) - 27)))


y <- ts(freq_dt[variable == 'Freq' & gender == 'total', value], 
   start = c(1988, 1)
   )

adf.test(y, alternative = "stationary")


x <- ts(freq_dt[variable == 'Bbp' & !is.na(value), value], 
        start = c(1995, 1)
)


## time invariant regression
buildModReg <- function(v) {
  dV <- exp(v[1])
  dW <- c(exp(v[2]), 0)
  m0 <- v[3:4]
  dlmModReg(x, dV = dV, dW = dW, m0 = m0)
}

varGuess <- var(diff(y), na.rm = TRUE)
mu0Guess <- as.numeric(y[1])
lambdaGuess <- mean(diff(y), na.rm = TRUE)

parm <- c(log(varGuess), log(varGuess/5), mu0Guess,
          lambdaGuess)

## local-level trend model
buildModPoly2 <- function(v) {
  dV <- exp(v[1])
  dW <- exp(v[2:3])
  m0 <- v[4:5]
  dlmModPoly(order = 2, dV = dV, dW = dW, m0 = m0)
}

varGuess <- var(diff(y), na.rm = TRUE)
mu0Guess <- as.numeric(y[1])
lambda0Guess <- 0

parm <- c(log(varGuess), log(varGuess), log(varGuess),
          mu0Guess, lambda0Guess)
mle <- dlmMLE(y, parm = parm, buildModPoly2)
if (mle$convergence != 0) stop(mle$message)


## time variant regression
buildModReg <- function(v) {
  dV <- exp(v[1])
  dW <- exp(v[2:3])
  m0 <- v[4:5]
  dlmModReg(x_qu, dV = dV, dW = dW, m0 = m0)
}

varGuess <- var(diff(y_qu), na.rm = TRUE)
mu0Guess <- as.numeric(y_qu[1])
lambda0Guess <- mean(diff(y_qu), na.rm = TRUE)
parm <- c(log(varGuess), log(varGuess/5), log(varGuess/5),
          mu0Guess, lambda0Guess)

## model
mle <- dlmMLE(y_qu, parm = parm, build = buildModReg)
if (mle$convergence != 0) stop(mle$message)
model <- buildModReg(mle$par)

## filter
filt <- dlmFilter(y_qu, model)
## The final, filtered data is this linear
## combination
filtered <- filt$m[-1, 1] + x_qu* filt$m[-1, 2]
both <- cbind(y = y_qu, filtered = filtered)


ggplot(data.table(time(both), both), aes(x=V1)) + 
  geom_line(aes(y=y), col = 'black') + 
  geom_line(aes(y=filtered), col = 'red') +
  theme_bw()


smooth <- dlmSmooth(y, model)
## smooth$s contains the smoothed values
## The final, smoothed data is this linear
## combination
smoothed <- smooth$s[-1, 1] + x * smooth$s[-1, 2]
both <- cbind(y = y, smoothed = smoothed)

ggplot(data.table(time(both), both), aes(x=V1)) + 
  geom_line(aes(y=y), col = 'black') + 
  geom_line(aes(y=smoothed), col = 'red') +
  theme_bw()


# simulate AR(1) process
# phi = .8, sig2 = .25
nobs <- 250
yt <- arima.sim(n=nobs,list(ar=.8,ma=0),sd=.5)


# estimate AR(1) for comparison
model10 <- Arima(yt,order=c(1,0,0),method="ML",include.mean=FALSE)
model10

# set parameter restrictions 
parm_rest <- function(parm){
  return( c(exp(parm[1])/(1+exp(parm[1])),exp(parm[2])) ) 
}

# set up SS model
ssm1 <- function(parm){
  parm <- parm_rest(parm)
  return( dlm(FF=1,V=0,GG=parm[1],W=parm[2],
              m0=0,C0=solve(1-parm[1]^2)*parm[2]) )
}
# estimate parameters
fit1 <- dlmMLE(y=yt,parm=c(0,1),build=ssm1,hessian=T)

# get estimates 
coef <- parm_rest(fit1$par)
# get standard errors using delta method
dg1 <- exp(fit1$par[1])/(1+exp(fit1$par[1]))^2
dg2 <- exp(fit1$par[2])
dg <- diag(c(dg1,dg2))
var <- dg%*%solve(fit1$hessian)%*%dg
# print results
coef; sqrt(diag(var))




freq_vec <- freq_dt[variable == 'Freq' & gender == 'total', value]
bbp_vec <- freq_dt[!is.na(bbp), bbp]

names(freq_vec) <- seq(1995, 2016)
names(bbp_vec) <- seq(1995, 2016)




############################################################################
## time varying parameters model
############################################################################
## run example from http://lalas.github.io/quantitativeThoughts/r/2014/09/01/dlmTutorial.html and http://past.rinfinance.com/agenda/2012/workshop/Zivot+Yollin.pdf


## check linear model with non-varying parameters
freq_fit <- lm(freq_vec ~ bbp_vec)

# Specifying a set model parameters
s2_obs = 1      # Variance of observations
s2_alpha = 0.01 # Variance of the alpha regression parameter
s2_beta = 0.01  # Variance of the beta regression parameter
# Construct a regression model
tvp.dlm = dlmModReg(X=bbp_vec, addInt=TRUE, dV=s2_obs, dW=c(s2_alpha, s2_beta))

# looking at the various component
tvp.dlm[c("FF","V","GG","W","m0","C0")]
tvp.dlm[c("JFF","JV","JGG","JW")]

start.vals = c(0,0,0)
# Names ln variance of: observation y, alpha and beta (corresponding intercept and slope of y (freq) with respect to X (bbp))
names(start.vals) = c("lns2_obs", "lns2_alpha", "lns2_beta")

# function to build Time Varying Parameter state space model
buildTVP <- function(parm, x.mat){
    parm <- exp(parm)
  return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2], parm[3])) )
}

# Estimate the model
TVP.mle = dlmMLE(y=freq_vec, parm=start.vals, x.mat=bbp_vec, build=buildTVP, hessian=T)

# get sd estimates
se2 <- sqrt(exp(TVP.mle$par))
names(se2) = c("s_obs", "s_alpha", "s_beta")
sqrt(se2)

# Build fitted ss model, passing to it bbp_vec as the matrix X in the model
TVP.dlm <- buildTVP(TVP.mle$par, bbp_vec)

## filtering and smoothing

# filtering: optimal estimates of θ_t given information available at time t.
TVP.f <- dlmFilter(y = freq_vec, mod = TVP.dlm)
class(TVP.f)
names(TVP.f)

# smoothing: optimal estimates of θ_t given information available at time T.
TVP.s <- dlmSmooth(TVP.f)
class(TVP.s)
names(TVP.s)

## rownames need to be interpreted as date
rownames(TVP.s$s) <- c(0, paste0(as.character(seq(1995, 2016)), '-01-01'))
## Plotting the results (smoothed values)

# extract smoothed states - intercept and slope coefs
alpha.s = xts(TVP.s$s[-1,1,drop=FALSE], as.Date(rownames(TVP.s$s[-1,])))
beta.s  = xts(TVP.s$s[-1,2,drop=FALSE], as.Date(rownames(TVP.s$s[-1,])))
colnames(alpha.s) = "alpha"
colnames(beta.s)  = "beta"


# extract std errors - dlmSvd2var gives list of MSE matrices
mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
se.xts = xts(se.mat[-1, ], index(beta.s))
colnames(se.xts) = c("alpha", "beta")
a.u = alpha.s + 1.96*se.xts[, "alpha"]
a.l = alpha.s - 1.96*se.xts[, "alpha"]
b.u = beta.s  + 1.96*se.xts[, "beta"]
b.l = beta.s  - 1.96*se.xts[, "beta"]


alpha.df <- data.frame(dateTime = index(se.xts), alpha = alpha.s, upr = a.u, lwr = a.l)
names(alpha.df) <- c("dateTime", "alpha", "upr", "lwr")
beta.df  <- data.frame(dateTime = index(se.xts), beta = beta.s, upr = b.u, lwr = b.l)
names(beta.df) <- c("dateTime", "beta", "upr", "lwr")

## Plotting alpha
ggplot(data = alpha.df, aes(dateTime, alpha) ) +
    geom_point () +
    geom_line() +
    geom_ribbon(data=alpha.df, aes(ymin=lwr,ymax=upr), alpha=0.3) +
    labs(x = "year", y = expression(alpha), title = expression(paste("State Space Values of ", alpha, " over Time")))


## Plotting beta
ggplot(data = beta.df, aes(dateTime, beta) ) +
    geom_point (data = beta.df, aes(dateTime, beta) ) +
    geom_line() +
    geom_ribbon(data=beta.df , aes(ymin=lwr,ymax=upr), alpha=0.3) +
    labs(x = "year", y = expression(beta), title = expression(paste("State Space Values of ", beta, " over Time")))


# Construct add 10 missing values to end of sample
new.xts = xts(rep(NA, 10), seq.Date(from=end(HAM1), by="months", length.out=11)[-1])
# Add this NA data to the original y (HAM1) series
HAM1.ext = merge(HAM1, new.xts)[,1]
# Filter extended y (HAM1) series
TVP.ext.f = dlmFilter(HAM1.ext, TVP.dlm)
# extract h-step ahead forecasts of state vector
TVP.ext.f$m[as.character(index(new.xts)),]
