#  Predict soil mosture from precipitation records
#  Applied to Maxfield Meadow Ipomopsis experiment
#  Model: Pan 2003, 2012
#  Optimization by genetic algorithm: Coopersmith et al. 2014, 2015, 2015
#  John Powers 2020-01-17

#TODO look at how continuous soil moisture data from Kettle Ponds S of Gothic, CO and Judd Falls compares to precip

library(tidyverse)
library(profvis)

load("data/maxfield_data.rda")
# data to work with
daily_all %>% select(date, EPA_NOAA_filled)
sm.subplot %>% select(date, snow, water, plotid, VWC)
waterdates
meltdates %>% select(year, plot, sun_date)
groundcover %>% select(first_snow)

# list of days that water was added to subplots
watering_dates <- waterdates %>% select(-day) %>% mutate(date=as.Date(date)) %>% 
  pivot_wider(names_from=precip_treatments, values_from=date) %>% 
  mutate(year=factor(year)) %>% rename(water_begin=started, water_end=ended) %>%
  mutate(watering_dates = map2(water_begin, water_end, seq, "2 days")) %>%
  unnest(watering_dates) %>% pull(watering_dates)

# list of indices since 2018 that were melt dates for each plot
sun_indices <- map(setNames(1:6,as.character(1:6)), ~ meltdates %>% filter(plot == .x) %>% 
                     mutate(sun_date=as.integer(as.Date(parse_date_time(paste0(year, "-", sun_date), orders="yj")) - 
                                                as.Date("2018-01-01"))) %>% pull(sun_date))

# make a fast vector of the precipitation data indexed by days since 2018 started
daily_precip_est <- daily_all %>% select(year, date, EPA_NOAA_filled, ground_covered) %>% 
  filter(year %in% 2018:2020) %>% arrange(date) %>% 
  mutate(year = fct_drop(year),
         droughted = findInterval(date, as.Date(waterdates$date)) %% 2 == 1,
         EPA_NOAA_filled  = replace_na(EPA_NOAA_filled, 0), #  fill a a few missing precip days with 0 mm - NOAA_billy data missing 
         est_Reduction = EPA_NOAA_filled * ifelse(droughted, 0.5, 1),
         est_Control   = EPA_NOAA_filled,
         est_Addition  = EPA_NOAA_filled + ifelse(date %in% watering_dates, 3.5, 0))

# average VWC across subplots for each water treatment and plot on each date
sm.waterplot <- sm.subplot %>% group_by(year, date, water, plot) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop") %>% 
  mutate(start = as.integer(as.Date(date) - daily_precip_est$date[1]), water=as.character(water)) %>% 
  filter(!is.nan(VWC), water=="Control")# drop 2019-07-29 plot 1 control subplots and try without water addition plots

# the moisture loss coefficient as a function of the day of year (ith day before measurement)
eta <- function(start, i, c_eta) {
  c_eta[1]  + c_eta[2] * sin(2 * pi * ((start - i) %% 365 + c_eta[3]) / 365)
}

# get the decayed precip contribution from the ith day before the morning soil measurement (i=1 means day prior)
# add big precip (20 mm) on sun_date to account for melt
# with depth of z = 120 mm / 2 (halfway down 12 cm probe - check this)
beta <- function(start, i, z, c_eta, water, plot) {
  meltwater <- 0
  if((start - i) %in% sun_indices[[plot]]) meltwater <- 20
  precip_term <- ((meltwater + daily_precip_est[[start-i, paste("est", as.character(water), sep="_")]]) / 
                    eta(start, i, c_eta)) * (1 - exp(-eta(start, i, c_eta)/z))
  ifelse(i == 1, precip_term, 
                 precip_term * exp(-sum(map_dbl(1:(i-1), eta, start=start, c_eta=c_eta)/z)))
}

betasum <- function(start, n, z, c_eta, water, plot, checkstart = TRUE, ...) {
  if(checkstart) {
    # can't give a value for dates before dataset starts
    if(start - n < 1) return(NA)
    # don't calculate betasum if n-1 window reaches before snowmelt date in plot
    yr <- as.character(daily_precip_est$year[start])
    plt <- plot
    plot_sun_date <- meltdates %>% filter(year==yr, plot == plt) %>% 
      pull(sun_date) %>% paste0("-",yr) %>% parse_date_time(orders="jy") %>% as.Date
    if(daily_precip_est$date[start] < plot_sun_date) return(NA) #TODO add back start - n
    # don't calculate betasum if start date is after first snow in fall
    first_snow_date <- groundcover %>% bind_rows(data.frame(first_snow=first_snow_2020[1])) %>% 
      filter(year(first_snow)==yr) %>% pull(first_snow) %>% as.Date
    if(daily_precip_est$date[start] >= first_snow_date) return(NA)
  }
  # sum the decayed precip contributions from the 1st to n-1th day before the measurement (n-1 days of precip)
  sum(map_dbl(1:(n-1), beta, start=start, z=z, c_eta=c_eta, water=water, plot=plot))
}

#  estimate of soil moisture at day "start" by converting betasum to VWC scale
theta_est = function(start, n, z, c_eta, c_theta, water, plot, checkstart=TRUE, ...) {
  c_theta[1] + (c_theta[2] - c_theta[1]) * (1 - exp(-c_theta[3] * betasum(start=start, n=n, z=z, c_eta=c_eta, water=water, plot=plot, checkstart=checkstart)))
}

# fitness function for c_eta coefficient = correlation between betasum and observed VWC
betasum_VWC_cor <- function(n, z, c_eta){
  betasum <- pmap_dbl(sm.waterplot, betasum, n=n, z=z, c_eta=c_eta, checkstart=F)
  return(list(cor=cor(betasum, sm.waterplot$VWC), betasum=betasum))
}

#fitness function for c_theta coefficients = RMSE between predicted and observed VWC
theta_VWC_rms <- function(n, z, c_eta, c_theta){
  VWC_pred <- pmap_dbl(sm.waterplot, theta_est, n=n, z=z, c_eta=c_eta, c_theta=c_theta, checkstart=F)
  return(list(rms=sqrt(sum((VWC_pred-sm.waterplot$VWC/100)^2)/length(VWC_pred)), VWC_pred=VWC_pred))
}

# genetic algorithm to maximize betasum_VWC_cor by optimizing c_eta parameters
# c1 = mean value of eta (mm/day), c2 = magnitude of variation(mm/day), c3 = phase (day)
# searching domain = {0 < c1 <  20 mm/day; 0 < c2 ≤ c1; 0 < c3 < 366}
library(GA)
GA.betasum <- ga(type = "real-valued", lower = c(0, 0, 0), upper = c(20, 20, 366), 
         fitness =  function(c_eta) betasum_VWC_cor(n=15, z=60, c_eta=c_eta)$cor - ifelse(c_eta[2]>c_eta[1],10,0), #subtract a penalty to force 0<c2≤c1
         popSize = 25, maxiter = 1, run = 1)
summary(GA.betasum)
daily_precip_est$eta <- map_dbl(1:nrow(daily_precip_est), eta, c_eta=GA.betasum@solution, start=1001)

#the eta solution is very flat through time and at its minimum (should be maximal in the summer) - simpler to just fix the eta value.
eta_optim_df <- tibble(c1=seq(from=13.8,to=14.2,by=0.05), cor=map_dbl(c1, ~ betasum_VWC_cor(n=15, z=60, c_eta=c(.x,0,0))$cor))
ggplot(eta_optim_df, aes(x=c1, y=cor))+ geom_point()+ geom_smooth(se=F) #eta= c1 = 13.9 mm/day is best loss parameter with n=15, z=60
c_eta_optim <- c(13.9, 0, 0)

# find the optimum window size n given a fixed eta - this is not what is described in Coopersmith papers
#n_optim_df <- tibble(n=seq(from=1,to=20,by=1), cor=map_dbl(n, ~ betasum_VWC_cor(n=.x, z=60, c_eta=c_eta_optim)$cor))
#ggplot(n_optim_df, aes(x=n, y=cor))+ geom_point()+ geom_line()

# optimize both n and eta at once
#n_eta_optim_df <- tibble(n=rep(seq(from=9,to=20,by=1), 12), c1=rep(seq(from=9,to=20,by=1), each=12), 
#                         cor=map2_dbl(n, c1, ~ betasum_VWC_cor(n=.x, z=60, c_eta=c(.y,0,0))$cor))
#ggplot(n_eta_optim_df, aes(x=n, y=c1, fill=cor))+ geom_tile() + geom_text(aes(label=round(cor,3))) + scale_fill_viridis_c()
#ggplot(n_eta_optim_df, aes(x=n, y=cor, color=c1))+ geom_point()+ geom_line(aes(group=c1))

# find suitable window size by correlations with large n - per Coopersmith
betasum.large_n <- betasum_VWC_cor(n=80, z=60, c_eta=c_eta_optim)$betasum
n_limit_df <- tibble(n=seq(from=1,to=20,by=1), betacor = map(n, ~ betasum_VWC_cor(n=.x, z=60, c_eta=c_eta_optim)$betasum) %>% 
                       map_dbl(cor, betasum.large_n)) 
ggplot(n_limit_df, aes(x=n, y=betacor))+ geom_point()+ geom_line() #n=15 has correlation >0.999 with n=80

# genetic algorithm to minimize theta_VWC_rms with c_theta
# c_theta[1] = theta_re = residual soil moisture
# c_theta[2] = phi_e = porosity
# c_theta[3] = c4 = soil hydraulic properties (conductance/drainage)
GA.theta <- ga(type = "real-valued", lower = c(0, 0, 0, 0), upper = c(0.03, 0.7, 5, 20), 
                 fitness =  function(x) -theta_VWC_rms(n=15, z=60, c_eta=c(x[4],0,0), c_theta=c(x[1],x[2],x[3]))$rms,
                 popSize = 25, maxiter = 1, run = 1)
summary(GA.theta)
c_theta_optim=GA.theta@solution[1:3]#c(0.02, 0.4, 1.2)
c_eta_optim=c(GA.theta@solution[4],0,0)#c(13.9,0,0)

# now that model is fit, predict soil moisture at each measuring date
sm.waterplot$VWC_pred <- theta_VWC_rms(n=15, z=60, c_eta=c_eta_optim, c_theta=c_theta_optim)$VWC_pred

ggplot(sm.waterplot, aes(x=date, color=water))+ facet_wrap(vars(year), scales="free_x") + scale_color_manual(values=water_pal) +
  geom_line(aes(y=VWC_pred)) +
  geom_point(aes(y=VWC/100)) + geom_smooth(aes(y=VWC/100), se=F)

ggplot(sm.waterplot, aes(y=VWC_pred, x=VWC/100, color=water))+ facet_wrap(vars(year), scales="free_x") + geom_point() + scale_color_manual(values=water_pal) + geom_smooth(method="lm", se=F) + geom_abline(intercept = 0, slope=1)

# predict soil moisture throughout the year
daily_precip_est$VWC_Control_4 <- map_dbl(1:nrow(daily_precip_est), 
                                           theta_est, n=15, z=60, c_eta= c_eta_optim, c_theta=c_theta_optim, 
                                           water="Control",plot="4")
ggplot(daily_precip_est, aes(x=date))+ facet_wrap(vars(year), scales="free_x", ncol=1) + 
  geom_line(aes(y=VWC_Control_3, color=droughted, group=1)) + scale_color_manual(values=c("black","blue"))+
  geom_col(data=daily_precip_est, aes(y=est_Control/100, fill=ground_covered))
  #geom_line(aes(y=eta/100, color=droughted, group=1))

# Look at just the summer and overlay actual soil moisture
ggplot(sm.waterplot, aes(x=yday(date), y=VWC/100, color=plot))+ facet_wrap(vars(year), ncol=1) +
  geom_col(data=daily_precip_est, aes(y=est_Control/100),color=NA, fill="grey20")+
  geom_line(data=daily_precip_est, aes(y=VWC_Control_4, color="4")) +
  geom_line(data=daily_precip_est, aes(y=VWC_Control_3, color="3")) +
  geom_point() + geom_line(size=1) + coord_cartesian(xlim=c(120,240)) + 
  geom_vline(data=meltdates, aes(xintercept=sun_date, color=plot), size=2)+
  labs(x="Day of year", y="Soil moisture", color="Plot", title=paste(paste(unique(sm.waterplot$water), collapse=" "), "treatment"))

#Show all water treatments
ggplot(sm.waterplot, aes(x=as.Date(date), color=water))+ facet_wrap(vars(year), scales="free_x") + scale_color_manual(values=water_pal) +
  geom_line(data=daily_precip_est, aes(y=VWC_Control_3, color="Control")) +
  geom_point(aes(y=VWC/100)) + geom_line(aes(y=VWC/100, group=plot))

save(daily_precip_est, file="data/daily_precip_est.rda")
