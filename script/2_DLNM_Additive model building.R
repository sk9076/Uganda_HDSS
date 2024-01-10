# Model building process

dat_comb_era2 <- rio::import(here::here("data", "Uganda", "clean data", "dat_comb_era2.xlsx"))
dat_comb <- dat_comb_era2
str(dat_comb)
dat_comb %<>% mutate_at(vars(year, week), as.numeric)

summary(dat_comb)
names(dat_comb)

# assign time sequence
dat_comb %<>% arrange(week_date, desc=F)
dat_comb$time = 1:nrow(dat_comb)

range_rain <-range(dat_comb$rain,na.rm=T)
range_tmin <-range(dat_comb$tmin,na.rm=T)
range_tmax <- range(dat_comb$tmax, na.rm=T)
range_temp <- range(dat_comb$temp, na.rm=T)

lag <- c(0,12) # 0 week to 12 weeks


## Best performing model for tmax
var_knots_rain = equal_knots(dat_comb$rain, 4)
var_knots_tmin = equal_knots(dat_comb$tmin, 3)
var_knots_tmax = equal_knots(dat_comb$tmax, 3)
var_knots_temp = equal_knots(dat_comb$temp, 3)
var_spline_rain = "ns"
var_spline_tmax = "ns"
var_spline_tmin = "ns"
var_spline_temp = "ns"
lag_knots = c(8)
time_df = 2*5


# BASIS FOR rainfall (ns, equal knots using 3 knots)
bound_rain <- range_rain
argvar_rain <- list(type= "ns",
                    knots= var_knots_rain,
                    cen=0
)

arglag_rain <- list(type="ns",
                    knots = lag_knots,
                    df=6)

##########*
# BASIS FOR Minimum TEMPERATURE (natural spline, equal space knots - the first 2 out of 3)
bound_tmin <- range_tmin

argvar_tmin <- list(type= var_spline_tmin,
                    knots= var_knots_tmin,
                    cen=0
)

arglag_tmin <- list(type= "ns",
                    knots = lag_knots
)

##########*
# BASIS FOR Maximum TEMPERATURE (ns, equal splits knots with 3 knots)
bound_tmax <- range_tmax

argvar_tmax <- list(type= var_spline_tmax,
                    knots= var_knots_tmax,
                    cen=0
)

arglag_tmax <- list(type= "ns",
                    knots = lag_knots
)

###########
# BASIS FOR average temperature (ns, equal splits knots with 3 knots)
bound_temp <- range_temp

argvar_temp <- list(type= var_spline_temp,
                    knots= var_knots_temp,
                    cen=0
)

arglag_temp <- list(type= "ns",
                    knots = lag_knots
)

## Build all CBs
suppressWarnings({
  cb_tmin <- crossbasis(dat_comb$tmin, lag= lag, argvar = argvar_tmin, arglag = arglag_tmin)
  cb_tmax <- crossbasis(dat_comb$tmax, lag= lag, argvar = argvar_tmax, arglag = arglag_tmax)
  cb_temp <- crossbasis(dat_comb$temp, lag= lag, argvar = argvar_temp, arglag = arglag_temp)
  cb_rain <- crossbasis(dat_comb$rain, lag = lag, argvar = argvar_rain, arglag = arglag_rain)
})

# Tmax only
model_tmax <- glm(n ~ cb_tmax +
               ns(time, df=time_df), # time as natural spline (df captures the season)
             family=quasipoisson(),
             dat_comb)

fqaic(model_tmax)
fqbic(model_tmax)
draw_lag_plot2(model=model_tmax,
               dat = dat_comb,
               tmax = T,
               cb_tmax = cb_tmax,
               cen_tmax = min(dat_comb$tmax)
)

# Tmin only
model_tmin <- glm(n ~ cb_tmin +
                    ns(time, df=time_df), # time as natural spline (df captures the season)
                  family=quasipoisson(),
                  dat_comb)

fqaic(model_tmin)
fqbic(model_tmin)
draw_lag_plot2(model=model_tmin,
               dat = dat_comb,
               tmin = T,
               cb_tmin = cb_tmin,
               cen_tmin = min(dat_comb$tmin)
               )

# Rain only
model_rain <- glm(n ~ cb_rain +
                    ns(time, df=time_df), # time as natural spline (df captures the season)
                  family=quasipoisson(),
                  dat_comb)

fqaic(model_rain)
fqbic(model_rain)
draw_lag_plot2(model=model_rain,
               dat = dat_comb,
               rain = T,
               cb_rain = cb_rain,
               cen_rain = min(dat_comb$rain))

## Temp only
model_temp <- glm(n ~ cb_temp +
                    ns(time, df=time_df), # time as natural spline (df captures the season)
                  family=quasipoisson(),
                  dat_comb)

fqaic(model_temp)
fqbic(model_temp)
draw_lag_plot2(model=model_temp,
               dat = dat_comb,
               temp = T,
               cb_temp = cb_temp,
               cen_temp = min(dat_comb$temp)
)

# Tmin + Rain
model_rain_tmin <- glm(n ~ cb_rain + cb_tmin+
                    ns(time, df=time_df), # time as natural spline (df captures the season)
                  family=quasipoisson(),
                  dat_comb)

fqaic(model_rain_tmin)
fqbic(model_rain_tmin)

draw_lag_plot2(model=model_rain_tmin,
               dat = dat_comb,
               rain = T, tmin=T,
               cb_rain = cb_rain, cb_tmin = cb_tmin,
               cen_rain = min(dat_comb$rain),
               cen_tmin = 20) # this is usually the threshold temperature.

# Tmax + Rain
model_rain_tmax <- glm(n ~ cb_rain + cb_tmax+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb)

fqaic(model_rain_tmax)
fqbic(model_rain_tmax)

draw_lag_plot2(model=model_rain_tmax,
              dat = dat_comb,
              rain = T, tmax=T,
              cb_rain = cb_rain, cb_tmax = cb_tmax,
              cen_rain = min(dat_comb$rain),
              cen_tmax = min(dat_comb$tmax))

# Temp + Rain
model_rain_temp <- glm(n ~ cb_rain + cb_temp+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb)

fqaic(model_rain_temp)
fqbic(model_rain_temp)

draw_lag_plot2(model=model_rain_temp,
               dat = dat_comb,
               rain = T, temp=T,
               cb_rain = cb_rain, cb_temp = cb_temp,
               cen_rain = min(dat_comb$rain),
               cen_temp = min(dat_comb$temp))


# Tmin + Tmax + Rain
model_rain_tmin_tmax <- glm(n ~ cb_rain + cb_tmax+ cb_tmin+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb)

fqaic(model_rain_tmin_tmax)
fqbic(model_rain_tmin_tmax)
draw_lag_plot2(model=model_rain_tmin_tmax,
               dat = dat_comb,
               rain = T, tmin=T, tmax =T,
               cb_rain = cb_rain, cb_tmin = cb_tmin, cb_tmax = cb_tmax,
               cen_rain = min(dat_comb$rain),
               cen_tmin = min(dat_comb$tmin),
               cen_tmax = min(dat_comb$tmax)
               )


# Tmin + Temp + Rain
model_rain_tmin_temp <- glm(n ~ cb_rain + cb_temp+ cb_tmin+
                              ns(time, df=time_df), # time as natural spline (df captures the season)
                            family=quasipoisson(),
                            dat_comb)

fqaic(model_rain_tmin_temp)
fqbic(model_rain_tmin_temp)

draw_lag_plot2(model=model_rain_tmin_temp,
               dat = dat_comb,
               rain = T, tmin=T, temp =T,
               cb_rain = cb_rain, cb_tmin = cb_tmin, cb_temp = cb_temp,
               cen_rain = min(dat_comb$rain),
               cen_tmin = min(dat_comb$tmin),
               cen_temp = min(dat_comb$temp)
)

# Tmax + Temp + Rain
model_rain_tmax_temp <- glm(n ~ cb_rain + cb_temp+ cb_tmax+
                              ns(time, df=time_df), # time as natural spline (df captures the season)
                            family=quasipoisson(),
                            dat_comb)

fqaic(model_rain_tmax_temp)
fqbic(model_rain_tmax_temp)

draw_lag_plot2(model=model_rain_tmax_temp,
               dat = dat_comb,
               rain = T, tmax=T, temp =T,
               cb_rain = cb_rain, cb_tmax = cb_tmax, cb_temp = cb_temp,
               cen_rain = min(dat_comb$rain),
               cen_tmax = min(dat_comb$tmax),
               cen_temp = min(dat_comb$temp)
)


# Tmax + Tmin + Temp + Rain
# Tmax + Temp + Rain
model_full <- glm(n ~ cb_rain + cb_temp+ cb_tmax+ cb_tmin+
                              ns(time, df=time_df), # time as natural spline (df captures the season)
                            family=quasipoisson(),
                            dat_comb)

fqaic(model_full)
fqbic(model_full)

draw_lag_plot2(model=model_full,
               dat = dat_comb,
               rain = T, tmin = T, tmax=T, temp =T,
               cb_rain = cb_rain, cb_tmin = cb_tmin, cb_tmax = cb_tmax, cb_temp = cb_temp,
               cen_rain = min(dat_comb$rain),
               cen_tmin = min(dat_comb$tmin),
               cen_tmax = min(dat_comb$tmax),
               cen_temp = min(dat_comb$temp)
)
