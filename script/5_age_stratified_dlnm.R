dat_comb_str <- rio::import(here::here("data", "Uganda", "clean data", "dat_comb_era2_str.xlsx"))

dat_comb_str %<>% mutate_at(vars(year, week), as.numeric)

summary(dat_comb_str)
names(dat_comb_str)

# assign time sequence
dat_comb_str %<>% arrange(week_date, desc=F)
dat_comb_str$time = 1:nrow(dat_comb_str)

range_rain <-range(dat_comb_str$rain,na.rm=T)
range_tmin <-range(dat_comb_str$tmin,na.rm=T)
range_tmax <- range(dat_comb_str$tmax, na.rm=T)
range_temp <- range(dat_comb_str$temp, na.rm=T)

lag <- c(0,12) # 0 week to 12 weeks


## Best performing model for tmax
var_knots_rain = equal_knots(dat_comb_str$rain, 4)
var_knots_tmin = equal_knots(dat_comb_str$tmin, 3)
var_knots_tmax = equal_knots(dat_comb_str$tmax, 3)
var_knots_temp = equal_knots(dat_comb_str$temp, 3)
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
  cb_tmin <- crossbasis(dat_comb_str$tmin, lag= lag, argvar = argvar_tmin, arglag = arglag_tmin)
  cb_tmax <- crossbasis(dat_comb_str$tmax, lag= lag, argvar = argvar_tmax, arglag = arglag_tmax)
  cb_temp <- crossbasis(dat_comb_str$temp, lag= lag, argvar = argvar_temp, arglag = arglag_temp)
  cb_rain <- crossbasis(dat_comb_str$rain, lag = lag, argvar = argvar_rain, arglag = arglag_rain)
})

############ School Age Children ################
# Tmin + Rain
model_rain_tmin_sac <- glm(School_age_children ~ cb_rain + cb_tmin+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb_str)

fqaic(model_rain_tmin_sac)
fqbic(model_rain_tmin_sac)

draw_lag_plot2(model=model_rain_tmin_sac,
               dat = dat_comb_str,
               rain = T, tmin=T,
               cb_rain = cb_rain, cb_tmin = cb_tmin,
               cen_rain = min(dat_comb_str$rain),
               cen_tmin = min(dat_comb_str$tmin))

# Tmax + Rain
model_rain_tmax_sac <- glm(School_age_children ~ cb_rain + cb_tmax+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb_str)

fqaic(model_rain_tmax_sac)
fqbic(model_rain_tmax_sac)

draw_lag_plot2(model=model_rain_tmax_sac,
               dat = dat_comb_str,
               rain = T, tmax=T,
               cb_rain = cb_rain, cb_tmax = cb_tmax,
               cen_rain = min(dat_comb_str$rain),
               cen_tmax = min(dat_comb_str$tmax))


############## Under 5 ######################
# Tmin + Rain
model_rain_tmin_u5 <- glm(Under_5 ~ cb_rain + cb_tmin+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb_str)

fqaic(model_rain_tmin_u5)
fqbic(model_rain_tmin_u5)

draw_lag_plot2(model=model_rain_tmin_u5,
               dat = dat_comb_str,
               rain = T, tmin=T,
               cb_rain = cb_rain, cb_tmin = cb_tmin,
               cen_rain = min(dat_comb_str$rain),
               cen_tmin = min(dat_comb_str$tmin))

# Tmax + Rain
model_rain_tmax_u5 <- glm(Under_5 ~ cb_rain + cb_tmax+
                         ns(time, df=time_df), # time as natural spline (df captures the season)
                       family=quasipoisson(),
                       dat_comb_str)

fqaic(model_rain_tmax_u5)
fqbic(model_rain_tmax_u5)

draw_lag_plot2(model=model_rain_tmax_u5,
               dat = dat_comb_str,
               rain = T, tmax=T,
               cb_rain = cb_rain, cb_tmax = cb_tmax,
               cen_rain = min(dat_comb_str$rain),
               cen_tmax = min(dat_comb_str$tmax))




############## Others ######################
# Tmin + Rain
model_rain_tmin_others <- glm(Others ~ cb_rain + cb_tmin+
                            ns(time, df=time_df), # time as natural spline (df captures the season)
                          family=quasipoisson(),
                          dat_comb_str)

fqaic(model_rain_tmin_others)
fqbic(model_rain_tmin_others)

draw_lag_plot2(model=model_rain_tmin_others,
               dat = dat_comb_str,
               rain = T, tmin=T,
               cb_rain = cb_rain, cb_tmin = cb_tmin,
               cen_rain = min(dat_comb_str$rain),
               cen_tmin = min(dat_comb_str$tmin))

# Tmax + Rain
model_rain_tmax_others <- glm(Others ~ cb_rain + cb_tmax+
                            ns(time, df=time_df), # time as natural spline (df captures the season)
                          family=quasipoisson(),
                          dat_comb_str)

fqaic(model_rain_tmax_others)
fqbic(model_rain_tmax_others)

draw_lag_plot2(model=model_rain_tmax_others,
               dat = dat_comb_str,
               rain = T, tmax=T,
               cb_rain = cb_rain, cb_tmax = cb_tmax,
               cen_rain = min(dat_comb_str$rain),
               cen_tmax = min(dat_comb_str$tmax))

########## Visualization
cen_rain <- bound_rain[1] # or 21? 80?
cen_temp <- bound_temp[1]
cen_tmax <- bound_tmax[1]
cen_tmin <- bound_tmin[1]


# Association between the rainfall and malaria incidence
cp_rain_max_sac <- crosspred(cb_rain,
                     model_rain_tmax_sac,
                     from=bound_rain[1],
                     to=bound_rain[2], by=1,
                     cen=cen_rain)

cp_rain_min_sac <- crosspred(cb_rain,
                             model_rain_tmin_sac,
                             from=bound_rain[1],
                             to=bound_rain[2], by=1,
                             cen=cen_rain)

cp_tmax_sac <- crosspred(cb_tmax,
                     model_rain_tmax_sac,
                     from=bound_tmax[1],
                     to=bound_tmax[2], by=.1,
                     cen=cen_tmax)

cp_tmin_sac <- crosspred(cb_tmin,
                         model_rain_tmin_sac,
                         from=bound_tmin[1],
                         to=bound_tmin[2], by=.1,
                         cen=cen_tmin)

cp_rain_max_u5 <- crosspred(cb_rain,
                             model_rain_tmax_u5,
                             from=bound_rain[1],
                             to=bound_rain[2], by=1,
                             cen=cen_rain)

cp_rain_min_u5 <- crosspred(cb_rain,
                             model_rain_tmin_u5,
                             from=bound_rain[1],
                             to=bound_rain[2], by=1,
                             cen=cen_rain)

cp_tmax_u5 <- crosspred(cb_tmax,
                         model_rain_tmax_u5,
                         from=bound_tmax[1],
                         to=bound_tmax[2], by=.1,
                         cen=cen_tmax)

cp_tmin_u5 <- crosspred(cb_tmin,
                         model_rain_tmin_u5,
                         from=bound_tmin[1],
                         to=bound_tmin[2], by=.1,
                         cen=cen_tmin)

cp_rain_max_others <- crosspred(cb_rain,
                            model_rain_tmax_others,
                            from=bound_rain[1],
                            to=bound_rain[2], by=1,
                            cen=cen_rain)

cp_rain_min_others <- crosspred(cb_rain,
                            model_rain_tmin_others,
                            from=bound_rain[1],
                            to=bound_rain[2], by=1,
                            cen=cen_rain)

cp_tmax_others <- crosspred(cb_tmax,
                        model_rain_tmax_others,
                        from=bound_tmax[1],
                        to=bound_tmax[2], by=.1,
                        cen=cen_tmax)

cp_tmin_others <- crosspred(cb_tmin,
                        model_rain_tmin_others,
                        from=bound_tmin[1],
                        to=bound_tmin[2], by=.1,
                        cen=cen_tmin)


### Tmin at Lag 12
lag_12_tmin_sac <- crossreduce(cb_tmin,model_rain_tmin_sac,type="lag",value=12, from=bound_tmin[1],
                           to=bound_tmin[2],bylag=0.1,cen=cen_tmin)
lag_12_tmin_u5 <- crossreduce(cb_tmin,model_rain_tmin_u5,type="lag",value=12, from=bound_tmin[1],
                               to=bound_tmin[2],bylag=0.1,cen=cen_tmin)
lag_12_tmin_others <- crossreduce(cb_tmin,model_rain_tmin_others,type="lag",value=12, from=bound_tmin[1],
                               to=bound_tmin[2],bylag=0.1,cen=cen_tmin)

png(here::here("results", "stratified", "Lag12_tmin_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_12_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"),
     cex.lab = 1.2
     )
lines(lag_12_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Aross all age-group association \n(Lag = 12 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_12_tmin_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmin_u5,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Children under 5 years\n(Lag = 12 weeks)",cex=.8)

plot(lag_12_tmin_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmin_sac,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="School-age Children\n(Lag = 12 weeks)",cex=.8)

plot(lag_12_tmin_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmin_others,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Others\n(Lag = 12 weeks)",cex=.8)
dev.off()


### Tmax at Lag 0
lag_0_tmax_sac <- crossreduce(cb_tmax,model_rain_tmax_sac,type="lag",value=0, from=bound_tmax[1],
                               to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_0_tmax_u5 <- crossreduce(cb_tmax,model_rain_tmax_u5,type="lag",value=0, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_0_tmax_others <- crossreduce(cb_tmax,model_rain_tmax_others,type="lag",value=0, from=bound_tmax[1],
                                  to=bound_tmax[2],bylag=0.1,cen=cen_tmax)



### Tmax at Lag 12
lag_12_tmax_sac <- crossreduce(cb_tmax,model_rain_tmax_sac,type="lag",value=12, from=bound_tmax[1],
                               to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_12_tmax_u5 <- crossreduce(cb_tmax,model_rain_tmax_u5,type="lag",value=12, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_12_tmax_others <- crossreduce(cb_tmax,model_rain_tmax_others,type="lag",value=12, from=bound_tmax[1],
                                  to=bound_tmax[2],bylag=0.1,cen=cen_tmax)

png(here::here("results", "stratified", "Lag12_tmax_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_12_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"),
     cex.lab = 1.2
)
lines(lag_12_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Aross all age-group association \n(Lag = 12 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_12_tmax_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmax_u5,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Children under 5 years\n(Lag = 12 weeks)",cex=.8)

plot(lag_12_tmax_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmax_sac,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="School-age Children\n(Lag = 12 weeks)",cex=.8)

plot(lag_12_tmax_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,1.2),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmax_others,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Others\n(Lag = 12 weeks)",cex=.8)
dev.off()


### Tmax at Lag 0
lag_0_tmax_sac <- crossreduce(cb_tmax,model_rain_tmax_sac,type="lag",value=0, from=bound_tmax[1],
                              to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_0_tmax_u5 <- crossreduce(cb_tmax,model_rain_tmax_u5,type="lag",value=0, from=bound_tmax[1],
                             to=bound_tmax[2],bylag=0.1,cen=cen_tmax)
lag_0_tmax_others <- crossreduce(cb_tmax,model_rain_tmax_others,type="lag",value=0, from=bound_tmax[1],
                                 to=bound_tmax[2],bylag=0.1,cen=cen_tmax)

png(here::here("results", "stratified", "Lag0_tmax_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_0_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,3.5),xlab=expression("Maximum Temperature ("*~degree*C*")"),
     cex.lab = 1.2
)
lines(lag_0_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Aross all age-group association \n(Lag = 0 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_0_tmax_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,3.5),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_0_tmax_u5,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Children under 5 years\n(Lag = 0 weeks)",cex=.8)

plot(lag_0_tmax_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,3.5),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_0_tmax_sac,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="School-age Children\n(Lag = 0 weeks)",cex=.8)

plot(lag_0_tmax_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,3.5),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_0_tmax_others,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Others\n(Lag = 0 weeks)",cex=.8)
dev.off()


## Rain at Lag 0 (using Rain tmax model)
lag_0_rain <- crossreduce(cb_rain,model_rain_tmax,type="lag",value=0, from=bound_rain[1],
                          to=bound_rain[2],bylag=0.2,cen=cen_rain)

lag_0_rain_sac <- crossreduce(cb_rain,model_rain_tmax_sac,type="lag",value=0, from=bound_rain[1],
                               to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_0_rain_u5 <- crossreduce(cb_rain,model_rain_tmax_u5,type="lag",value=0, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_0_rain_others <- crossreduce(cb_rain,model_rain_tmax_others,type="lag",value=0, from=bound_rain[1],
                                  to=bound_rain[2],bylag=0.1,cen=cen_rain)


png(here::here("results", "stratified", "Lag0_rain_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_0_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2
)
lines(lag_0_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Aross all age-group association \n(Lag = 0 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_0_rain_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_0_rain_u5,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Children under 5 years\n(Lag = 0 weeks)",cex=.8)

plot(lag_0_rain_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_0_rain_sac,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="School-age Children\n(Lag = 0 weeks)",cex=.8)

plot(lag_0_rain_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_0_rain_others,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Others\n(Lag = 0 weeks)",cex=.8)
dev.off()

## Rain at Lag 2 (using Rain tmax model)
lag_2_rain <- crossreduce(cb_rain,model_rain_tmax,type="lag",value=2, from=bound_rain[1],
                          to=bound_rain[2],bylag=0.2,cen=cen_rain)

lag_2_rain_sac <- crossreduce(cb_rain,model_rain_tmax_sac,type="lag",value=2, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_2_rain_u5 <- crossreduce(cb_rain,model_rain_tmax_u5,type="lag",value=2, from=bound_rain[1],
                             to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_2_rain_others <- crossreduce(cb_rain,model_rain_tmax_others,type="lag",value=2, from=bound_rain[1],
                                 to=bound_rain[2],bylag=0.1,cen=cen_rain)


png(here::here("results", "stratified", "Lag2_rain_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_2_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2
)
lines(lag_2_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Aross all age-group association \n(Lag = 2 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_2_rain_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_2_rain_u5,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Children under 5 years\n(Lag = 2 weeks)",cex=.8)

plot(lag_2_rain_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_2_rain_sac,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="School-age Children\n(Lag = 2 weeks)",cex=.8)

plot(lag_2_rain_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_2_rain_others,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Others\n(Lag = 2 weeks)",cex=.8)
dev.off()


## Rain at Lag 4 (using Rain tmax model)
lag_4_rain <- crossreduce(cb_rain,model_rain_tmax,type="lag",value=4, from=bound_rain[1],
                          to=bound_rain[2],bylag=0.2,cen=cen_rain)

lag_4_rain_sac <- crossreduce(cb_rain,model_rain_tmax_sac,type="lag",value=4, from=bound_rain[1],
                              to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_4_rain_u5 <- crossreduce(cb_rain,model_rain_tmax_u5,type="lag",value=4, from=bound_rain[1],
                             to=bound_rain[2],bylag=0.1,cen=cen_rain)
lag_4_rain_others <- crossreduce(cb_rain,model_rain_tmax_others,type="lag",value=4, from=bound_rain[1],
                                 to=bound_rain[2],bylag=0.1,cen=cen_rain)


png(here::here("results", "stratified", "Lag4_rain_stratified.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)

layout(mat = matrix(c(1,2,1,3,1,4), nrow=2, ncol=3),
       heights = c(2,1.5),
       widths = c(2,2,2))
par(mar=c(6,12,3,10))

plot(lag_4_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2
)
lines(lag_4_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Aross all age-group association \n(Lag = 4 weeks)",cex=1)

par(mar=c(4,4,4,2))
plot(lag_4_rain_u5,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_4_rain_u5,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Children under 5 years\n(Lag = 4 weeks)",cex=.8)

plot(lag_4_rain_sac,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_4_rain_sac,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="School-age Children\n(Lag = 4 weeks)",cex=.8)

plot(lag_4_rain_others,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.4,2.5),xlab=expression("Rainfall (mm)"),
     cex.lab = 1.2)
lines(lag_4_rain_others,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Others\n(Lag = 4 weeks)",cex=.8)
dev.off()
