#### Visualization based on the final built model

#model_rain_tmax
#model_rain_tmin
#model_rain

#cb_rain
#cb_temp
#cb_tmax
#cb_tmin

bound_rain <- range(dat_comb$rain, na.rm=T)
bound_temp <- range(dat_comb$temp, na.rm=T)
bound_tmin <- range(dat_comb$tmin, na.rm=T)
bound_tmax <- range(dat_comb$tmax, na.rm=T)

cen_rain <- bound_rain[1] # or 21? 80?
cen_temp <- bound_temp[1]
cen_tmax <- bound_tmax[1]
cen_tmin <- 20


#### Final model 1 : Rain + Tmax

model <- model_rain_tmax

# Association between the rainfall and malaria incidence
cp_rain <- crosspred(cb_rain,
                     model,
                     from=bound_rain[1],
                     to=bound_rain[2], by=1,
                     cen=cen_rain)

cr_rain <- crossreduce(cb_rain,
                       model,
                       from=bound_rain[1],
                       to=bound_rain[2], by=1,
                       cen=cen_rain)

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


## Overall contour plots side-by-side
dev.off()
png(here::here("results", "rain_tmax", "rain_tmax_3d.png"),
    height = 8,
    width = 16,
    units = "in",
    res = 200)

par(mfrow=c(1,2))
mar = c(1,1,1,1)
d3<- plot(cp_rain,
          xlab="Weekly Total Rainfall (mm)",cex=2,
          zlab="Relative Risk",cex=2,
          phi=30,theta=325,ltheta=110, # plot rotation
          shade=.2,
          lty = 0,
          #lwd = .2,
          col="#aad5f2")

lines(trans3d(x=cp_rain$predvar,
              y=6,
              z=cp_rain$matRRfit[,"lag0"],
              pmat=d3),lwd=2,col=2)

lines(trans3d(x=cp_rain$predvar,
              y=2,
              z=cp_rain$matRRfit[,"lag2"],
              pmat=d3),lwd=2,col=2)

lines(trans3d(x=cp_rain$predvar,
              y=4,
              z=cp_rain$matRRfit[,"lag4"],
              pmat=d3),lwd=2,col=2)


d3<- plot(cp_tmax,
          xlab="Maximum Temperature (°C)",cex=2,
          zlab="Relative Risk",cex=2,
          phi=30,theta=325,ltheta=110, # plot rotation
          shade=.2,
          lty = 0,
          #lwd = .2,
          col="#708ccc")
lines(trans3d(x=cp_tmax$predvar,
              y=0,
              z=cp_tmax$matRRfit[,"lag0"],
              pmat=d3),lwd=2,col=2)

dev.off()

png(here::here("results", "rain_tmax", "rain_tmax_2d_rain.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)
mar = c(3, 3, 3,3)
plot(cp_rain,"contour",xlab="Rainfall (mm)", key.title=title("RR"),
     ylab="lag (Weeks)")

dev.off()

png(here::here("results", "rain_tmax", "rain_tmax_2d_tmax.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)
plot(cp_tmax,"contour",xlab="Maximum Temperature (°C)", key.title=title("RR"),
     ylab="lag (Weeks)")
dev.off()


png(here::here("results", "rain_tmax", "rain_tmax_cum.png"),
    height = 6,
    width = 5,
    units = "in",
    res = 200)

par(mfrow = c(2,1),
    mar = c(4, 4, 3,1))

plot(cr_rain,type="n",ci="n",
     ylab="Relataive Risk (RR)",
     ylim=c(0,50),
     xlab="Weekly Rainfall (mm)")
lines(cr_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Overall Cumulative Association",line = 1.2,cex=1)

mar = c(4, 4, 1, 1)
plot(cr_tmax,type="n",ci="n",
     ylab="Relataive Risk (RR)",
     ylim=c(0.5,5),
     xlab="Maximum Temperature (°C)")
lines(cr_rain,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
#mtext(text="Overall Cumulative Association",line = 1.2,cex=1)
dev.off()


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


png(here::here("results", "rain_tmax", "rain_tmax_lag_rain.png"),
    height = 6,
    width = 10,
    units = "in",
    res = 200)
par(mfrow=c(2,3))
mar=c(0,0,0,0)

plot(lag_0_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_0_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 0 Weeks",cex=1)


plot(lag_2_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_2_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 2 Weeks",cex=1)


plot(lag_4_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_4_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 4 Weeks",cex=1)


plot(lag_6_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_6_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 6 Weeks",cex=1)


plot(lag_8_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_8_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 8 Weeks",cex=1)

plot(lag_12_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2),xlab=expression("Rainfall (mm)"))
lines(lag_12_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 12 Weeks",cex=1)
dev.off()

## Association between the temperature and the malaria incidence


# Lag-specific: Tmax

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

png(here::here("results", "rain_tmax", "rain_tmax_lag_tmax.png"),
    height = 6,
    width = 10,
    units = "in",
    res = 200)
par(mfrow=c(2,3))
mar=c(0,0,0,0)

plot(lag_0_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_0_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 0 Weeks",cex=1)


plot(lag_2_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_2_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 2 Weeks",cex=1)


plot(lag_4_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_4_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 4 Weeks",cex=1)


plot(lag_6_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_6_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 6 Weeks",cex=1)


plot(lag_8_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_8_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 8 Weeks",cex=1)

plot(lag_12_tmax,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0,3),xlab=expression("Maximum Temperature ("*~degree*C*")"))
lines(lag_12_tmax,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 12 Weeks",cex=1)
dev.off()






########## Final model 2 : Rain + Tmin

model <- model_rain_tmin
cen_tmin <- 20
# Association between the rainfall and malaria incidence
cp_rain <- crosspred(cb_rain,
                     model,
                     from=bound_rain[1],
                     to=bound_rain[2], by=1,
                     cen=cen_rain)

cr_rain <- crossreduce(cb_rain,
                       model,
                       from=bound_rain[1],
                       to=bound_rain[2], by=1,
                       cen=cen_rain)

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


## Overall contour plots side-by-side

png(here::here("results", "rain_tmin", "rain_tmin_3d.png"),
    height = 8,
    width = 16,
    units = "in",
    res = 200)

par(mfrow=c(1,2))
mar = c(1,1,1,1)
d3<- plot(cp_rain,
          xlab="Weekly Total Rainfall (mm)",cex=2,
          zlab="Relative Risk",cex=2,
          phi=30,theta=325,ltheta=110, # plot rotation
          shade=.2,
          lty = 0,
          #lwd = .2,
          col="#aad5f2")

lines(trans3d(x=cp_rain$predvar,
              y=0,
              z=cp_rain$matRRfit[,"lag0"],
              pmat=d3),lwd=2,col=2)

lines(trans3d(x=cp_rain$predvar,
              y=2,
              z=cp_rain$matRRfit[,"lag2"],
              pmat=d3),lwd=2,col=2)

lines(trans3d(x=cp_rain$predvar,
              y=4,
              z=cp_rain$matRRfit[,"lag4"],
              pmat=d3),lwd=2,col=2)


d3<- plot(cp_tmin,
          xlab="Minimum Temperature (°C)",cex=2,
          zlab="Relative Risk",cex=2,
          phi=30,theta=325,ltheta=110, # plot rotation
          shade=.2,
          lty = 0,
          #lwd = .2,
          col="#708ccc")
lines(trans3d(x=cp_tmin$predvar,
              y=12,
              z=cp_tmin$matRRfit[,"lag12"],
              pmat=d3),lwd=2,col=2)

dev.off()

png(here::here("results", "rain_tmin", "rain_tmin_2d_rain.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)
mar = c(3, 3, 3,3)
plot(cp_rain,"contour",xlab="Rainfall (mm)", key.title=title("RR"),
     ylab="lag (Weeks)")

dev.off()

png(here::here("results", "rain_tmin", "rain_tmin_2d_tmin.png"),
    height = 6,
    width = 7,
    units = "in",
    res = 200)
plot(cp_tmin,"contour",xlab="Minimum Temperature (°C)", key.title=title("RR"),
     ylab="lag (Weeks)")
dev.off()


png(here::here("results", "rain_tmin", "rain_tmin_cum.png"),
    height = 6,
    width = 5,
    units = "in",
    res = 200)

par(mfrow = c(2,1),
    mar = c(4, 4, 3,1))

plot(cr_rain,type="n",ci="n",
     ylab="Relataive Risk (RR)",
     ylim=c(0,50),
     xlab="Weekly Rainfall (mm)")
lines(cr_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Overall Cumulative Association",line = 1.2,cex=1)

mar = c(4, 4, 1, 1)
plot(cr_tmin,type="n",ci="n",
     ylab="Relataive Risk (RR)",
     ylim=c(0.5,2.5),
     xlab="Minimum Temperature (°C)")
lines(cr_rain,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
#mtext(text="Overall Cumulative Association",line = 1.2,cex=1)
dev.off()


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


png(here::here("results", "rain_tmin", "rain_tmin_lag_rain.png"),
    height = 6,
    width = 10,
    units = "in",
    res = 200)
par(mfrow=c(2,3))
mar=c(0,0,0,0)

plot(lag_0_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_0_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 0 Weeks",cex=1)


plot(lag_2_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_2_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 2 Weeks",cex=1)


plot(lag_4_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_4_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 4 Weeks",cex=1)


plot(lag_6_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_6_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 6 Weeks",cex=1)


plot(lag_8_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_8_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 8 Weeks",cex=1)

plot(lag_12_rain,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(.4,2.5),xlab=expression("Rainfall (mm)"))
lines(lag_12_rain,col="#0a1054",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=0.714, green=0.878, blue=0.949,0.5)))
mtext(text="Lag = 12 Weeks",cex=1)
dev.off()

## Association between the temperature and the malaria incidence


# Lag-specific: Tmax

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

png(here::here("results", "rain_tmin", "rain_tmin_lag_tmin.png"),
    height = 6,
    width = 10,
    units = "in",
    res = 200)
par(mfrow=c(2,3))
mar=c(0,0,0,0)

plot(lag_0_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_0_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 0 Weeks",cex=1)


plot(lag_2_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_2_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 2 Weeks",cex=1)


plot(lag_4_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_4_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 4 Weeks",cex=1)


plot(lag_6_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_6_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 6 Weeks",cex=1)


plot(lag_8_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_8_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 8 Weeks",cex=1)

plot(lag_12_tmin,type="n",ci="n",
     ylab="Relative Risk",
     ylim=c(0.5,2.0),xlab=expression("Minimum Temperature ("*~degree*C*")"))
lines(lag_12_tmin,col="#bf300f",lty=1,lwd=2,ci="area",
      ci.arg=list(col=rgb(red=1, green=0.835, blue=0.631,0.6)))
mtext(text="Lag = 12 Weeks",cex=1)
dev.off()
