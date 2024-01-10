#### ERA-5 hourly dataset

## Load data
dat_comb <- rio::import(here::here("data", "clean data", "dat_comb_era1.xlsx"))

# Plot to see if the joining worked
dat_comb2 <- dat_comb %<>% mutate(week_date = as.Date(week_date))
dat_comb2 <- dat_comb %>% arrange(week_date, desc=F)
dat_comb2$time = 1:nrow(dat_comb)
dat_comb2 %<>% mutate(n = ifelse(is.na(n), 0, n),
                      week_date = as.Date(week_date),
                      plot_label = format(week_date, "%b\n%Y"))

#dat_comb2 %<>% mutate(plot_label = ifelse(time %% 26==0, plot_label, ""))

#times <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(time) %>% unlist()
#xlabels <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(plot_label) %>% unlist()

times <- seq(range(dat_comb2$time)[1], range(dat_comb2$time)[2], by = 24)
xlabels <- dat_comb2$plot_label[which(dat_comb2$time %in% times)]

png(file = here::here("results", "p_descriptive.png"),
    width = 8,
    height = 4.5,
    res = 200,
    unit = "in")
with(dat_comb2 %>% arrange(by=time),
     {
       #plot(week_date, n, ylim = c(0, 2000))
       par(mar = c(5,5,3,12))
       barplot(n ~ time,
               names.arg = plot_label, cex.names = .5,
               cex.axis = .8,
               col = "#cfb04c", border = "#cfb04c",
               axes = F,
               #names.arg = "",
             ylim = c(0, 2000),
             xlab = "Week",
             ylab = "Weekly malaria incidence")
       #axis(side = 1, at = times, labels = xlabels, cex.axis = .5)
       axis(side = 2)
       par(new=T)
       plot(time, temp, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            ylim = c(5,35),
            col = "red")

       #axis.Date(side = 1, at = seq(min(week_date), max(week_date), by = "6 months"),
       #           format = "%b\n%Y", cex.axis = .5)
       axis(side = 4, at = pretty(range(temp)), cex.axis = .5
       )
       mtext("Average temperature (°C)", side=4, line = 2.5, at=23, cex = .8,
             col="red", las =1) # Rotated y axis label
       par(new=T)
       plot(time, rain, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            col = "blue",
            ylim = c(-0.6*1000, 0.15*1000))
       axis(side = 4, at = pretty(range(rain)), cex.axis = .5)
       mtext("Total rainfall (mm)", side=4, line = 2.5, at=0.06*1000, cex = .8,
             col="blue", las =1) # Rotated y axis label
     })
dev.off()




## MODIS daily

dat_comb <- rio::import(here::here("data", "clean data", "dat_comb_modis.xlsx"))

# Plot to see if the joining worked
dat_comb2 <- dat_comb %<>% mutate(week_date = as.Date(week_date))
dat_comb2 <- dat_comb %>% arrange(week_date, desc=F)
dat_comb2$time = 1:nrow(dat_comb)
dat_comb2 %<>% mutate(n = ifelse(is.na(n), 0, n),
                      week_date = as.Date(week_date),
                      plot_label = format(week_date, "%b\n%Y"))

#dat_comb2 %<>% mutate(plot_label = ifelse(time %% 26==0, plot_label, ""))

#times <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(time) %>% unlist()
#xlabels <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(plot_label) %>% unlist()

times <- seq(range(dat_comb2$time)[1], range(dat_comb2$time)[2], by = 24)
xlabels <- dat_comb2$plot_label[which(dat_comb2$time %in% times)]

png(file = here::here("results", "p_descriptive_modis.png"),
    width = 8,
    height = 4.5,
    res = 200,
    unit = "in")
with(dat_comb2 %>% arrange(by=time),
     {
       #plot(week_date, n, ylim = c(0, 2000))
       par(mar = c(5,5,3,12))
       barplot(n ~ time,
               names.arg = plot_label, cex.names = .5,
               cex.axis = .8,
               col = "#cfb04c", border = "#cfb04c",
               axes = F,
               #names.arg = "",
               ylim = c(0, 2000),
               xlab = "Week",
               ylab = "Weekly malaria incidence")
       #axis(side = 1, at = times, labels = xlabels, cex.axis = .5)
       axis(side = 2)
       par(new=T)
       plot(time, temp, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            ylim = c(0,60),
            col = "red")

       #axis.Date(side = 1, at = seq(min(week_date), max(week_date), by = "6 months"),
       #           format = "%b\n%Y", cex.axis = .5)
       axis(side = 4, at = pretty(range(temp, na.rm=T)), cex.axis = .5
       )
       mtext("Average temperature (°C)", side=4, line = 2.5, at=30, cex = .8,
             col="red", las =1) # Rotated y axis label
       par(new=T)
       plot(time, rain, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            col = "blue",
            ylim = c(-0.6*1000, 0.15*1000))
       axis(side = 4, at = pretty(range(rain)), cex.axis = .5)
       mtext("Total rainfall (mm)", side=4, line = 2.5, at=0.06*1000, cex = .8,
             col="blue", las =1) # Rotated y axis label
     })
dev.off()


#### ERA5 Land Daily

dat_comb <- rio::import(here::here("data", "clean data", "dat_comb_era2.xlsx"))

# Plot to see if the joining worked
dat_comb2 <- dat_comb %<>% mutate(week_date = as.Date(week_date))
dat_comb2 <- dat_comb %>% arrange(week_date, desc=F)
dat_comb2$time = 1:nrow(dat_comb)
dat_comb2 %<>% mutate(n = ifelse(is.na(n), 0, n),
                      week_date = as.Date(week_date),
                      plot_label = format(week_date, "%b\n%Y"))

#dat_comb2 %<>% mutate(plot_label = ifelse(time %% 26==0, plot_label, ""))

#times <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(time) %>% unlist()
#xlabels <- dat_comb2 %>% filter(plot_label!="") %>% dplyr::select(plot_label) %>% unlist()

times <- seq(range(dat_comb2$time)[1], range(dat_comb2$time)[2], by = 24)
xlabels <- dat_comb2$plot_label[which(dat_comb2$time %in% times)]

png(file = here::here("results", "p_descriptive_era2.png"),
    width = 8,
    height = 4.5,
    res = 200,
    unit = "in")
with(dat_comb2 %>% arrange(by=time),
     {
       #plot(week_date, n, ylim = c(0, 2000))
       par(mar = c(5,5,3,12))
       barplot(n ~ time,
               names.arg = plot_label, cex.names = .5,
               cex.axis = .8,
               col = "#cfb04c", border = "#cfb04c",
               axes = F,
               #names.arg = "",
               ylim = c(0, 2000),
               xlab = "Week",
               ylab = "Weekly malaria incidence")
       #axis(side = 1, at = times, labels = xlabels, cex.axis = .5)
       axis(side = 2)
       par(new=T)
       plot(time, temp, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            ylim = c(-15,50),
            col = "red")
       par(new=T)

       plot(time, tmin, type = "l", lty=2,
            axes = F, bty = "n", xlab = "", ylab = "",
            ylim = c(-15,50),
            col = "red")
       par(new=T)

       plot(time, tmax, type = "l", lty =2,
            axes = F, bty = "n", xlab = "", ylab = "",
            ylim = c(-15,50),
            col = "red")

       #axis.Date(side = 1, at = seq(min(week_date), max(week_date), by = "6 months"),
       #           format = "%b\n%Y", cex.axis = .5)
       axis(side = 4, at = pretty(c(min(tmin, na.rm=T), max(tmax, na.rm=T))), cex.axis = .5
       )
       mtext("Average (solid), minimum, \nand maximum (dashed) \ntemperature (°C)", side=4, line = 2.5, at=23, cex = .8,
             col="red", las =1) # Rotated y axis label
       par(new=T)
       plot(time, rain, type = "l",
            axes = F, bty = "n", xlab = "", ylab = "",
            col = "blue",
            ylim = c(-1.2*1000, 0.35*1000))
       axis(side = 4, at = pretty(range(rain)), cex.axis = .5)
       mtext("Total rainfall (mm)", side=4, line = 2.5, at=0.06*1000, cex = .8,
             col="blue", las =1) # Rotated y axis label
     })
dev.off()



###### Cumulative plot (using ERA daily merged data)
#dat_comb <- rio::import(here::here("data", "Uganda", "clean data", "dat_comb_era2.xlsx"))

str(dat_comb)
# aggregate by week
dat_cum <- dat_comb %>% mutate_at(c("week", "year"), as.numeric) %>%
  group_by(week) %>%
  mutate(n_cum = mean(n, na.rm=T),
         temp_cum = mean(temp, na.rm=T),
         tmin_cum = mean(tmin, na.rm=T),
         tmax_cum = mean(tmax, na.rm=T),
         rain_cum = mean(rain, na.rm=T),
         week_date = as.Date(paste0(2022, "-", week, "-1"), format = "%Y-%U-%u")
         )

ggthemr::ggthemr("fresh")
plot_cum <- ggpubr::ggarrange(
ggplot(dat_cum) +
  geom_ribbon(aes(x = week_date, ymin = tmin_cum, ymax = tmax_cum), fill = "red", alpha = 0.2)+
  geom_line(aes(week_date, temp_cum), color = "red") +
  scale_x_date(date_breaks = "1 month", date_labels="%b")+
  xlab("") + ylab("Temperature (°C)")
,
ggplot(dat_cum) +
  geom_area(aes(week_date, rain_cum), fill = "blue", alpha = 0.5)+
  scale_x_date(date_breaks = "1 month", date_labels="%b")+
  xlab("") + ylab("Rainfall (mm)")
,
ggplot(dat_cum) +
  geom_bar(aes(week_date, n_cum), stat = "identity", fill = "#cfb04c") +
  scale_x_date(date_breaks = "1 month", date_labels="%b")+
  xlab("Time within a given year") + ylab("Malaria incidence"),
ncol = 1
)

ggsave(here::here("results", "plot_cum.png"),
       plot_cum,
       width = 5,
       height = 7,
       dpi = 200)
