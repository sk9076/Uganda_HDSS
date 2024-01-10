### Age stratified (using ERA daily data)

uganda_weekly <- rio::import(here::here("data", "clean data", "dat_mal_weekly_str.xlsx"))


plot_by_age <- ggpubr::ggarrange(
uganda_weekly_str %>% ggplot(aes(x=week_date, y=n.x, fill = age_grp) ) +
  geom_bar(position = "fill", stat = "identity") + xlab("") + ylab("Percentage (%)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                     labels = paste0(seq(0, 100, by = 10), " %")) +
  scale_fill_discrete(name = "Age group") +
  theme(legend.position = "top",
        panel.grid.major = element_blank())
,
uganda_weekly_str %>% ggplot(aes(x=week_date, y=n.x, fill = age_grp) ) +
  geom_bar(stat = "identity")+
  #geom_point(size = .5) + geom_path() +
  ylab("Malaria incidence") + xlab("Date")+
  facet_wrap(.~age_grp,     ncol =1) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "bold"))
,
ncol = 1
)
ggsave(plot = plot_by_age,
       filename = here::here("results", "plot_by_age.png"),
       width = 6,
       height = 6,
       dpi = 200)


