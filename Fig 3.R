#### Double lollipop plot ######

# Load packages
library(tidyverse)
library(ggtext)
library(ragg)
library(cowplot)

# Set ggplot theme
theme_set(theme_minimal(base_size = 13))
theme_update(
  plot.background = element_rect(color = "#FFFCFC", fill = "#FFFCFC"),
  panel.grid.major.x = element_line(color = "grey94"),
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(color = "grey40"),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  legend.position = c(.07, .31), 
  legend.title = element_text(
    color = "grey40", 
    angle = 90, 
    hjust = .5
  ),
  legend.text = element_text(
    color = "grey40", 
    size = 12
  ),
  legend.box = "horizontal",
  legend.box.just = "bottom",
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing = unit(.6, "lines"),
  plot.title = element_text(
    face = "bold", 
    size = 17.45
  ),
  plot.subtitle = element_textbox_simple(
    color = "grey40", 
    size = 10.8,
    lineheight = 1.3, 
    margin = margin(t = 5, b = 30)
  ),
  plot.caption = element_text(
    color = "grey55", 
    size = 10.5, 
    margin = margin(t = 20, b = 0, r = 15)
  )
)


####### load in data #######

# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

# Adding full country names
waterfetch$country_full <- waterfetch$dhscc
waterfetch$country_full <- gsub("AO", "Angola", waterfetch$country_full)
waterfetch$country_full <- gsub("BF", "Burkina Faso", waterfetch$country_full)
waterfetch$country_full <- gsub("BJ", "Benin", waterfetch$country_full)
waterfetch$country_full <- gsub("BU", "Burundi", waterfetch$country_full)
waterfetch$country_full <- gsub("CD", "Dem. Rep. Congo", waterfetch$country_full)
waterfetch$country_full <- gsub("CF", "Central African Rep.", waterfetch$country_full)
waterfetch$country_full <- gsub("CI", "Cote d'Ivoire", waterfetch$country_full)
waterfetch$country_full <- gsub("CM", "Cameroon", waterfetch$country_full)
waterfetch$country_full <- gsub("GA", "Gabon", waterfetch$country_full)
waterfetch$country_full <- gsub("GH", "Ghana", waterfetch$country_full)
waterfetch$country_full <- gsub("GN", "Guinea", waterfetch$country_full)
waterfetch$country_full <- gsub("KE", "Kenya", waterfetch$country_full)
waterfetch$country_full <- gsub("KM", "Comoros", waterfetch$country_full)
waterfetch$country_full <- gsub("LB", "Liberia", waterfetch$country_full)
waterfetch$country_full <- gsub("LS", "Lesotho", waterfetch$country_full)
waterfetch$country_full <- gsub("MD", "Madagascar", waterfetch$country_full)
waterfetch$country_full <- gsub("ML", "Mali", waterfetch$country_full)
waterfetch$country_full <- gsub("MW", "Malawi", waterfetch$country_full)
waterfetch$country_full <- gsub("MZ", "Mozambique", waterfetch$country_full)
waterfetch$country_full <- gsub("NG", "Nigeria", waterfetch$country_full)
waterfetch$country_full <- gsub("NI", "Niger", waterfetch$country_full)
waterfetch$country_full <- gsub("NM", "Namibia", waterfetch$country_full)
waterfetch$country_full <- gsub("RW", "Rwanda", waterfetch$country_full)
waterfetch$country_full <- gsub("SL", "Sierra Leone", waterfetch$country_full)
waterfetch$country_full <- gsub("SN", "Senegal", waterfetch$country_full)
waterfetch$country_full <- gsub("TD", "Chad", waterfetch$country_full)
waterfetch$country_full <- gsub("TG", "Togo", waterfetch$country_full)
waterfetch$country_full <- gsub("TZ", "Tanzania", waterfetch$country_full)
waterfetch$country_full <- gsub("UG", "Uganda", waterfetch$country_full)
waterfetch$country_full <- gsub("ZM", "Zambia", waterfetch$country_full)
waterfetch$country_full <- gsub("ZW", "Zimbabwe", waterfetch$country_full)

##### Prep data #####
df_ru <- waterfetch %>%
  select(country_full,hv025,hv204) %>%
  group_by(country_full,hv025) %>%
  summarize(median_wt = median(hv204,na.rm=T)) %>%
  mutate(ru=ifelse(hv025==1,"urban","rural")) %>%
  pivot_wider(id_cols=country_full,values_from=median_wt,names_from=ru) %>%
  mutate(diff=rural-urban) %>%
  arrange(urban)

# in the `geom_*` calls when they are not explicitly specified.
p1 <- df_ru %>% 
  ggplot(aes(urban, country_full)) +
  # Segment when shortcut==no. Overlapped lineranges.
  geom_linerange(data = df_ru, aes(xmin = urban, xmax = rural), size = 2)+
  geom_point(aes(x = rural), size = 7, color = "#FFFCFC", fill = "grey45", shape = 21, stroke = .7) +
  # Point when shortcut==yes – latest record. 
  geom_point(fill = 'red', size = 7, color = "#FFFCFC", shape = 21, stroke = .7) +
  # Extend horizontal axis so trackl labels fit
  coord_cartesian(xlim = c(-5, 30)) +
  scale_x_continuous(
    breaks = seq(0, 30, by = 5), 
  ) +
  scale_y_discrete(limits=rev(df_ru$country_full),
                   expand = c(.02, .07)) + 
  geom_text(data = filter(df_ru, country_full == "AO"),
            aes(label = "Walk times\nurban"), 
            size = 3.5, color = "red", 
            lineheight = .8, vjust = 0, nudge_y = .8) +
  geom_text(data = filter(df_ru, country_full == "AO"),
            aes(x = rural, label = "Walk times\nrural"), 
            size = 3.5, color = "grey45", 
            lineheight = .8, vjust = 0, nudge_y = .8) +
  theme_minimal()+
  theme(axis.text.y=element_text(color='black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "grey40"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Walk Times by Rural/Urban Status")

df_elec <- waterfetch %>%
  select(country_full,new_hv206,hv204) %>%
  group_by(country_full,new_hv206) %>%
  summarize(median_wt = median(hv204,na.rm=T)) %>%
  mutate(elec=ifelse(new_hv206==1,"electricity","no_electricity")) %>%
  pivot_wider(id_cols=country_full,values_from=median_wt,names_from=elec) %>%
  mutate(diff=no_electricity-electricity) %>%
  arrange(electricity)

p2 <- df_elec %>% 
  ggplot(aes(electricity, country_full)) +
  # Segment when shortcut==no. Overlapped lineranges.
  geom_linerange(data = df_elec, aes(xmin = electricity, xmax = no_electricity), size = 2)+
  geom_point(aes(x = no_electricity), size = 7, color = "#FFFCFC", fill = "grey65", shape = 21, stroke = .7) +
  geom_point(fill = 'blue', size = 7, color = "#FFFCFC", shape = 21, stroke = .7) +
  # Extend horizontal axis so trackl labels fit
  coord_cartesian(xlim = c(-5, 30)) +
  scale_x_continuous(
    breaks = seq(0, 30, by = 5),
    expand=c(.07,.07)
  ) +
  scale_y_discrete(limits=rev(df_ru$country_full),
                   expand = c(.02, .07)) + 
  geom_text(data = filter(df_elec, country_full == "AO"),
            aes(label = "Walk times\nelectricity"), 
             size = 3.5, color = "#4a5a7b", 
            lineheight = .8, vjust = 0, nudge_y = .8) +
  geom_text(data = filter(df_elec, country_full == "AO"),
            aes(x = no_electricity, label = "Walk times\nno electricity"), 
            size = 3.5, color = "grey50", 
            lineheight = .8, vjust = 0, nudge_y = .8) +
  theme_minimal()+
  theme(axis.text.y=element_text(color='black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "grey40"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Walk Times by Electricity Access")
plot_grid(p1,p2,nrow=1)

ggsave('revision/ru_elec_lollipop.pdf')


