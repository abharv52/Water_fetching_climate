##### Creating maps ########
# Abby Harvey
# July 10, 2020

library(tidyverse)
library(sf)
library(cowplot)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale="medium",returnclass="sf")

iso2_africa <- c("AO", "BF", "BI", "BJ", "BW", "CD", "CF", "CG", "CI", "CM", "CV", "DJ", "DZ",
                 "EG", "EH","ER", "ET", "GA", "GH", "GM", "GN", "GQ", "GW", "KE", "KM", "LR",
                 "LS", "LY", "MA", "MG", "ML", "MR", "MU", "MW", "MZ", "NA", "NE", "NG", "RE",
                 "RW", "SC", "SD", "SH", "SL", "SN","SO","SS", "ST", "SZ", "TD", "TG", "TN", "TZ",
                 "UG", "YT", "ZA", "ZM", "ZW")
world <- filter(world,iso_a2%in%iso2_africa) #filter to only africa



# # Loading country border data
# world_map <- read_sf("world-administrative-boundaries/world-administrative-boundaries.shp")
# africa <- filter(world_map,continent=="Africa")
# # ggplot(africa) + geom_sf()

waterfetch <- readRDS('waterfetch.rds')

# Summarizing by cluster
waterfetch2 <- waterfetch %>%
  group_by(clusterid) %>%
  summarise(walk=mean(hv204,na.rm=T),
            tmax7=mean(tmax6d7dl,na.rm=T),
            precip7=mean(precip6d7dl,na.rm=T),
            latnum=mean(latnum),
            longnum=mean(longnum))

wflow <- filter(waterfetch2,walk<80)
wfhigh <- filter(waterfetch2,walk>=80)


# Map of all data points
worldmap <- ggplot(world) +
  # geom_sf(color=NA,fill=NA)+
  geom_point(data=waterfetch2,aes(x=longnum,y=latnum),size=0.5,alpha=0.5)+
  geom_sf(data=world,fill=NA,color="black")+
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
g1 <- ggplotGrob(worldmap)
ggsave('Paper figures/all_map.pdf',g1,width=6,height=6,useDingbats=F)

###### Weather data plots #####
waterfetch2$precip7 <- ifelse(waterfetch2$precip7>60,60,waterfetch2$precip7)

precip_plot <- ggplot(world)+
  # geom_sf(color=NA,fill=NA) +
  geom_point(data=waterfetch2,aes(x=longnum, y=latnum,color=precip7),size=1.05,shape=15) +
  geom_sf(data=world,fill=NA,color="black")+
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlim(-20, 53) + ylim(-35, 40)+
  scale_color_gradient(low = "#A1F8A3",  high = "#066507",breaks=c(20,40,60),labels=c(20,40,60)) +
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 title.position="top",
                                 title="Mean Weekly\nPrecipitation (cm)")) +
  ggtitle('Precipitation')


tmax_plot <- ggplot(world)+
  geom_point(data=waterfetch2,aes(x=longnum, y=latnum,color=tmax7),size=1.05,shape=15) +
  geom_sf(data=world,fill=NA,color="black")+
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlim(-20, 53) + ylim(-35, 40)+
  scale_color_gradient(low = "#FFFC00",  high = "#FF1600",breaks=c(15,25,35),labels=c(15,25,35)) +
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 title.position="top",
                                 title="Mean Daily Maximum\nTemperature")) +
  ggtitle("Daily Maximum Temperature")

# Walk time plot
wt_df <- waterfetch2 %>%
  filter(!is.na(walk)) %>%
  mutate(walk=ifelse(walk>100,100,walk))

wt_plot <- ggplot(world)+
  geom_point(data=wt_df,aes(x=longnum, y=latnum,color=walk),size=1.05,shape=15) +
  geom_sf(data=world,fill=NA,color="black")+
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlim(-20, 53) + ylim(-35, 40)+
  scale_color_gradient(low = "#0AFBFE",  high = "#0D5CC7")+#,breaks=c(20,40,60),labels=c(20,40,60)) +
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 title.position="top",
                                 title="Mean Walk\nTime (min)")) +
  ggtitle("One Way Walk Times")


### Climate Zone ####
wf_climate <- readRDS('waterfetch_climate.RDS')

wf_climate <- wf_climate %>% group_by(clusterid) %>% sample_n(.,1) %>%
  filter(!is.na(climate_cat)) %>%
  mutate(climate_cat=factor(climate_cat,labels=c("Arid","Temperature","Tropical"),levels=c("AR","TE","TR")))

climate_plot <- ggplot(world)+
  geom_point(data=wf_climate,aes(x=longnum, y=latnum,color=climate_cat),size=1.05,shape=15,alpha=0.2) +
  geom_sf(data=world,fill=NA,color="black")+
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlim(-20, 53) + ylim(-35, 40)+
  # scale_color_gradient(low = "#FFFC00",  high = "#FF1600",breaks=c(0.5,1,1.5,2),labels=c(0.5,1,1.5,2)) +
  guides(colour = guide_legend(title.position="top",
                                 title="Climate Zone",
                               override.aes=list(size=5,alpha=1))) +
  ggtitle("Climate Zones")
