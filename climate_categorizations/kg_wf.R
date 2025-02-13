# This script categorizes data according to Koppen-Geiger Climate Classifications


#### Pulling in climate info
library(raster)

kg <- raster("climate/Beck_KG_V1_present_0p5.tif")

# need to pull in waterfetch data

# Load in required packages
library('data.table')
library('sandwich')
library('lmtest')
library('dplyr')

# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

#### Combining data
waterfetch2 <- waterfetch
coordinates(waterfetch2) <- ~longnum + latnum

waterfetch$kg <- raster::extract(kg,waterfetch2)

tropical <- c(1:3)
arid <- c(4:7)
temperate <- c(8:15)

waterfetch$climate_cat <- factor(ifelse(waterfetch$kg%in%tropical,'TR',
                                 ifelse(waterfetch$kg%in%arid,'AR',
                                        ifelse(waterfetch$kg%in%temperate,'TE',NA))))

saveRDS(waterfetch,'climate/waterfetch_climate.RDS')

# Making a map of the classifications

library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
library(dplyr)

world <- ne_countries(scale="medium",returnclass="sf")

iso2_africa <- c("AO", "BF", "BI", "BJ", "BW", "CD", "CF", "CG", "CI", "CM", "CV", "DJ", "DZ",
                 "EG", "EH","ER", "ET", "GA", "GH", "GM", "GN", "GQ", "GW", "KE", "KM", "LR",
                 "LS", "LY", "MA", "MG", "ML", "MR", "MU", "MW", "MZ", "NA", "NE", "NG", "RE",
                 "RW", "SC", "SD", "SH", "SL", "SN","SO","SS", "ST", "SZ", "TD", "TG", "TN", "TZ",
                 "UG", "YT", "ZA", "ZM", "ZW")
world <- filter(world,iso_a2%in%iso2_africa) #filter to only africa


climatecat_plot <- ggplot(world) +
  geom_point(data=waterfetch,aes(x=longnum, y=latnum,color=climate_cat))+
  geom_sf(fill=NA) +
  coord_sf(xlim=c(-20,53),ylim=c(-35,40),expand=F) +
  theme(legend.position="bottom",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

