# setwd("/Users/abbyharvey/Documents/GitHub/water-fetch/code/active_r_code")

# Load in required packages
library('data.table')
library('sandwich')
library('lmtest')
library('dplyr')
library(cowplot)


# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

improved_sources <- c(1:5,25:29,7:11,19)
# add in water source binary variables
waterfetch <- waterfetch %>%
  mutate(groundwater=ifelse(water_source%in%c(1:5,20:29),1,ifelse(is.na(water_source),NA,0)),
         boreholes=ifelse(water_source%in%c(1:5),1,ifelse(is.na(water_source),NA,0)),
         wells=ifelse(water_source%in%c(20:29),1,ifelse(is.na(water_source),NA,0)),
         surface_water=ifelse(water_source==18,1,ifelse(is.na(water_source),NA,0)),
         vendor_water=ifelse(water_source%in%c(6,12:15),1,ifelse(is.na(water_source),NA,0)),
         unimproved_spring=ifelse(water_source==16,1,ifelse(is.na(water_source),NA,0)),
         rain_water = ifelse(water_source==19,1,ifelse(is.na(water_source),NA,0)),
         piped_water = ifelse(water_source%in%c(7:11),1,ifelse(is.na(water_source),NA,0)),
         protected_spring = ifelse(water_source==17,1,ifelse(is.na(water_source),NA,0)),
         other_source=ifelse(water_source==30,1,ifelse(is.na(water_source),NA,0)),
         missing=ifelse(is.na(water_source),1,0),
         improved=ifelse(water_source%in%improved_sources,1,ifelse(is.na(water_source),NA,0)))

lags <- c('6d7dl','29d30dl','89d90dl','179d180dl','364d365dl')

########### Scaling precip variable ####################
allprecip <- paste('precip',lags,sep="")
precipLengths <- c(7,30,90,180,365)
precipWeeks <- precipLengths/7

for (var in 1:length(allprecip)){
  waterfetch[,allprecip[var]] <- waterfetch[,allprecip[var]]/precipWeeks[var]
}

outcomes <- c("piped_water","surface_water","boreholes","wells","improved")
climatevars <- c(paste0('tmax',lags),paste0('precip',lags))

# make key variables factors
waterfetch <- waterfetch %>%
  mutate(year=factor(year),
         survey=factor(dhscc),
         month=factor(month))

analysis_vars <- expand.grid(outcomes,climatevars)
analysis_vars$Var1 <- as.character(analysis_vars$Var1)
analysis_vars$Var2 <- as.character(analysis_vars$Var2)

urban <- filter(waterfetch,hv025==1)
rural <- filter(waterfetch,hv025==0)

model_results <- data.frame(matrix(ncol=5,nrow=0))

for (i in 1:nrow(analysis_vars)){
  modeldf <- waterfetch[,c(analysis_vars[i,1],analysis_vars[i,2],"year","survey","month","clusterid")]
  names(modeldf)[1:2] <- c("var1","var2")
  
  m1 <- glm(var1 ~ var2, family=binomial,data = modeldf)
  m1_cl <- coeftest(m1, vcov = vcovCL, cluster = modeldf[c("year","survey","month")])
  
  model_out <- c(analysis_vars[i,1],analysis_vars[i,2],
                 m1$coefficients[2],m1_cl[2,4],m1_cl[2,2]*1.96)
  model_results <- rbind(model_results,model_out)
}
names(model_results) <- c("var1","var2",
                          "estimate","p","ci")
# convert to odds ratios
model_results <- model_results %>%
  mutate(across(c(estimate,ci),\(x) (as.numeric(x))))
model_results$ci[model_results$estimate<0] <- -1*model_results$ci[model_results$estimate<0]
model_results$lower_unadj <- model_results$estimate-model_results$ci
model_results$upper_unadj <- model_results$estimate+model_results$ci

model_results[,c(3,5:7)] <- exp(model_results[,c(3,5:7)])

model_results$combo <- paste0(model_results$var1,model_results$var2)

##### Rural #####
model_rural <- data.frame(matrix(ncol=5,nrow=0))

for (i in 1:nrow(analysis_vars)){
  modeldf <- rural[,c(analysis_vars[i,1],analysis_vars[i,2],"year","survey","month","clusterid")]
  names(modeldf)[1:2] <- c("var1","var2")
  
  m1 <- glm(var1 ~ var2, family=binomial,data = modeldf)
  m1_cl <- coeftest(m1, vcov = vcovCL, cluster = modeldf[c("year","survey","month")])
  
  model_out <- c(analysis_vars[i,1],analysis_vars[i,2],
                 m1$coefficients[2],m1_cl[2,4],m1_cl[2,2]*1.96)
  model_rural <- rbind(model_rural,model_out)
}
names(model_rural) <- c("var1","var2",
                          "estimate","p","ci")
# convert to odds ratios
model_rural <- model_rural %>%
  mutate(across(c(estimate,ci),\(x) (as.numeric(x))))
model_rural$ci[model_rural$estimate<0] <- -1*model_rural$ci[model_rural$estimate<0]
model_rural$lower_unadj <- model_rural$estimate-model_rural$ci
model_rural$upper_unadj <- model_rural$estimate+model_rural$ci

model_rural[,c(3,5:7)] <- exp(model_rural[,c(3,5:7)])

model_rural$combo <- paste0(model_rural$var1,model_rural$var2)

precip_rural <- filter(model_rural,var2%in%c(paste0('precip',lags)))
precip_rural$var2 <- factor(precip_rural$var2,levels=c("precip364d365dl","precip179d180dl","precip89d90dl",
                                                       "precip29d30dl","precip6d7dl"))
precip_rural$var1 <- factor(precip_rural$var1,levels=rev(c("surface_water","wells","boreholes","piped_water","improved")))

p1_rural <- ggplot(precip_rural,aes(x=estimate,y=var2,color=var2)) + geom_point() +
  facet_grid(rows=vars(var1))+
  theme_bw() + geom_vline(xintercept=1,linetype="dashed") + xlim(0.95,1.06) +
  geom_errorbar(aes(xmin=lower_unadj,xmax=upper_unadj)) + ggtitle("Precipitation Rural")+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=rev(c("#00B0F6","#A3A500","#E76BF3","#F8766D","#00BF7D")))


temp_rural <- filter(model_rural,var2%in%c(paste0('tmax',lags)))
temp_rural$var2 <- factor(temp_rural$var2,levels=c("tmax364d365dl","tmax179d180dl","tmax89d90dl",
                                                       "tmax29d30dl","tmax6d7dl"))
temp_rural$var1 <- factor(temp_rural$var1,levels=rev(c("surface_water","wells","boreholes","piped_water","improved")))

p2_rural <- ggplot(temp_rural,aes(x=estimate,y=var2,color=var2)) + geom_point() +
  facet_grid(rows=vars(var1))+
  theme_bw() + geom_vline(xintercept=1,linetype="dashed") +xlim(0.85,1.15) +
  geom_errorbar(aes(xmin=lower_unadj,xmax=upper_unadj)) + ggtitle("Temperature Rural")+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=rev(c("#00B0F6","#A3A500","#E76BF3","#F8766D","#00BF7D")))
# plot_grid(p1_rural,p2_rural,ncol=2)

##### Urban ######
model_urban <- data.frame(matrix(ncol=5,nrow=0))

for (i in 1:nrow(analysis_vars)){
  modeldf <- urban[,c(analysis_vars[i,1],analysis_vars[i,2],"year","survey","month","clusterid")]
  names(modeldf)[1:2] <- c("var1","var2")
  
  m1 <- glm(var1 ~ var2, family=binomial,data = modeldf)
  m1_cl <- coeftest(m1, vcov = vcovCL, cluster = modeldf[c("year","survey","month")])
  
  model_out <- c(analysis_vars[i,1],analysis_vars[i,2],
                 m1$coefficients[2],m1_cl[2,4],m1_cl[2,2]*1.96)
  model_urban <- rbind(model_urban,model_out)
}
names(model_urban) <- c("var1","var2",
                        "estimate","p","ci")
# convert to odds ratios
model_urban <- model_urban %>%
  mutate(across(c(estimate,ci),\(x) (as.numeric(x))))
model_urban$ci[model_urban$estimate<0] <- -1*model_urban$ci[model_urban$estimate<0]
model_urban$lower_unadj <- model_urban$estimate-model_urban$ci
model_urban$upper_unadj <- model_urban$estimate+model_urban$ci

model_urban[,c(3,5:7)] <- exp(model_urban[,c(3,5:7)])

model_urban$combo <- paste0(model_urban$var1,model_urban$var2)

precip_urban <- filter(model_urban,var2%in%c(paste0('precip',lags)))
precip_urban$var2 <- factor(precip_urban$var2,levels=c("precip364d365dl","precip179d180dl","precip89d90dl",
                                                       "precip29d30dl","precip6d7dl"))
precip_urban$var1 <- factor(precip_urban$var1,levels=rev(c("surface_water","wells","boreholes","piped_water","improved")))


p1_urban <- ggplot(precip_urban,aes(x=estimate,y=var2,color=var2)) + geom_point() +
  facet_grid(rows=vars(var1))+
  theme_bw() + geom_vline(xintercept=1,linetype="dashed") + xlim(0.95,1.06) +
  geom_errorbar(aes(xmin=lower_unadj,xmax=upper_unadj)) + ggtitle("Precipitation Urban")+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=rev(c("#00B0F6","#A3A500","#E76BF3","#F8766D","#00BF7D")))


temp_urban <- filter(model_urban,var2%in%c(paste0('tmax',lags)))
temp_urban$var2 <- factor(temp_urban$var2,levels=c("tmax364d365dl","tmax179d180dl","tmax89d90dl",
                                                   "tmax29d30dl","tmax6d7dl"))
temp_urban$var1 <- factor(temp_urban$var1,levels=rev(c("surface_water","wells","boreholes","piped_water","improved")))


p2_urban <- ggplot(temp_urban,aes(x=estimate,y=var2,color=var2)) + geom_point() +
  facet_grid(rows=vars(var1))+
  theme_bw() + geom_vline(xintercept=1,linetype="dashed") +xlim(0.85,1.2) +
  geom_errorbar(aes(xmin=lower_unadj,xmax=upper_unadj)) + ggtitle("Temperature Urban")+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=rev(c("#00B0F6","#A3A500","#E76BF3","#F8766D","#00BF7D")))

precip_plot <- plot_grid(p1_urban,p1_rural,ncol=2)
temp_plot <- plot_grid(p2_urban,p2_rural,ncol=2)
 
ggsave('sourceuse_OR_precip.pdf',precip_plot,width=10,height=4,useDingbats=F)
ggsave('sourceuse_OR_tmax.pdf',temp_plot,width=10,height=4,useDingbats=F)
