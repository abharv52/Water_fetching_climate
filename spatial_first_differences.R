#############################
# Code to run all spatial first differences regressions
# Abby Harvey
# Updated Dec 21, 2023


# Overall parameters to change
datasubset <- "all" # enter either urban, rural, elec, noelec, or all
secondsubset <- "all" # Should be different than datasubset, or "all"
reps <- 1000 #If using s1 method, specify number of replicates; otherwise, value does not matter

# This file will subset data as specified,
# run the SFD regressions,
# and output one csv file containing all results

# The filename used to save results will be:
# sfd_datasubset.csv

if (identical(secondsubset,"all")){
  filename <- paste("climatesfd",datasubset,reps,sep="_")
} else {
  filename <- paste("climatesfd",datasubset,secondsubset,"nit",reps,sep="_")  
}
filename <- paste0(filename,".csv")

# Load in required packages
library('data.table')
library('sandwich')
library('lmtest')
library('dplyr')

# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

############## Inputting SFD prep function ##############################
sfd_prep <- function(survey_year,waterfetch,preptype,direction,CPClat,CPClon){
  # Calculates first differences spatially, survey-by-survey
  # Running this function once will only calculate for one survey, not all
  # Need to use lapply to run this on all surveys and compile results
  
  # Subset to given survey
  waterfetch.survey <- subset(waterfetch,survey%in%survey_year)
  
  # Appending lat & lon grid number to data points
  waterfetch.survey$latcell <- cut(waterfetch.survey$latnum,CPClat,labels=FALSE)
  waterfetch.survey$loncell <- cut(waterfetch.survey$longnum,CPClon,labels=FALSE)
  
  # Specifying climate variables to keep
  climatevar <- c('tmax','precip')
  lags <- c('6d7dl','13d7dl','29d30dl','59d30dl','89d90dl','179d180dl','364d365dl')
  allvar.regress <- c("hv204",unlist(lapply(climatevar,paste,lags,sep="")))
  lagsc <- c('6d7dc_mean','13d7dc_mean','29d30dc_mean','59d30dc_mean','89d90dc_mean','364d365dc_mean')
  # Keeping all variables we need
  allvar <- c(allvar.regress,"latcell","loncell","survey","hv025","new_hv206","datetime",unlist(lapply(climatevar,paste,lagsc,sep="")))
  
  
  if (identical(direction,'lon')){
    #sort by lat then lon
    waterfetch.survey <- waterfetch.survey[order(waterfetch.survey$latcell,waterfetch.survey$loncell), ]
    
    # Subsetting to only needed variables
    waterfetch.survey <- waterfetch.survey[,names(waterfetch.survey)%in%allvar]
    
    if (identical(preptype,'mean')){
      # Averaging all climate variables by unique lat-lon combinations
      waterfetch.grids <- waterfetch.survey %>%
        group_by(latcell,loncell) %>%
        summarise_each(funs(mean(., na.rm=TRUE)))
    } else if (identical(preptype,'s1')){
      # Randomly select one data point per unique lat-lon combinations
      waterfetch.grids <- waterfetch.survey %>%
        group_by(latcell,loncell) %>%
        sample_n(.,1)
    } else {print('Invalid type - must be mean or s1')}
    
    #Creating variables to tell of each point's relationship to the preceding point
    waterfetch.grids$latcelldiff <- waterfetch.grids$latcell-lag(waterfetch.grids$latcell,1)
    waterfetch.grids$loncelldiff <- waterfetch.grids$loncell-lag(waterfetch.grids$loncell,1)
    waterfetch.grids$daydiff <- NA
    if (nrow(waterfetch.grids)>2){
      for (i in 2:nrow(waterfetch.grids)){
        if ((!is.na(waterfetch.grids$datetime[i]))&(!is.na(waterfetch.grids$datetime[i-1]))){
          waterfetch.grids$daydiff[i] <- daydiff(waterfetch.grids$datetime[i],waterfetch.grids$datetime[i-1])
        }
      }}
    diff <- waterfetch.grids[,names(waterfetch.grids)%in%allvar.regress]
    diff <- diff-mutate_all(diff,lag)
    colnames(diff) <- paste(colnames(diff),"sfd",sep=".")
    waterfetch.sfd <- data.frame(waterfetch.grids,diff)
    waterfetch.sfd <- filter(waterfetch.sfd,latcelldiff==0&loncelldiff==1)
  } else if (identical(direction,'lat')){
    # Subsetting to only needed variables
    waterfetch.survey <- waterfetch.survey[,names(waterfetch.survey)%in%allvar]
    waterfetch.survey <- waterfetch.survey[order(waterfetch.survey$loncell,waterfetch.survey$latcell), ]
    
    if (identical(preptype,'mean')){
      # Averaging all climate variables by unique lat-lon combinations
      waterfetch.grids <- waterfetch.survey %>%
        group_by(loncell,latcell) %>%
        summarise_each(funs(mean(., na.rm=TRUE)))
    } else if (identical(preptype,'s1')){
      # Randomly select one data point per unique lat-lon combinations
      waterfetch.grids <- waterfetch.survey %>%
        group_by(loncell,latcell) %>%
        sample_n(.,1)
    } else {print('Invalid type - must be mean or s1')}
    
    #Creating variables to tell of each point's relationship to the preceding point
    waterfetch.grids$latcelldiff <- waterfetch.grids$latcell-lag(waterfetch.grids$latcell,1)
    waterfetch.grids$loncelldiff <- waterfetch.grids$loncell-lag(waterfetch.grids$loncell,1)
    waterfetch.grids$daydiff <- NA
    if (nrow(waterfetch.grids)>2){
      for (i in 2:nrow(waterfetch.grids)){
        if ((!is.na(waterfetch.grids$datetime[i]))&(!is.na(waterfetch.grids$datetime[i-1]))){
          waterfetch.grids$daydiff[i] <- daydiff(waterfetch.grids$datetime[i],waterfetch.grids$datetime[i-1])
        }
      }}
    diff <- waterfetch.grids[,names(waterfetch.grids)%in%allvar.regress]
    diff <- diff-mutate_all(diff,lag)
    colnames(diff) <- paste(colnames(diff),"sfd",sep=".")
    waterfetch.sfd <- data.frame(waterfetch.grids,diff)
    waterfetch.sfd <- filter(waterfetch.sfd,latcelldiff==1&loncelldiff==0)
  }
  if (nrow(waterfetch.sfd)){
    waterfetch.sfd$survey <- survey_year
    waterfetch.sfd$dir <- direction}
  return(waterfetch.sfd)
}

############## Inputting SFD batchlm function ###########################
sfd.batchlm <- function(climatevar,dependent,data,timelag,direction,...){
  #This is a function to run batch regressions using lm on SFD data
  if (identical(direction,'lon')){errordir <- 'latcell'
  } else if (identical(direction,'lat')){errordir <- 'loncell'
  } else {print("Invalid direction - must be lat or lon")
    break}
  errors <- c("survey",errordir)
  lmresult <- data.frame()
  
  climate <- paste(climatevar,timelag,".sfd",sep="")
  rframe <- data.frame(dependent,data[,climate],...)#create data frame for regressions
  # Run lm
  model <- lm(dependent~., data=rframe)
  modelvcov <- coeftest(model, vcov=vcovHC(model,type="HC3",cluster=errors))
  # Save output as row in data frame
  rownames(modelvcov)[2] <- climate
  modelsum <- data.frame(row.names=1)
  count=1
  for (term in rownames(modelvcov)){
    temp <- data.frame(term,t(modelvcov[term,]))
    colnames(temp)[1] <- 'VarName'
    colnames(temp) <- paste(colnames(temp),count,sep="")
    if (count==2){
      temp$VarSD2 <- sd(data[,term],na.rm=TRUE)
    }
    modelsum <- cbind(modelsum,temp)
    count=count+1
  }
  temp2 <- cbind(summary(model)$df[2],summary(model)$r.squared,summary(model)$adj.r.squared,summary(model)$sigma,summary(model)$fstatistic[1],mean(abs(data$daydiff)),modelsum)
  colnames(temp2)[1:6] <- c('n','R_squared','Adj_r_squared','sigma','fstatistic','daydiff')
  lmresult <- structure(rbind(lmresult,temp2),.Names=names(temp2))
  lmresult$dir <- direction
  return(lmresult)
  #output as dataframe with stats
}


subsetWF <- function(waterfetch,subsetStr){
  # Function to subset data by urban, rural, elec, or no elec
  if (identical(subsetStr,"urban")){subsetWF <- subset(waterfetch,hv025==1)
  } else if (identical(subsetStr,"rural")){subsetWF <- subset(waterfetch,hv025==0)
  } else if (identical(subsetStr,"elec")){subsetWF <- subset(waterfetch,new_hv206==1)
  } else if (identical(subsetStr,"noelec")){subsetWF <- subset(waterfetch,new_hv206==0)
  } else {subsetWF <- waterfetch}
  return(subsetWF)
}
######## Subsetting to datasubset ####################
waterfetch <- subsetWF(waterfetch,datasubset)
######## Subsetting to second subset######################
waterfetch <- subsetWF(waterfetch,secondsubset)

# Specifying climate variables to run regression against
climatevar_CPC <- 'tmax'
climatevar_chirp <- 'precip'

############# Scaling precip ########################
lags <- c('6d7dl','29d30dl','89d90dl',
          '179d180dl','364d365dl')

allprecip <- paste('precip',lags,sep="")
allspei <- paste('surplus',lags,sep="")
precipLengths <- c(7,30,90,180,365)
precipWeeks <- precipLengths/7

for (var in 1:length(allprecip)){
  waterfetch[,allprecip[var]] <- waterfetch[,allprecip[var]]/precipWeeks[var]
  waterfetch[,allspei[var]] <- waterfetch[,allspei[var]]/precipWeeks[var]
}


# Preparing SFD data - for both chirp and CPC data
# Lat & lon spacing of CPC data (tmin, tmax, and CDD)
CPClat <- seq(-30.75,22.25,0.5)
CPClon <- seq(-17.75,50.75,0.5)

# Lat & lon spacing of chirp data (precip, surplus)
chirplat <- seq(-30.625,21.875,0.5)
chirplon <- seq(-17.625,50.625,0.5)

# allyear <- sort(unique(waterfetch$year))
allsurvey <- unique(waterfetch$survey)


sfd.data.cpc.lon <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,sfd_prep,waterfetch,'s1','lon',CPClat,CPClon)),simplify=FALSE)
sfd.data.cpc.lat <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,sfd_prep,waterfetch,'s1','lat',CPClat,CPClon)),simplify=FALSE)

sfd.data.chirp.lon <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,sfd_prep,waterfetch,'s1','lon',chirplat,chirplon)),simplify=FALSE)
sfd.data.chirp.lat <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,sfd_prep,waterfetch,'s1','lat',chirplat,chirplon)),simplify=FALSE)

#########################################
# Scaling climate variables
climatevar <- c("tmax","precip")
allclimatevar.sfd <- unlist(lapply(climatevar,paste,lags,".sfd",sep=""))

allresult.lat <- list()
allresult.lon <- list()
##################################
# Running regressions
# Running batchlm for each lag
for (rep in 1:reps){
  sfd.lon_cpc_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204.sfd,
                                                       data=sfd.data.cpc.lon[[rep]],timelag=lags[1],direction='lon'))
  sfd.lat_cpc_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204.sfd,
                                                       data=sfd.data.cpc.lat[[rep]],timelag=lags[1],direction='lat'))
  sfd.lon_chirp_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204.sfd,
                                                         data=sfd.data.chirp.lon[[rep]],timelag=lags[1],direction='lon'))
  sfd.lat_chirp_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204.sfd,
                                                         data=sfd.data.chirp.lat[[rep]],timelag=lags[1],direction='lat'))
  
  sfd.lon_cpc_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204.sfd,
                                                         data=sfd.data.cpc.lon[[rep]],timelag=lags[2],direction='lon'))
  sfd.lat_cpc_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204.sfd,
                                                         data=sfd.data.cpc.lat[[rep]],timelag=lags[2],direction='lat'))
  sfd.lon_chirp_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204.sfd,
                                                           data=sfd.data.chirp.lon[[rep]],timelag=lags[2],direction='lon'))
  sfd.lat_chirp_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204.sfd,
                                                           data=sfd.data.chirp.lat[[rep]],timelag=lags[2],direction='lat'))
  
  sfd.lon_cpc_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204.sfd,
                                                         data=sfd.data.cpc.lon[[rep]],timelag=lags[3],direction='lon'))
  sfd.lat_cpc_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204.sfd,
                                                         data=sfd.data.cpc.lat[[rep]],timelag=lags[3],direction='lat'))
  sfd.lon_chirp_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204.sfd,
                                                           data=sfd.data.chirp.lon[[rep]],timelag=lags[3],direction='lon'))
  sfd.lat_chirp_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204.sfd,
                                                           data=sfd.data.chirp.lat[[rep]],timelag=lags[3],direction='lat'))
  
  sfd.lon_cpc_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204.sfd,
                                                           data=sfd.data.cpc.lon[[rep]],timelag=lags[4],direction='lon'))
  sfd.lat_cpc_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204.sfd,
                                                           data=sfd.data.cpc.lat[[rep]],timelag=lags[4],direction='lat'))
  sfd.lon_chirp_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204.sfd,
                                                             data=sfd.data.chirp.lon[[rep]],timelag=lags[4],direction='lon'))
  sfd.lat_chirp_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204.sfd,
                                                             data=sfd.data.chirp.lat[[rep]],timelag=lags[4],direction='lat'))
  
  sfd.lon_cpc_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204.sfd,
                                                           data=sfd.data.cpc.lon[[rep]],timelag=lags[5],direction='lon'))
  sfd.lat_cpc_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204.sfd,
                                                           data=sfd.data.cpc.lat[[rep]],timelag=lags[5],direction='lat'))
  sfd.lon_chirp_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204.sfd,
                                                             data=sfd.data.chirp.lon[[rep]],timelag=lags[5],direction='lon'))
  sfd.lat_chirp_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204.sfd,
                                                             data=sfd.data.chirp.lat[[rep]],timelag=lags[5],direction='lat'))
  
  
  allsfd.lon <- rbind(sfd.lon_cpc_6d7dl, sfd.lon_chirp_6d7dl,
                      sfd.lon_cpc_29d30dl, sfd.lon_chirp_29d30dl,
                      sfd.lon_cpc_89d90dl, sfd.lon_chirp_89d90dl,
                      sfd.lon_cpc_179d180dl, sfd.lon_chirp_179d180dl,
                      sfd.lon_cpc_364d365dl, sfd.lon_chirp_364d365dl)
  
  allsfd.lat <- rbind(sfd.lat_cpc_6d7dl, sfd.lat_chirp_6d7dl,
                      sfd.lat_cpc_29d30dl, sfd.lat_chirp_29d30dl,
                      sfd.lat_cpc_89d90dl, sfd.lat_chirp_89d90dl,
                      sfd.lat_cpc_179d180dl, sfd.lat_chirp_179d180dl,
                      sfd.lat_cpc_364d365dl, sfd.lat_chirp_364d365dl)
  allresult.lon[[rep]] <- allsfd.lon
  allresult.lat[[rep]] <- allsfd.lat
}

allresult.lon <- do.call(rbind.data.frame,allresult.lon)
allresult.lat <- do.call(rbind.data.frame,allresult.lat)

allsfd.lon.summary <- allresult.lon %>%
  group_by(VarName2) %>%
  summarise_each(funs(mean(., na.rm=TRUE)))

allsfd.lon.summary$dir <- 'lon'

allsfd.lat.summary <- allresult.lat %>%
  group_by(VarName2) %>%
  summarise_each(funs(mean(., na.rm=TRUE)))

allsfd.lat.summary$dir <- 'lat'


lon.Estimate2.Std <- allresult.lon %>%
  group_by(VarName2) %>%
  summarise(sigma=sd(Estimate2,na.rm=TRUE))
allsfd.lon.summary$Estimate2.Std.sigma <- lon.Estimate2.Std$sigma

lat.Estimate2.Std <- allresult.lat %>%
  group_by(VarName2) %>%
  summarise(sigma=sd(Estimate2,na.rm=TRUE))
allsfd.lat.summary$Estimate2.Std.sigma <- lat.Estimate2.Std$sigma

# Write to csv
sink(filename)
write.table(allsfd.lon.summary,sep=",",row.names=F)
write.table(allsfd.lat.summary,sep=",", col.names = F,row.names=F)
sink()