#############################
# Run all SFD regressions
# Running regressions by randomly selecting one for each grid
# and iterating n times
# Feb 25, 2020 Updated Nov 2, 2020
# Abby Harvey

# Must run from same directory that contains waterfetch.rds file

rm(list=ls())
# Overall parameters to change
datasubset <- "all" # enter either urban, rural, elec, noelec, or all
secondsubset <- "all" # Should be different than datasubset, or "all"
chirpsize <- 0.50 #enter either 0.25 for default, or 0.5 for larger grid sizes
chirpmaxgap <- 1 #enter any number between 1 and ~20
cpcmaxgap <- 1 #enter any number between 1 and ~20, but generally 1 works best
reps <- 1 #If using s1 method, specify number of replicates; otherwise, value does not matter
type <- "s1" #don't change

# This file will subset data as specified,
# run the SFD regressions,
# and output one csv file containing all results

# The filename used to save results will be:
# sfd_type_datasubset_chirpsize_chirpmaxgap_cpcmaxgap.csv

if (identical(secondsubset,"all")){
  filename <- paste("maxpairs",datasubset,chirpsize,reps,sep="_")
} else {
  filename <- paste("maxpairs",datasubset,secondsubset,chirpsize,"nit",reps,sep="_")  
}
filename <- paste(filename,"csv",sep=".")
# By default, cluster standard errors at clusterid and lat or lon channel
# And standard error calculation method is using vcovHC from sandwich, HC3

# for parameter values, CASE MATTERS - DO NOT CAPITALIZE

# Load in required packages
library('data.table')
library('sandwich')
library('lmtest')
library('dplyr')
options(dplyr.summarise.inform = FALSE)

# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

##### Max pairs differencing ######
max_pairs_diff <- function(survey_year,waterfetch,direction,CPClat,CPClon){
  
  survey1 <- filter(waterfetch,survey==survey_year)
  
  # Specifying climate variables to keep
  climatevar <- c('tmax','precip')
  lags <- c('6d7dl','13d7dl','29d30dl','59d30dl','89d90dl','179d180dl','364d365dl')
  allvar.regress <- c("hv204",unlist(lapply(climatevar,paste,lags,sep="")),"hv025","new_hv206")
  # Keeping other descriptive vars
  keyvars <- c("latcell","loncell","survey","datetime")
  
  # mark cells
  survey1$latcell <- cut(survey1$latnum,CPClat,labels=FALSE)
  survey1$loncell <- cut(survey1$longnum,CPClon,labels=FALSE)
  
  if (direction=='lon'){
    # identify cells with adjacent neighbors
    test <- survey1 %>% group_by(latcell,loncell) %>% summarize(n=n())
    test$latcelldiff <- test$latcell-lag(test$latcell,1)
    test$loncelldiff <- test$loncell-lag(test$loncell,1)
    
    coords_diff <- test %>% filter(latcelldiff==0&loncelldiff==1) %>%
      select(latcell,loncell)
    
    if (nrow(coords_diff)==0){
      return()
    }
    
    wf_output <- data.frame()
    # build database that grows with every iteration?
    for (i in 1:nrow(coords_diff)){
      coords2 <- coords_diff[i,]
      coords1 <- c(coords2$latcell,coords2$loncell-1)
      
      wf_coords1 <- filter(survey1,latcell==coords2$latcell&loncell==coords2$loncell)
      wf_coords2 <- filter(survey1,latcell==coords1[1],loncell==coords1[2])
      
      npairs <- min(c(nrow(wf_coords1),nrow(wf_coords2)))
      
      key_vars <- wf_coords1[1,] %>% select(all_of(keyvars))
      wf1 <- wf_coords1 %>%
        sample_n(.,npairs)
      wf2 <- wf_coords2 %>%
        sample_n(.,npairs)
      
      wf_todiff1 <- wf1 %>% select(all_of(allvar.regress))
      wf_todiff2 <- wf2 %>% select(all_of(allvar.regress))
      
      
      # differencing
      wf_diff <- wf_todiff2 - wf_todiff1
      
      temp_out <- data.frame(key_vars,wf_diff,wf1$latnum,wf1$longnum,wf2$latnum,wf2$longnum)
      wf_output <- rbind(wf_output,temp_out)
    }
  } else if (direction=='lat'){
    # identify cells with adjacent neighbors
    test <- survey1 %>% group_by(loncell,latcell) %>% summarize(n=n())
    test$latcelldiff <- test$latcell-lag(test$latcell,1)
    test$loncelldiff <- test$loncell-lag(test$loncell,1)
    
    coords_diff <- test %>% filter(latcelldiff==1&loncelldiff==0) %>%
      select(latcell,loncell)
    
    if (nrow(coords_diff)==0){
      return()
    }
    wf_output <- data.frame()
    # build database that grows with every iteration?
    for (i in 1:nrow(coords_diff)){
      coords2 <- coords_diff[i,]
      coords1 <- c(coords2$latcell-1,coords2$loncell)
      
      wf_coords1 <- filter(survey1,latcell==coords2$latcell&loncell==coords2$loncell)
      wf_coords2 <- filter(survey1,latcell==coords1[1],loncell==coords1[2])
      
      npairs <- min(c(nrow(wf_coords1),nrow(wf_coords2)))
      
      key_vars <- wf_coords1[1,] %>% select(all_of(keyvars))
      wf1 <- wf_coords1 %>%
        sample_n(.,npairs)
      wf2 <- wf_coords2 %>%
        sample_n(.,npairs)
      
      wf_todiff1 <- wf1 %>% select(all_of(allvar.regress))
      wf_todiff2 <- wf2 %>% select(all_of(allvar.regress))
      
      
      # differencing
      wf_diff <- wf_todiff2 - wf_todiff1
      
      temp_out <- data.frame(key_vars,wf_diff,wf1$latnum,wf1$longnum,wf2$latnum,wf2$longnum)
      wf_output <- rbind(wf_output,temp_out)
    }
  }
  return(wf_output)
}

sfd.batchlm <- function(climatevar,dependent,data,timelag,direction){
  #This is a function to run batch regressions using lm on SFD data
  if (identical(direction,'lon')){errordir <- 'latcell'
  } else if (identical(direction,'lat')){errordir <- 'loncell'
  } else {print("Invalid direction - must be lat or lon")
    break}
  errors <- c("survey",errordir)
  lmresult <- data.frame()
  
  climate <- paste(climatevar,timelag,sep="")
  rframe <- data.frame(dependent,data[,climate])#create data frame for regressions
  names(rframe)[2] <- 'climatevar1'
  # Run lm
  model <- lm(dependent~climatevar1, data=rframe)
  modelvcov <- coeftest(model, vcov=vcovHC(model,type="HC3",cluster=errors))
  # Save output as row in data frame
  rownames(modelvcov)[2] <- climate
  # 
  # temp_stats <- cbind(summary(model)$df[2],summary(model)$r.squared,summary(model)$adj.r.squared,summary(model)$sigma,summary(model)$fstatistic[1],mean(abs(data$daydiff)))
  # 
  # sink(filename,append=T)
  # write.table(modelvcov,sep=",",row.names=F)
  # write.table(temp_stats)
  # sink()
  
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
  temp2 <- cbind(summary(model)$df[2],summary(model)$r.squared,summary(model)$adj.r.squared,summary(model)$sigma,summary(model)$fstatistic[1],modelsum)
  colnames(temp2)[1:6] <- c('n','R_squared','Adj_r_squared','sigma','fstatistic','name1')
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
chirplat <- seq(-30.625,21.875,chirpsize)
chirplon <- seq(-17.625,50.625,chirpsize)

# allyear <- sort(unique(waterfetch$year))
allsurvey <- unique(waterfetch$survey)

sfd.data.cpc.lon <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,max_pairs_diff,waterfetch,'lon',CPClat,CPClon)),simplify=FALSE)
sfd.data.cpc.lat <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,max_pairs_diff,waterfetch,'lat',CPClat,CPClon)),simplify=FALSE)

sfd.data.chirp.lon <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,max_pairs_diff,waterfetch,'lon',chirplat,chirplon)),simplify=FALSE)
sfd.data.chirp.lat <- replicate(reps,do.call(rbind.data.frame,lapply(allsurvey,max_pairs_diff,waterfetch,'lat',chirplat,chirplon)),simplify=FALSE)

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
  sfd.lon_cpc_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204,
                                                       data=sfd.data.cpc.lon[[rep]],timelag=lags[1],direction='lon'))
  sfd.lat_cpc_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204,
                                                       data=sfd.data.cpc.lat[[rep]],timelag=lags[1],direction='lat'))
  sfd.lon_chirp_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204,
                                                         data=sfd.data.chirp.lon[[rep]],timelag=lags[1],direction='lon'))
  sfd.lat_chirp_6d7dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204,
                                                         data=sfd.data.chirp.lat[[rep]],timelag=lags[1],direction='lat'))
  
  sfd.lon_cpc_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204,
                                                         data=sfd.data.cpc.lon[[rep]],timelag=lags[2],direction='lon'))
  sfd.lat_cpc_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204,
                                                         data=sfd.data.cpc.lat[[rep]],timelag=lags[2],direction='lat'))
  sfd.lon_chirp_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204,
                                                           data=sfd.data.chirp.lon[[rep]],timelag=lags[2],direction='lon'))
  sfd.lat_chirp_29d30dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204,
                                                           data=sfd.data.chirp.lat[[rep]],timelag=lags[2],direction='lat'))
  
  sfd.lon_cpc_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204,
                                                         data=sfd.data.cpc.lon[[rep]],timelag=lags[3],direction='lon'))
  sfd.lat_cpc_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204,
                                                         data=sfd.data.cpc.lat[[rep]],timelag=lags[3],direction='lat'))
  sfd.lon_chirp_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204,
                                                           data=sfd.data.chirp.lon[[rep]],timelag=lags[3],direction='lon'))
  sfd.lat_chirp_89d90dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204,
                                                           data=sfd.data.chirp.lat[[rep]],timelag=lags[3],direction='lat'))
  
  sfd.lon_cpc_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204,
                                                           data=sfd.data.cpc.lon[[rep]],timelag=lags[4],direction='lon'))
  sfd.lat_cpc_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204,
                                                           data=sfd.data.cpc.lat[[rep]],timelag=lags[4],direction='lat'))
  sfd.lon_chirp_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204,
                                                             data=sfd.data.chirp.lon[[rep]],timelag=lags[4],direction='lon'))
  sfd.lat_chirp_179d180dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204,
                                                             data=sfd.data.chirp.lat[[rep]],timelag=lags[4],direction='lat'))
  
  sfd.lon_cpc_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lon[[rep]]$hv204,
                                                           data=sfd.data.cpc.lon[[rep]],timelag=lags[5],direction='lon'))
  sfd.lat_cpc_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_CPC,sfd.batchlm,dependent=sfd.data.cpc.lat[[rep]]$hv204,
                                                           data=sfd.data.cpc.lat[[rep]],timelag=lags[5],direction='lat'))
  sfd.lon_chirp_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lon[[rep]]$hv204,
                                                             data=sfd.data.chirp.lon[[rep]],timelag=lags[5],direction='lon'))
  sfd.lat_chirp_364d365dl <- do.call(rbind.data.frame,lapply(climatevar_chirp,sfd.batchlm,dependent=sfd.data.chirp.lat[[rep]]$hv204,
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
  summarize(n=mean(n),
            R_squared=mean(R_squared),
            Adj_r_squared=mean(Adj_r_squared),
            sigma=mean(sigma),
            fstatistic=mean(fstatistic),
            intercept=mean(Estimate1,na.rm=T),
            interceptsd=sd(Estimate1,na.rm=T),
            poly1=(mean(Estimate2, na.rm=T)),
            poly1sd=sd(Estimate2,na.rm=T),
            poly2=mean(Estimate3,na.rm=T),
            poly2sd=sd(Estimate3,na.rm=T),
            poly3=mean(Estimate4,na.rm=T),
            poly3sd=sd(Estimate4,na.rm=T))

allsfd.lon.summary$dir <- 'lon'

allsfd.lat.summary <- allresult.lat %>%
  group_by(VarName2) %>%
  summarize(n=mean(n),
            R_squared=mean(R_squared),
            Adj_r_squared=mean(Adj_r_squared),
            sigma=mean(sigma),
            fstatistic=mean(fstatistic),
            intercept=mean(Estimate1,na.rm=T),
            interceptsd=sd(Estimate1,na.rm=T),
            poly1=(mean(Estimate2, na.rm=T)),
            poly1sd=sd(Estimate2,na.rm=T),
            poly2=mean(Estimate3,na.rm=T),
            poly2sd=sd(Estimate3,na.rm=T),
            poly3=mean(Estimate4,na.rm=T),
            poly3sd=sd(Estimate4,na.rm=T))

allsfd.lat.summary$dir <- 'lat'

# Write to csv
write.table(allsfd.lon.summary,filename,sep=",",row.names=F)
write.table(allsfd.lat.summary,filename,sep=",", col.names = F,row.names=F,append=TRUE)