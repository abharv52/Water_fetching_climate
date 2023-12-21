#############################
# Code to run all spatial first differences regressions
# Abby Harvey
# Updated Dec 21, 2023

# Overall parameters to change
datasubset <- "all" # enter either urban, rural, elec, noelec, or all
secondsubset <- "all" # Should be different than datasubset, or "all"
numcores <- 1 # specifying number of cores to run on - default should be 1 core

# This file will subset data as specified,
# run the FE regressions,
# and output one csv file containing all results

if (!identical(secondsubset,"all")){
  filename <- paste("FE_",datasubset,"_",secondsubset,".csv",sep="")
} else {
  filename <- paste("FE_",datasubset,".csv",sep="")
}

# Load in required packages
library('data.table')
library('sandwich')
library('lmtest')
library('dplyr')
library(doParallel)

# Load in water-fetch data
waterfetch <- readRDS('waterfetch.rds')

lags <- c('6d7dl','29d30dl','89d90dl','179d180dl','364d365dl')

#### Adding in function to generate dummy (0,1) variables for all quantitative data
dummyvar <- function(input_var,colremove=FALSE){
  ## This is function to create a data frame of dummy variables from an input column
  # Makes each unique value in input_var into a column
  # If colremove=TRUE then the function will remove one column for use in Fixed Effects Regressions
  uniq_val <- unique(input_var) #Finding unique values in input_var
  dummy_mat <- matrix(0,nrow=length(input_var),ncol=length(uniq_val)) #initializing matrix to store dummy variables
  count = 1
  for (val in uniq_val){
    dummy_mat[,count] <- as.numeric(input_var == val)
    count = count+1
  }
  dummy_frame <- data.frame(dummy_mat)
  colnames(dummy_frame) <- uniq_val
  
  if (colremove){
    dummy_frame <- dummy_frame[,1:ncol(dummy_frame)-1]
  }
  
  return(dummy_frame)
}

# Defining function to run linear regression for fixed effects
# Use the ... to enter as many fixed effects as needed
batchlm <- function(climatevar,dependent,data,timelag,...){
  #This is a function to run batch regressions using lm
  lmresult <- data.frame()
  
  climate <- paste(climatevar,timelag,sep="")
  rframe <- data.frame(dependent,data[,climate],...)#create data frame for regressions
  # Run lm
  model <- lm(dependent~., data=rframe)
  modelvcov <- coeftest(model, vcov=vcovHC(model,type="HC3",cluster="clusterid"))
  # Save output as row in data frame
  rownames(modelvcov)[2] <- climate
  modelsum <- data.frame(row.names=1)
  count=1
  for (term in rownames(modelvcov)){
    temp <- data.frame(term,t(modelvcov[term,]))
    colnames(temp)[1] <- 'VarName'
    colnames(temp) <- paste(colnames(temp),count,sep="")
    modelsum <- cbind(modelsum,temp)
    count=count+1
  }
  temp2 <- cbind(summary(model)$df[2],summary(model)$r.squared,summary(model)$adj.r.squared,summary(model)$sigma,summary(model)$fstatistic[1],modelsum)
  colnames(temp2)[1:5] <- c('n','R_squared','Adj_r_squared','sigma','fstatistic')
  lmresult <- structure(rbind(lmresult,temp2),.Names=names(temp2))
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


########### Scaling precip variable ####################
allprecip <- paste('precip',lags,sep="")
precipLengths <- c(7,30,90,180,365)
precipWeeks <- precipLengths/7

for (var in 1:length(allprecip)){
  waterfetch[,allprecip[var]] <- waterfetch[,allprecip[var]]/precipWeeks[var]
}


# Scaling climate variables
climatevar <- c("tmax","precip")

# Run dummyvar on country to get dummy variables for each country
dummy_country <- dummyvar(waterfetch$dhscc,colremove=TRUE)
# Run dummyvar on month to get dummy variables for each month
dummy_month <- dummyvar(waterfetch$month,colremove=TRUE)
# Run dummyvar on year to get dummy variables for each year
dummy_year <- dummyvar(waterfetch$year,colremove=TRUE)

# Running batchlm for each lag
climate_6d7dl <- do.call(rbind.data.frame,mclapply(climatevar,batchlm,dependent=waterfetch$hv204,data=waterfetch,
                                                   timelag=lags[1],dummy_country,dummy_month,dummy_year,mc.cores=numcores))
climate_29d30dl <- do.call(rbind.data.frame,mclapply(climatevar,batchlm,dependent=waterfetch$hv204,data=waterfetch,
                                                     timelag=lags[2],dummy_country,dummy_month,dummy_year,mc.cores=numcores))
climate_89d90dl <- do.call(rbind.data.frame,mclapply(climatevar,batchlm,dependent=waterfetch$hv204,data=waterfetch,
                                                     timelag=lags[3],dummy_country,dummy_month,dummy_year,mc.cores=numcores))
climate_179d180dl <- do.call(rbind.data.frame,mclapply(climatevar,batchlm,dependent=waterfetch$hv204,data=waterfetch,
                                                       timelag=lags[4],dummy_country,dummy_month,dummy_year,mc.cores=numcores))
climate_364d365dl <- do.call(rbind.data.frame,mclapply(climatevar,batchlm,dependent=waterfetch$hv204,data=waterfetch,
                                                       timelag=lags[5],dummy_country,dummy_month,dummy_year,mc.cores=numcores))

allclimatelags <- rbind(climate_6d7dl,climate_29d30dl,climate_89d90dl,
                        climate_179d180dl,climate_364d365dl)
write.csv(allclimatelags,file=filename)
