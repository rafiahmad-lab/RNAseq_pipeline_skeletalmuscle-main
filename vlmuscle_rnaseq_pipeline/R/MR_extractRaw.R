

#### Extract MR CSA estimation #############

### MR data  ####
# Extensor muscle CSA was determined from 4 sequential images.
# A 2nd degree polynominal model was fitted to CSA as a function of location normalized to the patella
# Pre- and post-CSA was determined from the model fit at the same relative location

source("./R/lib_fun.R")

dat <- read_excel("./data/MR.xlsx")
# extract all unique participants
subj<-unique(dat$code)

## create a sample based on code, leg and timepoint
dat$sample<-paste(dat$code,dat$leg,dat$timepoint)

## list all samples in data set
samples<-unique(dat$sample)

## extract maximal locations for all legs and timepoints ## 
locations<-data.frame(sample=rep(NA, length(samples)),
                      ml=rep(NA, length(samples)))

dat<-data.frame(dat)

# Store all maximum locations
for(i in 1:length(samples)){
  
  locations[i,1]<-samples[i]
  locations[i,2]<-max(dat[dat$sample==samples[i],11])
}


## maximum location per leg -- the maximum common location ##
temp<-str_split_fixed(locations$sample, " ", 3)
temp<-paste(temp[,1], temp[,2])
locations$leg<-temp

legs<-unique(locations$leg)

locations2<-data.frame(legs=legs,
                       loc=rep(NA, length(legs)))

for(i in 1:length(legs)){
  locations2[i,1]<-legs[i]
  locations2[i,2]<-min(locations[locations[,3]==legs[i],2])
}

## join maximum locations with original data set
dat$legs<-paste(dat$code, dat$leg)

dat<-inner_join(locations2, dat)

## loop to create for model fit and results extraction for all data ## 
mr.results<-data.frame(sample=rep(NA,length(samples)),
                       CSA.point=rep(NA,length(samples)),
                       CSA.avg=rep(NA,length(samples)),
                       CSA.femur=rep(NA,length(samples)),
                       CSA.raw.avg=rep(NA,length(samples)))



for(i in 1:length(samples)){
  temp<-dat%>%
    filter(sample==samples[i])  # filter out a sample
  
  model<-lm(CSA.ext~poly(normalized.location,3), data=temp) # a 2nd degree polynominal fit
  femur.model<-lm(CSA.femur~poly(normalized.location,2), data=temp) # a 2nd degree polynominal fit

  newdata<-data.frame(normalized.location=c(as.numeric(temp$loc[1])-5,
                                            as.numeric(temp$loc[1])-10,
                                            as.numeric(temp$loc[1])-15, 
                                            as.numeric(temp$loc[1])-20,
                                            as.numeric(temp$loc[1])-25,
                                            as.numeric(temp$loc[1])-30)) # average of 6 
  
  mr.results[i,1]<-samples[i]
  mr.results[i,2]<-predict(model, newdata)[2]
  mr.results[i,3]<-mean(predict(model, newdata))
  mr.results[i,4]<-mean(predict(femur.model, newdata))
  mr.results[i,5]<-mean(temp$CSA.ext, na.rm=T)
  
} 


mr.results<-cbind(str_split_fixed(mr.results$sample, " ", 3), mr.results)
colnames(mr.results)<-c("code", "leg", "timepoint", "sample", "CSA.point", "CSA.avg", "CSA.femur", "CSA.raw.avg")


## read code for identifications ##
code<-read_excel("./data/MR_codeSubject.xlsx")
code<-data.frame(code[,c(2,3)])
code$subject<-gsub("X", "FP", code$subject)
code$code<-toupper(code$code)

mr.results<-data.frame(mr.results)
mr.results$subject<-rep(NA, length(nrow(mr.results)))

for(i in 1:nrow(code)){
  mr.results[mr.results$code==code[i,2],]$subject<-code[i,1]
}

mr.results<-mr.results%>%
  dplyr::select(subject, leg, timepoint, CSA.point, CSA.avg, CSA.femur, CSA.raw.avg)

rm(code);rm(dat);rm(locations);rm(locations2); rm(newdata); rm(temp); rm(femur.model); rm(i);rm(legs);rm(model);rm(samples);rm(subj)

