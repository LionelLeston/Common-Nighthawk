library(boot)
library(dplyr)
library(lme4)
library(lubridate)
library(MuMIn)
library(tidyr)
library(unmarked)
library(rsample)
library(broom)
library(purrr)
library(viridis)
set.seed(27)

# function to obtain regression weights
fSampleNint <- function(dat, intervals, n) {
  intervals <- enquo(intervals)
  dat %>%
    group_by(sm_id)%>%
    filter(UQ(intervals) %in% sample(unique(UQ(intervals)), n)) %>%
    slice(sample(row_number()))
}#remove season if intervals from nights without temperature are
#excluded, and because we're limiting samples to one narrow window
#of time in the summer (Julian = 160-190)


#function to filter data to dates/times of higher activity
#rates, process a sample data into occupancy model format, 
#run occupancy models, and save model coefficients for generating 
#predictions
bs <- function(data, indices) {
  #filter data to best times and dates
  #LAST BIT OF FILTERING: REMOVE ANY SITES THAT HAD NO PEENT DETECTIONS WHATSOEVER
  #SO THAT ALL SITES IN ANALYSIS WERE ACTUALLY OCCUPIED
  data$halfmonth<-ifelse(data$day<16,"1","2")
  data$season<-as.factor(paste0(data$month,"_",data$halfmonth))
  data$intervals<-as.factor(data$intervals)
  data$julian<-yday(as.Date(data$date))
  dataB<-data[which(data$julian<185&data$julian>175),]
  dataC<-dataB[dataB$TSSScorr<3600,]
  maxdur<-max(data$duration)
  dataD<-dataC[dataC$duration==maxdur,]#gets rid of any "remainder" intervals
  #visits sampled from 2-week periods, minimum of 4 visits for occupancy models

  d<-fSampleNint(dataD, intervals, n)#1 sample per season=4 visits sampled

  d$PEENTS<-ifelse(is.na(d$eventID),0,1)
  
  d.withround<-d %>%
    group_by(sm_id) %>% 
    mutate(round= droplevels(intervals) %>% as.numeric)#, season
  #write.csv(d.withround, file="0_data/processed/5_NSamplePEENTS_OccupancyModels/temp/PEENTSintwithround.csv")
  
  peents.perround<-d.withround %>%
    group_by(sm_id, round) %>% 
    summarize(numdetect= sum(PEENTS))#, season
  #at this point, I can convert all detections > 0 to 1's to run
  #occupancy models rather than mixture models
  
  #peents.perround$seasonround<-paste0(peents.perround$season,".",peents.perround$round)
  #despite the name "seasonround" I am not treating the occupancy models I run as 
  #"multi-season" occupancy models
  
  #Do I still want to use tapply or do I use recast or dplyr?
  #peents.ofs.df<-tapply(peents.perround$numdetect, list(peents.perround$sm_id, peents.perround$seasonround), max, na.rm=TRUE)
  #write.csv(peents.ofs.df, file="0_data/processed/5_NSamplePEENTS_OccupancyModels/temp/peents.ofs.df.Nmix.csv")
  
  #output check. At this point we have abundance per species per station visit.
  #We could use this output for mixed models or N-mixture models.
  #Since we are doing single-season occupancy models, we convert abundance to detected
  #(0 or 1) 
  
  peents.perround$detectYN<-ifelse(peents.perround$numdetect>0,1,0)
  peents.ofs.df<-tapply(peents.perround$detectYN, list(peents.perround$sm_id, peents.perround$round), max, na.rm=TRUE)

  #y<-peents.ofs.df[,1:ncol(peents.ofs.df)]#>1 visit
  y<-peents.ofs.df#1 visit
  #n<-rowSums(y[,c(1:ncol(y))], na.rm = TRUE)
  #n.bin<-ifelse(n>1,1,n)
  #siteswithPEENTS<-sum(n.bin)
  #probdet<-siteswithPEENTS/nrow(y)
  #samplenum<-ncol(peents.ofs.df)
  #duration<-mean(peents.perround$duration)

  umf_all <- unmarkedFrameOccu(y=y, siteCovs=NULL,  obsCovs=NULL)
  try(Null <- occu(~1 ~1, data = umf_all))
  try(coefs <- coef(Null))
  try(return(coefs))
}

visits<-c(1)#2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20
duration<-c(2,4,6,8,10,12,14,16,18,20)
for (i in duration) {
  peent.all<-read.csv(paste0("0_data/processed/3_PeentsMapped_SunAndMoon/TempSunMoonPeentDetections.",i,"min.int.csv"), header=TRUE)
  peent.allB<-peent.all[!peent.all$sm_id %in% c("7689","9990","12011","12122","13588"),]#removes sites with zero peent detections all season
  #placeholder
  for (j in visits) {
    n<-j
    # bootstrapping with 100 replications: each drawn data set has j random visits per site
    peentresults <- boot(data=peent.allB,
                         statistic=bs, parallel="multicore",
                         R=100)
    print(paste0("bootstraps successfully run for ",j," ",i,"-minute visits per site"))
 
    CoefB <-t(peentresults$t)
    #get rid of columns with errors
    CoefB<-data.frame(CoefB)
    CoefB[] <- lapply(CoefB, function(x) as.numeric(as.character(x)))
    #delete columns with NA values
    all_na <- function(x) any(!is.na(x))
    CoefB<-CoefB %>% select_if(all_na)
    CoefB<-as.matrix(CoefB)
    ncol(CoefB)
    CoefBCI50 <- t(apply(CoefB, 1, quantile, c(0.5, 0.25, 0.75), na.rm=TRUE))
    CoefBCI90 <- t(apply(CoefB, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
    save(peentresults, CoefB, CoefBCI50, CoefBCI90, file=paste0("3_output/data/10A_PeentBootstraps_Null/boot",i,"minutes",j,"visits.RData"))
    #save as csv
    df<-data.frame(CoefB)
    row.names(df) = c("psi(Int)","p(Int)")
    df$duration<-i
    df$visits<-j
    write.csv(df, file=paste0("3_output/data/5A_PeentBootstraps_NullOccModel/coef",i,"minutes",j,"visits.csv"))
    print(paste0("bootstrap coefficients saved for ",i," duration, ",j," visits"))
  }
}
#Note: failed at 20 visits from 20-minute samples
#I think there weren't twenty 20-minute intervals from some sites after filtering data

#combine results into 1 file for graphing
MASTERLIST =list.files("3_output/data/5A_PeentBootstraps_NullOccModel",pattern="coef")
modelresults<-list()#empty list each time we draw a new number of samples

for (j in MASTERLIST){
  spptdf<-read.csv(paste0("3_output/data/5A_PeentBootstraps_NullOccModel/",j), header=TRUE) #temporary data frame
  modelresults[[j]]<-spptdf#append temporary data frame at each loop iteration to the list
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="3_output/data/5A_PeentBootstraps_NullOccModel/allpeentoccupancymodelresults.csv")

#after back-transforming occupancy and detection probabilities
probs<-read.csv("3_output/data/5A_PeentBootstraps_NullOccModel/allpeentoccupancymodelresults.detandoccprobs.csv", header=TRUE)
str(probs)
probs<-probs[probs$visits>2,]
probs.just<-probs[,2:101]
probs.ci <- t(apply(probs.just, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
probs.plusci<-cbind(probs,probs.ci)
probs.plusci<-data.frame(probs.plusci)

occ<-probs.plusci[probs.plusci$X=="psi(Int)",]
occ$visits.f<-as.factor(occ$visits)
occ$duration.f<-as.factor(occ$duration)

det<-probs.plusci[probs.plusci$X=="p(Int)",]
det$visits.f<-as.factor(det$visits)
det$duration.f<-as.factor(det$duration)

#Graph occupancy probability (script 5C)

