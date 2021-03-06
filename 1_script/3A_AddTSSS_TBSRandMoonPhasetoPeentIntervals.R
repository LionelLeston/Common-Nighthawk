library(lubridate)
library(intrval)
library(dplyr)
library(tidyr)
library(suncalc)

#example
peent.10min.int<-read.csv("0_data/processed/2_PeentDetectionsMapped/PeentDetections.10min.int.csv", header=TRUE)
sites<-read.csv("0_data/raw/SitesVisits.csv", header=TRUE)
nrow(sites)#489
sites<-unique(sites[,c("sm_id","location","latitude","longitude")])
nrow(sites)#23
m1<-merge(peent.10min.int, sites[,c("sm_id","latitude","longitude")], by=c("sm_id"), all.x=TRUE)
str(m1)

#get TSSS
m1$SS.time<-paste0(m1$date," ",m1$SS)
m1$SS.time.posix<-as.POSIXct(m1$SS.time, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
#m1$TSSS<-ifelse(as.POSIXct(m1$start_time_new)<m1$SS.time.posix,
#                -1*(as.POSIXct(m1$start_time_new)-m1$SS.time.posix),
#                as.POSIXct(m1$start_time_new)-m1$SS.time.posix)
if(as.POSIXct(m1$start_time_new)<m1$SS.time.posix){
  m1$TSSS<-(as.POSIXct(m1$start_time_new)-m1$SS.time.posix)*-1
} else {
  m1$TSSS<-as.POSIXct(m1$start_time_new)-m1$SS.time.posix
}
#negative values indicate a "visit" happened before sunset

#get TBSR (time before sun rise)
m1$SR.time.posix<-m1$SS.time.posix+m1$SStoSRend

if(as.POSIXct(m1$start_time_new)<m1$SR.time.posix){
  m1$TBSR<-(m1$SR.time.posix-as.POSIXct(m1$start_time_new))
} else {
  m1$TBSR<-(m1$SR.time.posix-as.POSIXct(m1$start_time_new))*-1
}
#negative values indicate a "visit" happened after sunrise

#get midpoint of night (halfway point between sunset and sunrise)
#We can use this variable to filter "visits" that occur before the midpoint to look at
#effect of TSSS and visits after the midpoint to look at effect of TBSR

m1$nightmidpoint<-m1$SS.time.posix+(m1$SStoSRend/2)

#Once we store night midpoint as a posixCT object, we can compare start time of visit to
#the midpoint and determine if a visit occurred closer to sunset or to sunrise
#Alternatively, I could simply compare TSSS to TBSR and use TSSS if TSSS < TBSR
m1$useTSSS<-ifelse(as.POSIXct(m1$start_time_new)<m1$nightmidpoint, 1, 0)

#get moon illumination
moonphase<-getMoonIllumination(date = as.POSIXct(m1$date))
m1<-cbind(m1,moonphase[,c("fraction","phase","angle")])
#get sun angle
data<-data.frame(date = as.POSIXct(m1$start_time_new), 
                 lat=m1$latitude, 
                 lon=m1$longitude)
sunangle<-getSunlightPosition(data=data)
# sunangle$start_time_new<-as.factor(sunangle$date)
# sunangle$latitude<-sunangle$lat
# sunangle$longitude<-sunangle$lon
# sunangle$date<-NULL
# sunangle$lat<-NULL
# sunangle$lon<-NULL
m2<-cbind(m1,sunangle[,c("altitude","azimuth")])#do not use merge
#write.csv(m2, file="0_data/processed/m2.csv")

#get twilight period, 
#using starttime.posix, 
#c1, n1, a1, nightcalc, 
#a2, n2, c2, endCivil2SRstart
m2$civil1.start<-as.POSIXct(m2$startime)+m2$c1
m2$nautical1.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1
m2$astro1.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1
m2$night.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc
m2$astro2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2
m2$nautical2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2+m2$n2
m2$civil2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2+m2$n2+m2$c2
m2$civil2.end<-as.POSIXct(m2$starttime)+3600+m2$SStoSRend-m2$ms8
  
m2$twilightperiod<-ifelse(as.POSIXct(m2$start_time_new)<m2$civil1.start,"BEFORE",
                          ifelse(as.POSIXct(m2$start_time_new)<m2$nautical1.start,"CIVIL",
                                 ifelse(as.POSIXct(m2$start_time_new)<m2$astro1.start,"NAUTICAL",
                                        ifelse(as.POSIXct(m2$start_time_new)<m2$night.start,"ASTRONOMICAL",
                                               ifelse(as.POSIXct(m2$start_time_new)<m2$astro2.start,"NIGHT",
                                                      ifelse(as.POSIXct(m2$start_time_new)<m2$nautical2.start,"ASTRONOMICAL",
                                                             ifelse(as.POSIXct(m2$start_time_new)<m2$civil2.start,"NAUTICAL",
                                                                    ifelse(as.POSIXct(m2$start_time_new)<m2$civil2.end,"CIVIL","AFTER"))))))))

write.csv(m2, file="0_data/processed/m2b.csv")

#now repeat for all files
ROOT<-getwd()
fl <- list.files(path=paste0(ROOT,"/0_data/processed/2_PeentDetectionsMapped"), pattern=".csv") # filename must have habitat0HF
for (i in 1:length(fl)) {
  fn <- fl[i]
  peent.det<-read.csv(paste0("0_data/processed/2_PeentDetectionsMapped/",fn), header=TRUE)
  print(paste0(fn," read into R"))
  sites<-read.csv("0_data/raw/SitesVisits.csv", header=TRUE)
  sites<-unique(sites[,c("sm_id","latitude","longitude")])#locations written differently
  m1<-merge(peent.det, sites[,c("sm_id","latitude","longitude")], by=c("sm_id"), all.x=TRUE)
  print(paste0(fn," merged with site data"))
  
  #get TSSS
  m1$SS.time<-paste0(m1$date," ",m1$SS)
  m1$SS.time.posix<-as.POSIXct(m1$SS.time, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  m1$start_time.posix<-as.POSIXct(m1$start_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  m1$BeforeSS<-ifelse(m1$start_time.posix<m1$SS.time.posix,1,0)
  m1$TSSS<-ifelse(m1$start_time.posix<m1$SS.time.posix,
                  (m1$SS.time.posix-m1$start_time.posix),
                  (m1$start_time.posix-m1$SS.time.posix))
  m1$TSSScorr<-ifelse(m1$BeforeSS==1,m1$TSSS*-1,m1$TSSS)
  print(paste0("time since sunset added to ",fn))

  #get TBSR (time before sun rise)
  m1$SR.time.posix<-m1$SS.time.posix+m1$SStoSRend
  
  if(as.POSIXct(m1$start_time_new)<m1$SR.time.posix){
    m1$TBSR<-(m1$SR.time.posix-as.POSIXct(m1$start_time_new))
  } else {
    m1$TBSR<-(m1$SR.time.posix-as.POSIXct(m1$start_time_new))*-1
  }
  print(paste0("time before sunrise added to ",fn))
  
  m1$nightmidpoint<-m1$SS.time.posix+(m1$SStoSRend/2)
  m1$useTSSS<-ifelse(as.POSIXct(m1$start_time_new)<m1$nightmidpoint, 1, 0)
  print(paste0("night midpoint filter added to ",fn))
  
  #get moon illumination
  moonphase<-getMoonIllumination(date = as.POSIXct(m1$date))
  m1b<-cbind(m1,moonphase[,c("fraction","phase","angle")])#do not use merge
  print(paste0("moon position added to ",fn))
  
  #get sun angle
  data<-data.frame(date = as.POSIXct(m1$start_time_new), 
                   lat=m1$latitude, 
                   lon=m1$longitude)
  sunangle<-getSunlightPosition(data=data)
  
  m2<-cbind(m1b,sunangle[,c("altitude","azimuth")])#do not use merge
  print(paste0("sun position added to ",fn))
  
  m2$civil1.start<-as.POSIXct(m2$startime)+m2$c1
  m2$nautical1.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1
  m2$astro1.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1
  m2$night.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc
  m2$astro2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2
  m2$nautical2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2+m2$n2
  m2$civil2.start<-as.POSIXct(m2$startime)+m2$c1+m2$n1+m2$a1+m2$nightcalc+m2$a2+m2$n2+m2$c2
  m2$civil2.end<-as.POSIXct(m2$starttime)+3600+m2$SStoSRend-m2$ms8
  
  m2$twilightperiod<-ifelse(as.POSIXct(m2$start_time_new)<m2$civil1.start,"BEFORE",
                            ifelse(as.POSIXct(m2$start_time_new)<m2$nautical1.start,"CIVIL",
                                   ifelse(as.POSIXct(m2$start_time_new)<m2$astro1.start,"NAUTICAL",
                                          ifelse(as.POSIXct(m2$start_time_new)<m2$night.start,"ASTRONOMICAL",
                                                 ifelse(as.POSIXct(m2$start_time_new)<m2$astro2.start,"NIGHT",
                                                        ifelse(as.POSIXct(m2$start_time_new)<m2$nautical2.start,"ASTRONOMICAL",
                                                               ifelse(as.POSIXct(m2$start_time_new)<m2$civil2.start,"NAUTICAL",
                                                                      ifelse(as.POSIXct(m2$start_time_new)<m2$civil2.end,"CIVIL","AFTER"))))))))
  print(paste0("twilight period added to ",fn))
  write.csv(m2, file=paste0("0_data/processed/3_PeentsMapped_SunAndMoon/SunMoon",fn))
}


#So what I've found is that TSSS seems to be correctly calculated, but for some reason,
#the estimate in some files is in minutes, while the estimate in other files is in seconds.
#Also I'm finding that due to how the start times for different twilight periods were 
#calculated previously using suncalc before this contract, civil twilight is calculated
#as starting before sunset for some intervals, while some intervals after sunset
#are calculated as occurring before civil twilight. 

#So some recalculation is necessary. There is no consistency in predicting which
#files will have TSSS in minutes or seconds or which intervals are incorrectly classified 
#to twilight period, so I specified the recalculations on specific files below.

#First, the files for which TSSScorrected needs to be recalculated in seconds, by 
#multiplying by 60, along with the reclassification of twilight periods:

fl2<-c("PeentDetections.7min.int.csv",
       "PeentDetections.8min.int.csv",
       "PeentDetections.9min.int.csv",
       "PeentDetections.11min.int.csv",
       "PeentDetections.13min.int.csv",
       "PeentDetections.14min.int.csv",
       "PeentDetections.16min.int.csv",
       "PeentDetections.17min.int.csv",
       "PeentDetections.18min.int.csv",
       "PeentDetections.19min.int.csv")
for (i in 1:length(fl2)) {
  fn <- fl2[i]
  m3<-read.csv(paste0("0_data/processed/3_PeentsMapped_SunAndMoon/SunMoon",fn), header=TRUE)
  m3$TSSScorr<-m3$TSSScorr*60
  m3$twilightperiod[which(m3$BeforeSS==1&m3$twilightperiod=="CIVIL")]<-"BEFORE"
  m3$twilightperiod[which(m3$BeforeSS==0&m3$twilightperiod=="BEFORE")]<-"CIVIL"
  print(paste0("TSSS and twilight period recalculated for ",fn))
  write.csv(m3, file=paste0("0_data/processed/3_PeentsMapped_SunAndMoon/SunMoon",fn))
}

#Next, the files for which TSSScorrected doesn't need to be recalculated but some
#twilight periods need to be reclassified.
fl3<-c("PeentDetections.1min.int.csv",
       "PeentDetections.2min.int.csv",
       "PeentDetections.3min.int.csv",
       "PeentDetections.4min.int.csv",
       "PeentDetections.5min.int.csv",
       "PeentDetections.6min.int.csv",
       "PeentDetections.10min.int.csv",
       "PeentDetections.12min.int.csv",
       "PeentDetections.15min.int.csv",
       "PeentDetections.20min.int.csv",
       "PeentDetections.1hour.int.csv")
for (i in 1:length(fl3)) {
  fn2 <- fl3[i]
  m3<-read.csv(paste0("0_data/processed/3_PeentsMapped_SunAndMoon/SunMoon",fn2), header=TRUE)
  m3$twilightperiod[which(m3$BeforeSS==1&m3$twilightperiod=="CIVIL")]<-"BEFORE"
  m3$twilightperiod[which(m3$BeforeSS==0&m3$twilightperiod=="BEFORE")]<-"CIVIL"
  print(paste0("TSSS and twilight period recalculated for ",fn2))
  write.csv(m3, file=paste0("0_data/processed/3_PeentsMapped_SunAndMoon/SunMoon",fn2))
}

OldFile<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/OldSunMoonPeentDetections.10min.int.csv", header=TRUE)
NewFile<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/NewSunMoonPeentDetections.10min.int.csv", header=TRUE)

OldFileDetections<-OldFile[!is.na(OldFile$eventTime),]
nrow(OldFileDetections)#61597

NewFileDetections<-NewFile[!is.na(NewFile$eventTime),]
nrow(NewFileDetections)#61597

hist(OldFileDetections$TSSScorr)
hist(NewFileDetections$TSSScorr)#TSSScorr histograms look the same since I added TBSR and useTSSS calculations