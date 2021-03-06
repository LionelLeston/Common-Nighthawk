library(boot)
library(dplyr)
library(gam)
library(gamm4)
library(ggplot2)
library(grid)
library(gridExtra)
library(jpeg)
library(lme4)
library(lubridate)
library(mgcv)
library(MuMIn)
library(plotly)
library(pROC)
library(tidyr)
library(unmarked)
library(viridis)


my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(margin=margin(16,0,0,0)),
        axis.title.y=element_text(margin=margin(0,16,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

# function to obtain random samples stratified by site 
#cannot stratify by 2-week periods because after excluding sites
#without temperature readings one site (13577)
#only had samples from 2 nights within 1 2-week period
fSampleNint <- function(dat, intervals, n) {
  intervals <- enquo(intervals)
  dat %>%
    group_by(sm_id)%>%
    filter(UQ(intervals) %in% sample(unique(UQ(intervals)), n)) %>%
    slice(sample(row_number()))
}#dropped season because after removing nights without temperature
#at least one site only had 2 nights available.
#I could add twilight period to ensure even distribution of samples
#throughout the night. If I do so, then 6 sites will have fewer
#samples because they lack intervals from "Night" and "Astronomical"
#periods on nights with temperature data. All 6 sites
#are northern sites

latGroupA<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/latGroupA.csv", header=TRUE)
hist(latGroupA$TSSScorr)
latGroupC<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/latGroupC.csv", header=TRUE)
hist(latGroupC$TSSScorr)
#Check how number of booms counted per interval varies
boomsint<-read.csv("0_data/processed/3_BoomsMapped_SunAndMoon/TempSunMoonBoomDetections.10min.int.csv", header=TRUE)

LATS<-unique(boomsint$latitude)
hist(LATS)
library(BAMMtools)
#2 groups
getJenksBreaks(LATS, 1, subset = NULL)
#45.35199: 1 below, 22 above

#3 groups
getJenksBreaks(LATS, 2, subset = NULL)
#45.35199 62.73417

#4 groups
getJenksBreaks(LATS, 3, subset = NULL)
#45.35199 54.32880 62.73417

getJenksBreaks(LATS, 4, subset = NULL)
#45.35199 48.25271 56.73500 62.73417

boomsint$BOOMS<-ifelse(is.na(boomsint$eventID),0,1)
boomsint$intervals<-as.factor(boomsint$intervals)
boomsint$latitude.F<-ifelse(boomsint$latitude<50,"A",
                             ifelse(boomsint$latitude<55,"B",
                                    ifelse(boomsint$latitude<60,"C","D")))
boomsint$halfmonth<-ifelse(boomsint$day<16,"1","2")
boomsint$season<-as.factor(paste0(boomsint$month,"_",boomsint$halfmonth))

boomsint$start_time_new<- ymd_hms(as.character(boomsint$start_time_new))
boomsint$julian<-yday(as.Date(boomsint$date))
boomsint$hour<-hour(boomsint$start_time_new)
boomsint$minute<-minute(boomsint$start_time_new)
boomsint$start_time_numeric<-boomsint$hour+(boomsint$minute/60)

boomsint$TSSS.s<-scale(boomsint$TSSScorr, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
#use TSSScorr as there were inaccuracies in TSSS
boomsint$ordinal.s<-scale(as.vector(boomsint$julian), center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$latitude.s<-scale(boomsint$latitude, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$longitude.s<-scale(boomsint$longitude, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$altitude.s<-scale(boomsint$altitude, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$meantemp.s<-scale(boomsint$meantemp, center=FALSE, scale=TRUE)#not currently being used in these GAMMs
#boomsint$moon.s<-scale(boomsint$moon.fraction, center=FALSE, scale=TRUE)#not currently being used in these GAMMs


boomsint$sm_id<-as.factor(boomsint$sm_id)
boomsint$intervals<-as.factor(boomsint$intervals)
boomsintB<-boomsint[!boomsint$month==8,]#drop August
boomsintC<-boomsintB[!is.na(boomsintB$meantemp),]#drop intervals missing temperature boomsint
maxdur<-max(boomsintC$duration)
boomsintD<-boomsintC[boomsintC$duration==maxdur,]#gets rid of any "remainder" intervals

#this is the point where I could split into training and testing data for purposes of validation
#training data is what I draw samples from
#the test data is then summarized separately afterwards
#and predictions are generated by creating a test data 
#model matrix and matrix-multiplying it with each column of
#bootstrapped coefficients to get 1 column of predictions
#I then graph the observed:predicted relationship and get
#R-squared or correlation coefficient

#I then get a column of predicted values for each bootstrap that I can compare to actual numbers counted
nrow(boomsintD)
boomsintD$sm_id.interval<-as.factor(paste0(boomsintD$sm_id,boomsintD$intervals))
boomsintD.test <- fSampleNint(boomsintD, intervals, 20)
nrow(boomsintD.test)
write.csv(boomsintD.test, file="0_data/processed/4B_BoomActivityRates_GAMMs/boomsintD.test.csv")
boomsintD.train  <- boomsintD[!boomsintD$sm_id.interval %in% boomsintD.test$sm_id.interval,]
nrow(boomsintD.train)
write.csv(boomsintD.train, file="0_data/processed/4B_BoomActivityRates_GAMMs/boomsintD.train.csv")
#GAMM run once before bootstrap to get model matrix

boomsintD.train<-read.csv("0_data/processed/4B_BoomActivityRates_GAMMs/boomsintD.train.csv", header=TRUE)
boomsintD.train$nightlength<-as.POSIXct(boomsintD.train$SR.time.posix)-as.POSIXct(boomsintD.train$SS.time.posix)
boomsintD.train$TSSS.SR<-boomsintD.train$nightlength*3600

d<-fSampleNint(boomsintD.train, intervals, 20)

booms.pervisit<-d %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(BOOMS), 
            duration=mean(duration_new),
            TSSS=mean(TSSScorr), 
            TSSS.s=mean(TSSS.s), 
            meantemp=mean(meantemp),
            meantemp.s=mean(meantemp.s),
            twilightperiod=names(which.max(table(twilightperiod))),
            latitude.F=names(which.max(table(latitude.F))),
            latitude=mean(latitude),
            latitude.s=mean(latitude.s),
            longitude=mean(longitude),
            longitude.s=mean(longitude.s),
            altitude=mean(altitude),
            altitude.s=mean(altitude.s),
            #moon.fraction=mean(fraction),
            date=names(which.max(table(date))),
            ordinal.day=mean(julian),
            ordinal.s=mean(ordinal.s),
            time=mean(start_time_numeric),
            civil1.start=mean(civil1.start),
            nautical1.start=mean(nautical1.start),
            astro1.start=mean(astro1.start),
            night.start=mean(night.start),
            astro2.start=mean(astro2.start),
            nautical2.start=mean(nautical2.start),
            civil2.start=mean(civil2.start),
            TBSR=mean(TBSR))
write.csv(booms.pervisit, file="0_data/processed/4B_BoomActivityRates_GAMMs/temp/boomspervisit.csv")
write.csv(boomsintD.train, file="0_data/processed/4B_BoomActivityRates_GAMMs/temp/boomspervisittraining.csv")

#GAMM

#basic model for 10-minute boom intervals: run once to get model matrix for making predictions
booms.pervisit$latitude.F<-as.factor(booms.pervisit$latitude.F)

mean(booms.pervisit$numdetect)  
var(booms.pervisit$numdetect)  
thetaval=(var(booms.pervisit$numdetect)/mean(booms.pervisit$numdetect))

#full.GAM.SS <- gamm4(numdetect ~latgroupN + s(meantemp.s, bs = "cs", by=latgroupN) + s(ordinal.s, bs = "cs", by=latgroupN) + s(TSSS.s, bs = "cs", by=latgroupN), random=~(1|sm_id), data=booms.pervisit, family=negbin(theta=thetaval, link="log"))
#full.GAM.SS <- gamm4(numdetect ~latitude.F + s(meantemp.s, k = 4, bs = "cs", by=latitude.F) + s(ordinal.s, k = 4, bs = "cs", by=latitude.F) + s(TSSS.s, k = 4, bs = "cs", by=latitude.F), random=~(1|sm_id), data=booms.pervisit, family=negbin(theta=thetaval, link="log"))
#plot.gam(full.GAM.SS$gam, main = "Smoothed Untransformed Boom Activity Estimates")
#coef.gamm<-coef(full.GAM.SS$gam)#40 coefficients
full.GAM.SS <- gamm4(numdetect ~latitude.F + s(ordinal.s, k = 3, bs = "cs", by=latitude.F) + s(TSSS.s, k = 5, bs = "cs", by=latitude.F), random=~(1|sm_id), data=booms.pervisit, family=negbin(theta=thetaval, link="log"))
coef.gammB<-coef(full.GAM.SS$gam)
model.matrix(full.GAM.SS$gam)
#"(Intercept)","latitude.FB","latitude.FC","latitude.FD", 
#"s(ordinal.s):latitude.FA.1","s(ordinal.s):latitude.FA.2","s(ordinal.s):latitude.FB.1","s(ordinal.s):latitude.FB.2", 
#"s(ordinal.s):latitude.FC.1","s(ordinal.s):latitude.FC.2","s(ordinal.s):latitude.FD.1","s(ordinal.s):latitude.FD.2", 
#"s(TSSS.s):latitude.FA.1","s(TSSS.s):latitude.FA.2","s(TSSS.s):latitude.FA.3","s(TSSS.s):latitude.FA.4", 
#"s(TSSS.s):latitude.FB.1","s(TSSS.s):latitude.FB.2","s(TSSS.s):latitude.FB.3","s(TSSS.s):latitude.FB.4", 
#"s(TSSS.s):latitude.FC.1","s(TSSS.s):latitude.FC.2","s(TSSS.s):latitude.FC.3","s(TSSS.s):latitude.FC.4", 
#"s(TSSS.s):latitude.FD.1","s(TSSS.s):latitude.FD.2","s(TSSS.s):latitude.FD.3","s(TSSS.s):latitude.FD.4" 

df<-data.frame(coef.gamm)#28 coefficients
#write.csv(df, file="3_output/data/4A_BoomGAMMS/GAMMbooms.10min.int.csv")

#Now create the bootstrap
bs <- function(data, indices){
  d<-fSampleNint(data, intervals, 20)
  
  booms.pervisit<-d %>%
    group_by(sm_id, intervalsP) %>% 
    summarize(numdetect= sum(BOOMS), 
              duration=mean(duration_new),
              TSSS=mean(TSSScorr), 
              TSSS.s=mean(TSSS.s), 
              meantemp=mean(meantemp),
              meantemp.s=mean(meantemp.s),
              twilightperiod=names(which.max(table(twilightperiod))),
              latitude=mean(latitude),
              latitude.F=names(which.max(table(latitude.F))),
              latitude.s=mean(latitude.s),
              longitude=mean(longitude),
              longitude.s=mean(longitude.s),
              altitude=mean(altitude),
              altitude.s=mean(altitude.s),
              #moon.fraction=mean(fraction),
              date=names(which.max(table(date))),
              ordinal.day=mean(julian),
              ordinal.s=mean(ordinal.s),
              time=mean(start_time_numeric))
  
  #GAMM

  #basic model for 10-minute boom intervals
  booms.pervisit$latitude.F<-as.factor(booms.pervisit$latitude.F)
  mean(booms.pervisit$numdetect)  
  var(booms.pervisit$numdetect)  
  thetaval=(var(booms.pervisit$numdetect)/mean(booms.pervisit$numdetect))
  
  try(full.GAM.SS <- gamm4(numdetect ~latitude.F  + s(ordinal.s, k = 3, bs = "cs", by=latitude.F) + s(TSSS.s, k = 5, bs = "cs", by=latitude.F), random=~(1|sm_id), data=booms.pervisit, family=negbin(theta=thetaval, link="log")))
  #plot.gam(full.GAM.SS$gam, main = "Smoothed Untransformed Boom Activity Estimates")
  try(coef.gamm<-coef(full.GAM.SS$gam))
  try(return(coef.gamm))
}

#boom.all<-boomsintD.train#read.csv("0_data/processed/3_BoomsMapped_SunAndMoon/TempSunMoonBoomDetections.10min.int.csv", header=TRUE)

boomresults <- boot(data=boomsintD.train,
                    statistic=bs, parallel="multicore",
                    R=100)
#print(paste0("bootstraps successfully run for ",fn))

CoefB <-t(boomresults$t)
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
save(boomresults, CoefB, CoefBCI50, CoefBCI90, file="3_output/data/4B_March2021_BootstrappedBoomActivityGAMMs/boot20samples10minuteintBooms.RData")
#save as csv
df<-data.frame(CoefB)
row.names(df) = c("(Intercept)","latitude.FB","latitude.FC","latitude.FD", 
                  "s(ordinal.s):latitude.FA.1","s(ordinal.s):latitude.FA.2","s(ordinal.s):latitude.FB.1","s(ordinal.s):latitude.FB.2", 
                  "s(ordinal.s):latitude.FC.1","s(ordinal.s):latitude.FC.2","s(ordinal.s):latitude.FD.1","s(ordinal.s):latitude.FD.2", 
                  "s(TSSS.s):latitude.FA.1","s(TSSS.s):latitude.FA.2","s(TSSS.s):latitude.FA.3","s(TSSS.s):latitude.FA.4", 
                  "s(TSSS.s):latitude.FB.1","s(TSSS.s):latitude.FB.2","s(TSSS.s):latitude.FB.3","s(TSSS.s):latitude.FB.4", 
                  "s(TSSS.s):latitude.FC.1","s(TSSS.s):latitude.FC.2","s(TSSS.s):latitude.FC.3","s(TSSS.s):latitude.FC.4", 
                  "s(TSSS.s):latitude.FD.1","s(TSSS.s):latitude.FD.2","s(TSSS.s):latitude.FD.3","s(TSSS.s):latitude.FD.4" 
)
write.csv(df, file="3_output/data/4B_March2021_BootstrappedBoomActivityGAMMs/boot20samples10minuteintBooms.csv")
#print(paste0("bootstrap coefficients saved for ",fn))

## box plot of variables with confidence intervals
prednames = c("(Intercept)","latitude.FB","latitude.FC","latitude.FD", 
              "s(ordinal.s):latitude.FA.1","s(ordinal.s):latitude.FA.2","s(ordinal.s):latitude.FB.1","s(ordinal.s):latitude.FB.2", 
              "s(ordinal.s):latitude.FC.1","s(ordinal.s):latitude.FC.2","s(ordinal.s):latitude.FD.1","s(ordinal.s):latitude.FD.2", 
              "s(TSSS.s):latitude.FA.1","s(TSSS.s):latitude.FA.2","s(TSSS.s):latitude.FA.3","s(TSSS.s):latitude.FA.4", 
              "s(TSSS.s):latitude.FB.1","s(TSSS.s):latitude.FB.2","s(TSSS.s):latitude.FB.3","s(TSSS.s):latitude.FB.4", 
              "s(TSSS.s):latitude.FC.1","s(TSSS.s):latitude.FC.2","s(TSSS.s):latitude.FC.3","s(TSSS.s):latitude.FC.4", 
              "s(TSSS.s):latitude.FD.1","s(TSSS.s):latitude.FD.2","s(TSSS.s):latitude.FD.3","s(TSSS.s):latitude.FD.4" 
)
var.bci<-CoefBCI90
var.bci<-data.frame(var.bci)
var.bci$median<-as.numeric(var.bci$X50.)
var.bci$lcl<-as.numeric(var.bci$X5.)
var.bci$ucl<-as.numeric(var.bci$X95.)
var.bci$prednames<-prednames

GG<-ggplot(var.bci, aes(x=prednames, y=median))+
  geom_point(aes(x=prednames, y=median))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl))+
  geom_hline(yintercept=0)+
  xlab("Predictor")+
  ylab("Effect on booms counted")+coord_flip()+my.theme
ggsave("3_output/figures/4B_March2021_BootstrappedBoomActivityGAMMs/boot20visits10minuteintervalBoom.jpeg", plot=GG, width=12, height=6, units=c("in"), dpi=300)

#To get GAMM predictions within a bootstrap:
#https://gist.github.com/davharris/3da3a43f3d45ce01ee20

## make data for prediction - TSSS
load("3_output/data/4B_March2021_BootstrappedBoomActivityGAMMs/boot20samples10minuteintBooms.RData")
#gives us CoefB

library(ggplot2)
library(viridis)

#Creating datasets, then matrix plots of predicted activity (one curve per latitude and date)
booms<-read.csv("0_data/processed/4B_BoomActivityRates_GAMMs/temp/boomspervisittraining.csv", header=TRUE)
dates<-unique(booms[,c("date","julian","ordinal.s")])#"ordinal.day",
dates<-data.frame(dates)
#gives us a list of dates, ordinal days, and scaled values
#to loop through

sunrisetimes<-booms%>%
  group_by(date,latitude.F)%>%
  summarise(TSSS.SR=mean(TSSS.SR))
sunrisetimes<-data.frame(sunrisetimes)
sunrisetimes$latitude.F<-as.factor(sunrisetimes$latitude.F)
#gives us a list of mean sunrise times on different dates for
#sites in different latitude groups

range(booms$julian)#152 212
range(booms$TSSScorr)#-3600 34200_

for (i in 1:nrow(dates)){
  Day<-dates[i,1]
  ordinal.day<-dates[i,2]
  boom.tsss.jd<-expand.grid(intercept=1,
                             latitude.F=c("A","B","C","D"), 
                             TSSS=seq(from=-3600,to=34800,by=(34800--3600)/100),
                             OrdinalDay=seq(from=152,to=212,by=2)
  )
  boom.tsss.jd$latitude.FB<-ifelse(boom.tsss.jd$latitude.F=="B",1,0)
  boom.tsss.jd$latitude.FC<-ifelse(boom.tsss.jd$latitude.F=="C",1,0)
  boom.tsss.jd$latitude.FD<-ifelse(boom.tsss.jd$latitude.F=="D",1,0)
  Latitude<-boom.tsss.jd$latitude.F
  #boom.tsss.jd$ordinal.day<-dates[i,2]#
  boom.tsss.jd$ordinal.s<-scale(boom.tsss.jd$OrdinalDay, center=FALSE, scale=TRUE)#scaled value associated with ordinal day
  boom.tsss.jd$TSSS.s<-scale(boom.tsss.jd$TSSS, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
  #boom.tsss.jd$LxO<-boom.tsss.jd$latitude.s*boom.tsss.jd$Ordinal.orth1#boom.tsss.jd$sm_id<-"7524"
  #boom.tsss.jd$LxT<-boom.tsss.jd$latitude.s*boom.tsss.jd$TSSS.orth1#boom.tsss.jd$sm_id<-"7524"
  #boom.tsss.jd$LxT2<-boom.tsss.jd$latitude.s*boom.tsss.jd$TSSS.orth2#boom.tsss.jd$sm_id<-"7524"
  #boom.tsss.jd$LxO2<-boom.tsss.jd$latitude.s*boom.tsss.jd$Ordinal.orth2#boom.tsss.jd$sm_id<-"7524"
  #boom.tsss.jd$sm_id<-as.factor(boom.tsss.jd$sm_id)
  
  #Now limit predictions to those from a single date
  boom.tsss.jd<-boom.tsss.jd[boom.tsss.jd$OrdinalDay==ordinal.day,]
  
  #set unscaled latitude and TSSS aside for now and remove these
  #variables from prediction dataset prior to creating fixed effects
  #design matrix
  OrdinalDay<-boom.tsss.jd$OrdinalDay
  #ORDINAL.S<-boom.tsss.jd$ordinal.s
  #boom.tsss.jd$OrdinalDay<-NULL
  #boom.tsss.jd$ordinal.s<-NULL
  #boom.tsss.jd$ordinal.day<-NULL
  Latitude<-boom.tsss.jd$latitude.F
  #boom.tsss.jd$latitude.F<-NULL
  TSSS<-boom.tsss.jd$TSSS
  #boom.tsss.jd$TSSS<-NULL
  #boom.tsss.jd$TSSS.s<-NULL
  str(boom.tsss.jd)
  names(boom.tsss.jd)
  
  #boom.tsss.jdmm<-as.matrix(boom.tsss.jd)
  boom.tsss.jdmm <- model.matrix(full.GAM.SS$gam, boom.tsss.jd)#meantemp.s+moon.s+
  
  ## predict and apply inverse link function: this gives an n_new x B+1 matrix
  bootpreds100 <- (exp(boom.tsss.jdmm %*% CoefB))#take antilog?
  
  ## calculate median and 95% CIs for predictions
  bp100.CI95 <- t(apply(bootpreds100, 1, quantile, seq(from=0.025, to=0.975, by=((0.975-0.025)/95)), na.rm=TRUE))
  
  ## Use these bootstrapped predictions to generate heat maps, 3d plots, etc.
  
  boom.tsss.jd.bci<-cbind(boom.tsss.jd, bp100.CI95)
  boom.tsss.jd.F<-cbind(boom.tsss.jd.bci,Latitude,TSSS,OrdinalDay)#put Latitude and TSSS back in
  boom.tsss.jd.F$TSSS.min<-boom.tsss.jd.F$TSSS/60
  boom.tsss.jd.F$TSSS.hr<-boom.tsss.jd.F$TSSS/3600
  boom.tsss.jd.F$Latitude.F<-as.factor(boom.tsss.jd.F$Latitude)
  boom.tsss.jd.F$date<-Day

  levels(boom.tsss.jd.F$Latitude.F)
  latlevels<-c("A","B","C","D")
  for (j in 1:length(latlevels)){
    preds.sitelevel<-boom.tsss.jd.F[boom.tsss.jd.F$Latitude.F==latlevels[j],]
    preds.025.975<-preds.sitelevel[,10:96]
    for (k in 1:ncol(preds.025.975)){
      preds.025.975[,k]<-ifelse(preds.025.975[,k]>600,600,preds.025.975[,k])
    }
    sunrisetimesB<-sunrisetimes[sunrisetimes$latitude.F==latlevels[j],]
    sunrisetimesC<-sunrisetimesB[sunrisetimesB$date==Day,]
    if (nrow(sunrisetimesC)>0){
      preds.sitelevel$TSSS.SR<-sunrisetimesC$TSSS.SR
      preds.sitelevel$TSSS.SR.hr<-preds.sitelevel$TSSS.SR/3600
    }
    if (nrow(sunrisetimesC)==0){
      preds.sitelevel$TSSS.SR.hr<-0#if no sunrise times are available for that combination of date and latitude
    }
    tiff(paste0("3_output/figures/4B_Bootstrapped BoomGAMMs/Matrix Plots booms X TSSS X lat different days March 15/",ifelse(j==1,"A_",ifelse(j==2,"B_",ifelse(j==3,"C_","D_"))),"BoomvsTSSS-",preds.sitelevel$Latitude.F,"-",preds.sitelevel$OrdinalDay,".tiff"), units="in", width=6, height=4, res=300)
    #X <- c(0,unique(preds.sitelevel$TSSS.SR.hr))
    #Y <- c(0,max(preds.025.975))
    matplot(preds.sitelevel$TSSS.hr, 
            preds.025.975, 
            type="l", 
            col="blue",#cividis(100),
            xlab=paste0("Hours After Sunset ", Day),
            #xlim=rev(range(preds.sitelevel$TSSS.hr)),
            ylab=paste0("Boom Activity (",ifelse(j==1,"45-50 deg. N ", ifelse(j==2,"50-55 deg. N",ifelse(j==3,"55-60 deg. N","60-65 deg. N"))),")"))   
    abline(v=c(0,unique(preds.sitelevel$TSSS.SR.hr)), col=c("red", "red"), lty=c(1,2), lwd=c(1, 3))
    #lim <- par("usr")
    #rect(X[1], Y[1], X[2], Y[2], border = "red", col = NA)
    #axis(1) ## add axes back
    #axis(2)
    #box()
    dev.off()
  }
}


#Get multipanel plot of peents vs. TSSS
library(magick)
library(dplyr)
library(tidyr)
library(magrittr)
# read images and then create a montage
# tile =2 , means arrange the images in 2 columns
# geometry controls the pixel sixe, spacing between each image in the collage output. 
# read the the png files into a list
pngfiles <-
  list.files(
    path = "3_output/figures/4B_Bootstrapped BoomGAMMs/Matrix Plots booms X TSSS X lat different days March 15",
    recursive = TRUE,
    pattern = "\\.png$",
    full.names = T
  )

magick::image_read(pngfiles) %>%
  magick::image_montage(tile = "3", geometry = "x500+10+5") %>%
  magick::image_convert("jpg") %>%
  magick::image_write(
    format = ".png", path = "3_output/figures/4B_Bootstrapped BoomGAMMs/Figure5_boomsXtsss_collage.png",
    quality = 100
  )

#Matrix plots: Booms X different days 1 hour after sunset at different latitudes
boom.tsss.jd<-expand.grid(intercept=1,
                           latitude.F=c("A","B","C","D"), 
                           TSSS=3600,
                           OrdinalDay=seq(from=152,to=212,by=2)
)
boom.tsss.jd$latitude.FB<-ifelse(boom.tsss.jd$latitude.F=="B",1,0)
boom.tsss.jd$latitude.FC<-ifelse(boom.tsss.jd$latitude.F=="C",1,0)
boom.tsss.jd$latitude.FD<-ifelse(boom.tsss.jd$latitude.F=="D",1,0)
Latitude<-boom.tsss.jd$latitude.F
boom.tsss.jd$ordinal.s<-scale(boom.tsss.jd$OrdinalDay, center=FALSE, scale=TRUE)#scaled value associated with ordinal day
boom.tsss.jd$TSSS.s<-0.216293523#value associated with TSSS=3600 in training data

#set unscaled latitude and TSSS aside for now and remove these
#variables from prediction dataset prior to creating fixed effects
#design matrix
OrdinalDay<-boom.tsss.jd$OrdinalDay
Latitude<-boom.tsss.jd$latitude.F
TSSS<-boom.tsss.jd$TSSS
str(boom.tsss.jd)
names(boom.tsss.jd)

#boom.tsss.jdmm<-as.matrix(boom.tsss.jd)
boom.tsss.jdmm <- model.matrix(full.GAM.SS$gam, boom.tsss.jd)#meantemp.s+moon.s+

## predict and apply inverse link function: this gives an n_new x B+1 matrix
bootpreds100 <- (exp(boom.tsss.jdmm %*% CoefB))#take antilog?

## calculate median and 95% CIs for predictions
bp100.CI95 <- t(apply(bootpreds100, 1, quantile, seq(from=0.025, to=0.975, by=((0.975-0.025)/95)), na.rm=TRUE))

## Use these bootstrapped predictions to generate heat maps, 3d plots, etc.

boom.tsss.jd.bci<-cbind(boom.tsss.jd, bp100.CI95)
boom.tsss.jd.F<-cbind(boom.tsss.jd.bci,Latitude,TSSS,OrdinalDay)#put Latitude and TSSS back in
boom.tsss.jd.F$TSSS.min<-boom.tsss.jd.F$TSSS/60
boom.tsss.jd.F$TSSS.hr<-boom.tsss.jd.F$TSSS/3600
boom.tsss.jd.F$Latitude.F<-as.factor(boom.tsss.jd.F$Latitude)

levels(boom.tsss.jd.F$Latitude.F)
latlevels<-c("A","B","C","D")
for (j in 1:length(latlevels)){
  preds.sitelevel<-boom.tsss.jd.F[boom.tsss.jd.F$Latitude.F==latlevels[j],]
  preds.025.975<-preds.sitelevel[,10:96]
  for (k in 1:ncol(preds.025.975)){
    preds.025.975[,k]<-ifelse(preds.025.975[,k]>600,600,preds.025.975[,k])
  }
  tiff(paste0("3_output/figures/4B_Bootstrapped BoomGAMMs/Matrix Plots booms X ordinal day TSSS 3600/",ifelse(j==1,"A",ifelse(j==2,"B_",ifelse(j==3,"C","D"))),"BoomvsDate-",preds.sitelevel$Latitude.F,"-",preds.sitelevel$OrdinalDay,".tiff"), units="in", width=6, height=4, res=300)
  Y <- c(0,max(preds.025.975))
  matplot(preds.sitelevel$OrdinalDay, 
          preds.025.975, 
          type="l", 
          col="blue",#cividis(100),
          xlab=paste0("Ordinal Day"),
          ylab=paste0("Boom Activity (",ifelse(j==1,"45-50 deg. N ", ifelse(j==2,"50-55 deg. N",ifelse(j==3,"55-60 deg. N","60-65 deg. N"))),")"))   
  dev.off()
}


#Get multipanel plot of peents vs. ordinal day
library(magick)
library(dplyr)
library(tidyr)
library(magrittr)
# read images and then create a montage
# tile =2 , means arrange the images in 2 columns
# geometry controls the pixel sixe, spacing between each image in the collage output. 
# read the the png files into a list
pngfiles <-
  list.files(
    path = "3_output/figures/4B_Bootstrapped BoomGAMMs/Matrix Plots booms X ordinal day TSSS 3600",
    recursive = TRUE,
    pattern = "\\.png$",
    full.names = T
  )

magick::image_read(pngfiles) %>%
  magick::image_montage(tile = "2", geometry = "x500+10+5") %>%
  magick::image_convert("jpg") %>%
  magick::image_write(
    format = ".png", path = "3_output/figures/4B_Bootstrapped BoomGAMMs/Figure7_boomsXordinalday_collage.png",
    quality = 100
  )


#Validating the GAMMS, using the test data from earlier
#i.e. boomsintD.test, which still must be summarized by interval
boomsintD.test<-read.csv("0_data/processed/4B_BoomActivityRates_GAMMs/boomsintD.test.csv") 
booms.pervisit.test<-boomsintD.test %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(BOOMS), 
            duration=mean(duration_new),
            TSSS=mean(TSSScorr), 
            TSSS.s=mean(TSSS.s), 
            meantemp=mean(meantemp),
            meantemp.s=mean(meantemp.s),
            twilightperiod=names(which.max(table(twilightperiod))),
            latitude=mean(latitude),
            latitude.s=mean(latitude.s),
            latitude.F=names(which.max(table(latitude.F))),
            longitude=mean(longitude),
            #moon.fraction=mean(fraction),
            date=names(which.max(table(date))),
            ordinal.day=mean(julian),
            ordinal.s=mean(ordinal.s),
            time=mean(start_time_numeric))

booms.pervisit.test$latitude.FB<-ifelse(booms.pervisit.test$latitude.F=="B",1,0)
booms.pervisit.test$latitude.FC<-ifelse(booms.pervisit.test$latitude.F=="C",1,0)
booms.pervisit.test$latitude.FD<-ifelse(booms.pervisit.test$latitude.F=="D",1,0)

## model matrix that'll match the coefficients
Xnew.Test.mm <- model.matrix(full.GAM.SS$gam, booms.pervisit.test)

## predict and apply inverse link function: this gives an n_new x B+1 matrix
Preds <- (exp(Xnew.Test.mm %*% CoefB))

booms.pervisit.test.df<-data.frame(booms.pervisit.test)
Preds.df<-data.frame(Preds)
TestDatapreds<-cbind(booms.pervisit.test.df, Preds.df)

str(TestDatapreds)
write.csv(TestDatapreds, file="3_output/data/4B_March2021_BootstrappedBoomActivityGAMMs/TestDatapreds.csv")

spearmancorlist<-list()
names<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
         "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
         "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30",
         "X31","X32","X33","X34","X35","X36","X37","X38","X39","X40",
         "X41","X42","X43","X44","X45","X46","X47","X48","X49","X50",
         "X51","X52","X53","X54","X55","X56","X57","X58","X59","X60",
         "X61","X62","X63","X64","X65","X66","X67","X68","X69","X70",
         "X71","X72","X73","X74","X75","X76","X77","X78","X79","X80",
         "X81","X82","X83","X84","X85","X86","X87","X88","X89","X90",
         "X91","X92","X93","X94","X95","X96","X97","X98","X99","X100")
for (i in names){
  TestDatapreds$Predicted<-TestDatapreds[,i]
  spearmancorlist[[i]]<-cor(TestDatapreds$numdetect, TestDatapreds$Predicted, method="spearman")
}
spearmancorvec<-unlist(spearmancorlist)
spearmandf<-data.frame(spearmancorvec)
## calculate median and 90% CIs
RhoN <-quantile(spearmancorvec, c(0.5, 0.05, 0.95), na.rm=TRUE)
#       50%         5%        95% 
#0.13488234 0.02502464 0.22827643 (in 2020 using 2 latitude groups)
#0.2520381  0.1848156  0.2981680  (in 2021 using 4 latitude groups)




