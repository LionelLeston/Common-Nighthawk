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
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(margin=margin(16,0,0,0)),
        axis.title.y=element_text(margin=margin(0,16,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

#One of Elly’s comments in a previous draft recommended including 
#a brief, general summary of peent and boom activity at the 
#beginning of the results section. I agree that this would be a 
#useful contribution to the paper and could serve as a relative 
#measure of activity for future surveys/research. Lionel, since 
#you have the most recent and updated version of the database, 
#could you summarize a few of these basic metrics? I think total 
#number of peents overall, total number booms overall, average 
#number peents/booms per site +/- SD, max and min number of peents/
#booms per site would be sufficient. I can embed these results in 
#a few sentences at the start of the results section.

#Check how number of peents counted per interval varies
peentsint<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/TempSunMoonPeentDetections.10min.int.csv", header=TRUE)
peentsint$PEENTS<-ifelse(is.na(peentsint$eventID),0,1)

#Total number of peents overall
sum(peentsint$PEENTS)#61597 peents

#average, sd, min, max
peentsint.sitesumm<-peentsint %>% 
  group_by(sm_id) %>%
  summarise(NumberPerSite=sum(PEENTS))
mean(peentsint.sitesumm$NumberPerSite)#2678.13
sd(peentsint.sitesumm$NumberPerSite)#5727.791
min(peentsint.sitesumm$NumberPerSite)#0
max(peentsint.sitesumm$NumberPerSite)#25096

#Check how number of booms counted per interval varies
boomsint<-read.csv("0_data/processed/3_BoomsMapped_SunAndMoon/TempSunMoonBoomDetections.10min.int.csv", header=TRUE)
boomsint$BOOMS<-ifelse(is.na(boomsint$eventID),0,1)

#Total number of booms overall
sum(boomsint$BOOMS)#12825 booms

#average, sd, min, max
boomsint.sitesumm<-boomsint %>% 
  group_by(sm_id) %>%
  summarise(NumberPerSite=sum(BOOMS))
mean(boomsint.sitesumm$NumberPerSite)#557.6087
sd(boomsint.sitesumm$NumberPerSite)#713.6503
min(boomsint.sitesumm$NumberPerSite)#0
max(boomsint.sitesumm$NumberPerSite)#3087

