library(boot)
library(dplyr)
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

#Graph predicted peent counts by site and twilight period
TestDatapreds.P<-read.csv("3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/TestDatapredsNolat.csv", header=TRUE)
TestDatapreds.P$twilightperiodB<-as.factor(TestDatapreds.P$twilightperiodB)
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="BEFORE"] <- "Before Sunset"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="CIVIL.1"] <- "Civil - PM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="CIVIL.2"] <- "Civil - AM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="ASTRONOMICAL.1"] <- "Astro - PM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="ASTRONOMICAL.2"] <- "Astro - AM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="NAUTICAL.1"] <- "Nautical - PM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="NAUTICAL.2"] <- "Nautical - AM"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="AFTER"] <- "After Sunrise"
levels(TestDatapreds.P$twilightperiodB)[levels(TestDatapreds.P$twilightperiodB)=="NIGHT"] <- "Night"

TestDatapreds.P.g<-TestDatapreds.P%>%
  gather(X7:X100, key="bootstrap", value="predictedcount")

TestDatapreds.P.g$twilightperiodB<-as.factor(TestDatapreds.P.g$twilightperiodB)

TestDatapreds.P.g$twilightperiodB<-factor(TestDatapreds.P.g$twilightperiodB, levels=c("Before Sunset","Civil - PM","Nautical - PM","Astro - PM","Night","Astro - AM","Nautical - AM","Civil - AM","After Sunrise"))
levels(TestDatapreds.P.g$twilightperiodB)

TestDatapreds.P.g$predictedcount.r<-round(TestDatapreds.P.g$predictedcount, digits=1)
# TestDatapreds.P.g$predictedcount.r[TestDatapreds.P.g$twilightperiodB=="Before Sunset"]<-0.005
# TestDatapreds.P.g$predictedcount.r[TestDatapreds.P.g$twilightperiodB=="Civil - PM"]<-0.005
# TestDatapreds.P.g$predictedcount.r[TestDatapreds.P.g$twilightperiodB=="After Sunrise"]<-0.005
# Basic violin plot - peents
peentplot <- ggplot(TestDatapreds.P.g, aes(x=twilightperiodB, y=predictedcount.r, col=twilightperiodB, fill=twilightperiodB)) + 
  geom_violin()+
  my.theme+
  ylab("Predicted # peents per 10 min")+
  ylim(0,0.6)+
  xlab("Time Period")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_viridis(discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  theme(legend.position = "none")+ 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4)


#Graph predicted boom counts by site and twilight period
TestDatapreds.B<-read.csv("3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/TestDatapredsNolat.csv", header=TRUE)
TestDatapreds.B$twilightperiodB<-as.factor(TestDatapreds.B$twilightperiodB)
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="BEFORE"] <- "Before Sunset"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="CIVIL.1"] <- "Civil - PM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="CIVIL.2"] <- "Civil - AM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="ASTRONOMICAL.1"] <- "Astro - PM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="ASTRONOMICAL.2"] <- "Astro - AM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="NAUTICAL.1"] <- "Nautical - PM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="NAUTICAL.2"] <- "Nautical - AM"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="AFTER"] <- "After Sunrise"
levels(TestDatapreds.B$twilightperiodB)[levels(TestDatapreds.B$twilightperiodB)=="NIGHT"] <- "Night"


TestDatapreds.B.g<-TestDatapreds.B%>%
  gather(X5:X100, key="bootstrap", value="predictedcount")

TestDatapreds.B.g$twilightperiodB<-as.factor(TestDatapreds.B.g$twilightperiodB)

TestDatapreds.B.g$twilightperiodB<-factor(TestDatapreds.B.g$twilightperiodB, levels=c("Before Sunset","Civil - PM","Nautical - PM","Astro - PM","Night","Astro - AM","Nautical - AM","Civil - AM","After Sunrise"))
levels(TestDatapreds.B.g$twilightperiodB)

TestDatapreds.B.g$predictedcount.r<-round(TestDatapreds.B.g$predictedcount, digits=1)
# TestDatapreds.B.g$predictedcount.r[TestDatapreds.B.g$twilightperiodB=="Astro - PM"]<-0.005
# TestDatapreds.B.g$predictedcount.r[TestDatapreds.B.g$twilightperiodB=="Night"]<-0.005
# TestDatapreds.B.g$predictedcount.r[TestDatapreds.B.g$twilightperiodB=="After Sunrise"]<-0.005

# Basic violin plot - booms
boomplot <- ggplot(TestDatapreds.B.g, aes(x=twilightperiodB, y=predictedcount.r, col=twilightperiodB, fill=twilightperiodB)) + 
  geom_violin()+
  my.theme+
  ylab("Predicted # booms per 10 min")+
  ylim(0,0.6)+
  xlab("Time Period")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_viridis(discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  theme(legend.position = "none")+ 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4)


tiff("3_output/figures/6_Twilight Period Mixed Models/peentandbootviolinsidebyside_26April2021.tiff", units="in", width=16, height=8, res=300)
grid.arrange(peentplot, boomplot, 
                        ncol=2, nrow=1)
dev.off()