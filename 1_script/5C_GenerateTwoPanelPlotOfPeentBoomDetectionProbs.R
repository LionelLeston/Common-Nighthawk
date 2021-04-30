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
library(rsample)
library(broom)
library(purrr)
library(viridis)
set.seed(27)

my.theme <- theme_classic() +
  theme(text=element_text(size=20, family="Arial"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))


#peent detections after back-transforming occupancy and detection probabilities
probs.P<-read.csv("3_output/data/5A_PeentBootstraps_NullOccModel/allpeentoccupancymodelresults.detandoccprobs.csv", header=TRUE)
str(probs.P)
probs.P<-probs.P[probs.P$visits>2,]
probs.P.just<-probs.P[,2:101]
probs.P.ci <- t(apply(probs.P.just, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
probs.P.plusci<-cbind(probs.P,probs.P.ci)
probs.P.plusci<-data.frame(probs.P.plusci)

occ.P<-probs.P.plusci[probs.P.plusci$X=="psi(Int)",]
occ.P$visits.f<-as.factor(occ.P$visits)
occ.P$duration.f<-as.factor(occ.P$duration)

det.P<-probs.P.plusci[probs.P.plusci$X=="p(Int)",]
det.P$visits.f<-as.factor(det.P$visits)
det.P$duration.f<-as.factor(det.P$duration)

peent.detectplot<-ggplot(det.P, aes(x=duration, y=X50., color=visits.f)) +
  geom_point(size=2)+my.theme+
  xlab("Duration of visit")+
  ylab("Peent detection probability")+
  ylim(0,0.8)+
  geom_hline(yintercept=0.5)+
  scale_color_viridis("Visits per site", discrete=TRUE)+ 
  geom_smooth(aes(group = visits.f), method = "glm", formula=y~log(x), se = FALSE)

peent.occprobplot<-ggplot(occ.P, aes(x=duration, y=X50., color=visits.f)) +
  geom_point(size=2)+my.theme+
  xlab("Duration of visit")+
  ylab("Peent occupancy probability")+
  ylim(0,1)+
  geom_hline(yintercept=0.5)+
  scale_color_viridis("Visits per site", discrete=TRUE)+ 
  geom_smooth(aes(group = visits.f), method = "glm", formula=y~log(x), se = FALSE)

peent.detectplotB<-ggplot(det.P, aes(x=duration, y=X50., group=visits.f, color=visits.f)) +
  geom_line(aes(col=visits.f)) +
  geom_hline(yintercept=0.5)+
  xlab("Duration of visit")+
  ylab("Peent detection probability")+
  ylim(0,0.8)+
  scale_color_viridis("Visits per site", discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  my.theme+ 
  theme(legend.position = "none")+
  geom_ribbon(data=det.P,aes(ymin=X5., ymax=X95.,
                           fill=visits.f), alpha=0.3)

peent.occprobplotB<-ggplot(occ.P, aes(x=duration, y=X50., group=visits.f, color=visits.f)) +
  geom_line(aes(col=visits.f)) +
  geom_hline(yintercept=0.5)+
  xlab("Duration of visit")+
  ylab("Peent occupancy probability")+
  ylim(0,1)+
  scale_color_viridis("Visits per site", discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  my.theme+ 
  theme(legend.position = "none")+
  geom_ribbon(data=occ.P,aes(ymin=X5., ymax=X95.,
                             fill=visits.f), alpha=0.3)

#and boom detections
probs.B<-read.csv("3_output/data/5B_BoomBootstraps_NullOccModel/allBoomoccupancymodelresults.detandoccprobs.csv", header=TRUE)
str(probs.B)
probs.B<-probs.B[probs.B$visits>2,]
probs.B.just<-probs.B[,2:101]
probs.B.ci <- t(apply(probs.B.just, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
probs.B.plusci<-cbind(probs.B,probs.B.ci)
probs.B.plusci<-data.frame(probs.B.plusci)

occ.B<-probs.B.plusci[probs.B.plusci$X=="psi(Int)",]
occ.B$visits.f<-as.factor(occ.B$visits)
occ.B$duration.f<-as.factor(occ.B$duration)

det.B<-probs.B.plusci[probs.B.plusci$X=="p(Int)",]
det.B$visits.f<-as.factor(det.B$visits)
det.B$duration.f<-as.factor(det.B$duration)


boom.detectplot<-ggplot(det.B, aes(x=duration, y=X50., color=visits.f)) +
  geom_point(size=2)+my.theme+
  xlab("Duration of visit")+
  ylab("Boom detection probability")+
  ylim(0,0.8)+
  geom_hline(yintercept=0.5)+
  scale_color_viridis("Visits per site", discrete=TRUE)+ 
  geom_smooth(aes(group = visits.f), method = "glm", formula=y~log(x), se = FALSE)

boom.occprobplot<-ggplot(occ.B, aes(x=duration, y=X50., color=visits.f)) +
  geom_point(size=2)+my.theme+
  xlab("Duration of visit")+
  ylab("Boom occupancy probability")+
  ylim(0,1)+
  geom_hline(yintercept=0.5)+
  scale_color_viridis("Visits per site", discrete=TRUE)+ 
  geom_smooth(aes(group = visits.f), method = "glm", formula=y~log(x), se = FALSE)

boom.detectplotB<-ggplot(det.B, aes(x=duration, y=X50., group=visits.f, color=visits.f)) +
  geom_line(aes(col=visits.f)) +
  geom_ribbon(data=det.B,aes(ymin=X5., ymax=X95.,
                           fill=visits.f), alpha=0.3)+
  geom_hline(yintercept=0.5)+
  xlab("Duration of visit")+
  ylab("Boom detection probability")+
  ylim(0,0.8)+
  scale_color_viridis("Visits per site", discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  my.theme+ 
  theme(legend.position = "none")

boom.occprobplotB<-ggplot(occ.B, aes(x=duration, y=X50., group=visits.f, color=visits.f)) +
  geom_line(aes(col=visits.f)) +
  geom_ribbon(data=occ.B,aes(ymin=X5., ymax=X95.,
                             fill=visits.f), alpha=0.3)+
  geom_hline(yintercept=0.5)+
  xlab("Duration of visit")+
  ylab("Boom occupancy probability")+
  ylim(0,1)+
  scale_color_viridis("Visits per site", discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)+ 
  my.theme+ 
  theme(legend.position = "none")

tiff("3_output/figures/5_Occupancy Models/peentandbootdetectionsidebyside_S3.tiff", units="in", width=16, height=8, res=300)
grid.arrange(peent.detectplot, boom.detectplot, 
             ncol=2, nrow=1)
dev.off()

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(peent.detectplot)

p1<-arrangeGrob(peent.detectplot+ theme(legend.position="none"), 
                boom.detectplot+ theme(legend.position="none"), 
                peent.detectplotB, 
                boom.detectplotB, 
                ncol=2, nrow=2)
p2<-arrangeGrob(mylegend)
tiff("3_output/figures/5_Occupancy Models/peentandbootdetection2x2sidebyside_S3.tiff", units="in", width=20, height=16, res=300)
grid.arrange(p1, p2, widths=c(10,1))
dev.off()


mylegend<-g_legend(peent.occprobplot)

g1<-arrangeGrob(peent.occprobplot+ theme(legend.position="none"), 
                boom.occprobplot+ theme(legend.position="none"), 
                peent.occprobplotB, 
                boom.occprobplotB, 
                ncol=2, nrow=2)
g2<-arrangeGrob(mylegend)
tiff("3_output/figures/5_Occupancy Models/peentandbootoccprob2x2sidebyside_26April2021.tiff", units="in", width=20, height=16, res=300)
grid.arrange(g1, g2, widths=c(10,1))
dev.off()