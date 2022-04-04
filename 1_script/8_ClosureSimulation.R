library(tidyverse)
library(unmarked)
library(data.table)
library(gridExtra)

sites <- 100 
visits <- 20

#A. RANDOM----

#1. Create data----
#100% occupancy

detection <- sample.int(n=visits, size=sites, replace=TRUE)

history <- data.frame()
for(i in 1:length(detection)){
  int <- detection[i]
  history.i <- c(rep(0,int-1), 1, rep(0, visits-int))
  history <- rbind(history, history.i)
  
}


#2. Model----
pred <- list()
for(j in 2:visits){
  
  history.j <- history[,1:j]
  dat.occ <- unmarkedFrameOccu(history.j)
  mod <- occu(~1 ~1, data=dat.occ)
  
  pred[[j]] <- history.j %>% 
#  pred.test <- history.j %>%   
    cbind(predict(mod, type="det", newdata=dat.occ)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
    cbind(predict(mod, type="state", newdata=dat.occ)) %>% 
    dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    unique() %>% 
    summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
              occu.mn.sd=sd(Occu, na.rm=TRUE),
              occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
              occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
              det.mn.mn=mean(Det, na.rm=TRUE),
              det.mn.sd=sd(Det, na.rm=TRUE),
              det.lwr.mn=mean(DetLower, na.rm=TRUE),
              det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
    mutate(visits=j)
  
}

predictions.random <- rbindlist(pred)

#3. Plot----
plot.random.det <- ggplot(predictions.random) +
  geom_point(aes(x=visits, y=det.mn.mn)) +
  ggtitle("Random") +
  ylab("p(detectability)") +
  xlab("")

plot.random.occu <- ggplot(predictions.random) +
  geom_point(aes(x=visits, y=occu.mn.mn)) +
  ylab("p(occupancy)")


#B. ORDERED----

#1. Create data----
#100% occupancy

detection <- sample.int(n=20, size=sites, replace=TRUE)

history <- data.frame()
for(i in 1:length(detection)){
  int <- detection[i]
  history.i <- c(rep(1,int-1), 1, rep(0, visits-int))
  history <- rbind(history, history.i)
  
}


#2. Model----
pred <- list()
for(j in 2:visits){
  
  cols <- data.frame(col=c(1:20)) %>% 
    sample_n(j)
  history.j <- history[,cols$col]
  dat.occ <- unmarkedFrameOccu(history.j)
  mod <- occu(~1 ~1, data=dat.occ)
  
  pred[[j]] <- history.j %>% 
    #  pred.test <- history.j %>%   
    cbind(predict(mod, type="det", newdata=dat.occ)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
    cbind(predict(mod, type="state", newdata=dat.occ)) %>% 
    dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    unique() %>% 
    summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
              occu.mn.sd=sd(Occu, na.rm=TRUE),
              occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
              occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
              det.mn.mn=mean(Det, na.rm=TRUE),
              det.mn.sd=sd(Det, na.rm=TRUE),
              det.lwr.mn=mean(DetLower, na.rm=TRUE),
              det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
    mutate(visits=j)
  
}

predictions.order <- rbindlist(pred)

#3. Plot----
plot.order.det <- ggplot(predictions.order) +
  geom_point(aes(x=visits, y=det.mn.mn)) +
  ggtitle("Ordered") +
  ylab("") +
  xlab("")

plot.order.occu <- ggplot(predictions.order) +
  geom_point(aes(x=visits, y=occu.mn.mn)) +
  ylab("")

grid.arrange(plot.order.det, plot.order.occu)

grid.arrange(plot.random.det, plot.order.det, plot.random.occu, plot.order.occu, nrow=2, ncol=2)
