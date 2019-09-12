#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%                     SUMMARY OF SIMULATION RESULTS               %#
#% MI for CACE latent models                                       %#
#%                        UPDATE: JuL  2018                        %#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
rm(list=ls())

# Upload necessary packages
library(foreign)
#install.packages("ggplot2","gridExtra","Cairo","extrafont")

library(tidyverse)
library(dplyr)

library(ggplot2)
library(lattice)
library(grid)
library(gridExtra)
library(Cairo)
library(extrafont)
setwd("~/Google Drive/CACE MI")


results_cont<- read.csv("Results/cont_results.csv")
colnames(results_cont)

results_cont$Miss.Y<-factor(results_cont$Miss.Y, labels = c("fully obs", "MAR"))

results_cont$Miss<-interaction(results_cont$Miss.Y,results_cont$Sample.Size)
#results_cont$Scenario<-interaction(results_cont$Noncompliance,results_cont$Effect.size)
results_cont$Scenario<-interaction(results_cont$Sample.Size,results_cont$Noncompliance)
vlabel=c("N=200,30% Noncmp","N=200,50% Noncmp","N=1000,30% Noncmp","N=1000, 50% Noncmp")

results_cont$scenario<-factor(results_cont$Scenario, levels=c("200.30","200.50","1000.30", "1000.50"), labels=vlabel)


data1 <- results_cont %>% select(-(c("Outcome","Unobs.conf")))
data1<- 
  mutate(data1, Noncompliance = as.factor(Noncompliance))

data1<- data1%>% 
mutate(Sample.Size = as.factor(Sample.Size)) %>%
mutate(Effect.size = as.factor(Effect.size))%>%
mutate(Miss.Y = as.factor(Miss.Y)) 
  
data_c_obs<- subset(data1, data1$Miss.Y==0)
data_c_Ymiss<- subset(data1, data1$Miss.Y==1)




#cbPalette <- c("#56B4E9", "#0072B2", "#D55E00", "#000000","darkred")
#Palette5 <- c("#56B4E9", "green", "#D55E00", "#000000","darkred")

Palette <- c("dark grey", "black")


#all together

bias.cont <-
  ggplot(data1, aes(x=Method,y=bias, shape=Miss.Y, colour=Effect.size)) +
  geom_errorbar(aes(ymin=bias-1.96*mc.error, ymax=bias+1.96*mc.error), width=0.4, size=.4, position=position_dodge(0.75)) + 
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_hline(yintercept= 0, linetype="dotted") +  
  labs(  x= "", y = "Mean bias", shape="Continuous outcome")




coverage.cont <-
  ggplot(data1, aes(x=Method,y=coverage, shape=Miss.Y, colour=Effect.size, ymin=0.88,ymax=1)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_y_continuous(breaks=c(.9,.95,1), labels = c("0.90","0.95","1.00")) +
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") + 
  theme_classic() +
  geom_hline(yintercept= 0.95, linetype="dotted") +  
  geom_hline( aes(yintercept= round(.95-1.96*sqrt((.95*.05)/1999),3)), linetype="longdash", color="black", size=0.5)+ 
  geom_hline(aes(yintercept=round(.95+1.96*sqrt((.95*.05)/1999),3)), linetype="longdash", color="black", size=0.5)+ 
  labs(  x= "", y = "Coverage", shape="Continuous outcome")

  

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(bias.cont)


grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom","right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}




#RMSE
rmse.c <-
  ggplot(data=data1, aes(x=Method,y=RMSE, 
                        shape=Miss.Y, colour=Effect.size)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7.5), axis.title=element_text(size=8.2)) +
  xlab("Method") + ylab("RMSE") +
  theme(plot.title = element_text(hjust=0, size=8.2)) 

#CI width
width.c <-
  ggplot(data=data1, aes(x=Method,y=CI.width, 
                         shape=Miss.Y, colour=Effect.size)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7.5), axis.title=element_text(size=8.2)) +
  xlab("Method") + ylab("CI width") +
  theme(plot.title = element_text(hjust=0, size=8.2)) 



#pdf("SMC MIC TeX documents/plot_rmse_contY.pdf")
#grid_arrange_shared_legend(rmse.c, nrow = 1, ncol = 1)
#dev.off()

plot.c <- grid_arrange_shared_legend(bias.cont, coverage.cont, rmse.c, width.c, nrow = 4, ncol = 1)


pdf("SMC MIC TeX documents/plot_Ycont.pdf")
grid_arrange_shared_legend(bias.cont, coverage.cont, width.c,rmse.c, nrow = 4, ncol = 1)
dev.off()


##### Bin results ####


results_bin<- read.csv("Results/bin_results_wc0.csv")
colnames(results_bin)

results_bin$Miss.Y<-factor(results_bin$Miss.Y, labels = c("fully obs", "MAR"))
results_bin$Miss<-interaction(results_bin$Miss.Y,results_bin$Sample.Size)

results_bin$Scenario<-interaction(results_bin$Sample.Size,results_bin$Noncompliance)
vlabel=c("N=200,30% Noncmp","N=200,50% Noncmp","N=1000,30% Noncmp","N=1000, 50% Noncmp")

results_bin$scenario<-factor(results_bin$Scenario, levels=c("200.30","200.50","1000.30", "1000.50"), labels=vlabel)


data2 <- results_bin %>% select(-(c("Outcome","Unobs.conf")))

data2<- data2 %>% 
  mutate(Noncompliance = as.factor(Noncompliance)) %>% 
mutate(Sample.Size = as.factor(Sample.Size)) %>%
  mutate(Effect.size = as.factor(Effect.size))%>%
  mutate(Miss.Y = as.factor(Miss.Y)) 

data_b_obs<- subset(data2, data2$Miss.Y==0)
data_b_Ymiss<- subset(data2, data2$Miss.Y==1)




#cbPalette <- c("#56B4E9", "#0072B2", "#D55E00", "#000000","darkred")
#Palette5 <- c("#56B4E9", "green", "#D55E00", "#000000","darkred")


#bias.bin.Yobs <-  ggplot(data_b_obs, aes(x=Method,y=bias, shape=Noncompliance, colour=Effect.size,)) +
#  geom_errorbar(aes(ymin=bias-1.96*mc.error, ymax=bias+1.96*mc.error), width=0.4, size=.4, position=position_dodge(0.5)) + 
#  geom_point(size=2,position=position_dodge(0.5))+
#  facet_grid(.~ Sample.Size, scales = "free") + 
#  theme_classic() +
#  geom_hline(yintercept= 0, linetype="dotted") +  
#  labs(  x= "Method", y = "Mean bias", shape="% noncompliance")



#all together

bias.b <-  ggplot(data2, aes(x=Method,y=bias, shape=Miss.Y, colour=Effect.size)) +
  geom_errorbar(aes(ymin=bias-1.96*mc.error, ymax=bias+1.96*mc.error), width=0.4, size=.4, position=position_dodge(0.75)) + 
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_hline(yintercept= 0, linetype="dotted") +  
  labs(  x= "", y = "Mean bias", shape="Binary outcome")




coverage.b <-
  ggplot(data2, aes(x=Method,y=coverage, shape=Miss.Y, colour=Effect.size, ymin=0.88,ymax=0.99)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") + 
  theme_classic() +
  geom_hline(yintercept= 0.95, linetype="dotted") +  
  geom_hline( aes(yintercept= round(.95-1.96*sqrt((.95*.05)/1999),3)), linetype="longdash", color="black", size=0.5)+ 
  geom_hline(aes(yintercept=round(.95+1.96*sqrt((.95*.05)/1999),3)), linetype="longdash", color="black", size=0.5)+ 
  labs(  x= "", y = "Coverage", shape="Binary outcome")


#RMSE
rmse.b <-
  ggplot(data=data2, aes(x=Method,y=RMSE, 
                         shape=Miss.Y, colour=Effect.size)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7.5), axis.title=element_text(size=8.2)) +
  xlab("Method") + ylab("RMSE") +
  theme(plot.title = element_text(hjust=0, size=8.2)) 


#CI width
width.b <-
  ggplot(data=data2, aes(x=Method,y=CI.width, 
                         shape=Miss.Y, colour=Effect.size)) +
  coord_cartesian(ylim = c(0, 16)) +
  geom_point(size=2,position=position_dodge(0.75))+
  scale_color_manual(values=Palette) + scale_shape_manual(values=c(15,22))+
  facet_grid(.~ scenario, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7.5), axis.title=element_text(size=8.2)) +
  xlab("Method") + ylab("CI width") +
  theme(plot.title = element_text(hjust=0, size=8.2)) 


grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom","right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


plot.b <- grid_arrange_shared_legend(bias.b, coverage.b, rmse.b, width.b, nrow = 4, ncol = 1)


pdf("SMC MIC TeX documents/plot_Ybin_wc0.pdf")
grid_arrange_shared_legend(bias.b, coverage.b,width.b, rmse.b,  nrow = 4, ncol = 1)
dev.off()


