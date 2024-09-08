library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(grid)
library(ggpmisc)
library(lme4)
library(lmerTest)

make_pct <- function(x) (exp(x) - 1) * 100
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12))#11.6inches

a <- read.csv("Root biomass2.csv")
a$study <- as.factor(a$study)

root <- lmer(rrroot~1+(1|study),a)
root_boot <- confint(root,method="boot")

root_type <- lmer(rrroot~type-1+(1|study),a)
root_type_boot <- confint(root_type,method="boot")

root.df <- rbind.data.frame(root_boot[3,],root_type_boot[3:5,])
root.df$mean <- (root.df$`2.5 %`+root.df$`97.5 %`)/2
root.df <- make_pct(root.df)
root.df$type <- c("total","legume","mixture","non-legume")
root.df$n <- c(nrow(a),nrow(a%>%filter(type=="legume")),
               nrow(a%>%filter(type=="mixture")),nrow(a%>%filter(type=="non-legume")))

###Plot###
root.df$type <- factor(root.df$type,levels = c("total","legume","non-legume","mixture"),
                       labels = c("Total\n(152)","Legume\n(35)","Non-legume\n(41)","Mixture\n(76)"))
names(root.df) <- c("CI1","CI2","mean","type","n")
root.df$sig <- c("b","b","a","b")

ggplot(root.df)+
  geom_point(aes(x=mean,y=type,color=type),size=3.4)+
  geom_errorbar(aes(x=mean,xmin=CI1,xmax=CI2,y=type, col=type), width=0, size=0.75)+
  geom_vline(aes(xintercept=0),linetype="dashed",colour="black")+
  scale_color_manual(values = c("black","#a3d393","#f0d175","#f8984e"))+
  scale_y_discrete(limits=rev(levels(root.df$type)))+
  scale_x_continuous(limits=c(-15,45))+
  xlab("CC effect on root biomass (%)")+
  theme_cowplot()+
  theme(legend.position = "none",axis.title.y = element_blank())+
  annotate("text",label="n.s.",x=41,y=4)+
  annotate("text",label="n.s.",x=41,y=3)+
  annotate("text",label="n.s.",x=41,y=2)+
  annotate("text",label=expression(paste(italic(P)," < 0.05",sep="")),x=41,y=1)+font


