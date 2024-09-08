library(ggplot2)
library(ggpmisc)
library(ggplot2)
library(metafor)
library(boot)
library(dplyr)
library(lemon)
library(cowplot)

font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12))#11.6inches

a <- read.csv("SOC and clay interaction.csv")
a$ID <- as.factor(a$ID)
a$study <- as.factor(a$study)
a$clay <- log(a$clay)
a$SOC <- log(a$SOC)
###clay###
Yield.mod5 <- rma(rryield,varyield,data=a,mods=~clay)
SOC.mod5 <- rma(rrSOC,varSOC,data=a,mods=~clay)

Yield.Weight <- 1/(a$varyield+Yield.mod5$tau2)
SOC.Weight <- 1/(a$varSOC+SOC.mod5$tau2)
a <- cbind.data.frame(a,Yield.Weight,SOC.Weight)

p1 <- ggplot(a)+
  geom_point(aes(clay,rryield, size = Yield.Weight), color="#e4615d",alpha=0.2) +
  geom_smooth(aes(clay,rryield, weight = Yield.Weight), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), color = "black",se=T)+
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  #scale_x_continuous(limits = c(0,20))+
  #scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(clay,rryield,
                   label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 2, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  theme_cowplot()+
  xlab("Ln-Clay (%)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of yield",sep="")))+font
p1

p2 <- ggplot(a)+
  geom_point(aes(clay,rrSOC, size = SOC.Weight), color="#fdc58f",alpha=0.2) +
  geom_smooth(aes(clay,rrSOC, weight = SOC.Weight), method = "lm", 
              formula = y ~ poly(x, 1, raw = TRUE), color = "black",se=T)+
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  #scale_x_continuous(limits = c(0,20))+
  scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(clay,rrSOC,
                   label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  theme_cowplot()+
  xlab("Ln-Clay (%)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of ",SOC[stock],sep="")))+font
p2

###SOC###
Yield.mod6 <- rma(rryield,varyield,data=a,mods=~SOC)
SOC.mod6 <- rma(rrSOC,varSOC,data=a,mods=~SOC)

Yield.Weight2 <- 1/(a$varyield+Yield.mod6$tau2)
SOC.Weight2 <- 1/(a$varSOC+SOC.mod6$tau2)
a <- cbind.data.frame(a,Yield.Weight2,SOC.Weight2)

p3 <- ggplot(a)+
  geom_point(aes(SOC,rryield, size = Yield.Weight2), color="#e4615d",alpha=0.2) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  #scale_x_continuous(limits = c(0,20))+
  #scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(SOC,rryield,
                   label = paste(..p.value.label..)),
               formula = y ~ poly(x, 2, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  theme_cowplot()+
  xlab("Ln-Initial SOC (kg/ha)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of yield",sep="")))+font
p3

p4 <- ggplot(a)+
  geom_point(aes(SOC,rrSOC, size = SOC.Weight2), color="#fdc58f",alpha=0.2) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  #scale_x_continuous(limits = c(0,20))+
  scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(SOC,rrSOC,
                   label = paste(..p.value.label..)),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  theme_cowplot()+
  xlab("Ln-Initial SOC (kg/ha)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of ",SOC[stock],sep="")))+font
p4

###SOC-sand###
a$SOC_categorical <- factor(a$SOC_categorical,levels = c("low","moderate","high"),
                            labels = c("Low","Moderate","High"))

p5 <- ggplot(a)+
  geom_point(aes(clay,rryield, size = Yield.Weight), color="#e4615d",alpha=0.2) +
  geom_smooth(aes(clay,rryield, weight = Yield.Weight), method = "lm", 
              formula = y ~ poly(x, 1, raw = TRUE), color = "black",se=T)+
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  facet_wrap(~SOC_categorical,scales = "free_x")+
  #scale_x_continuous(limits = c(0,20))+
  #scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(clay,rryield,
                   label = paste(..eq.label..)),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  stat_poly_eq(aes(clay,rryield,
                   label = paste(..p.value.label..)),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "right",label.y = "bottom")+ 
  theme_cowplot()+
  xlab("Ln-Clay (%)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of yield",sep="")))+font
p5

p6 <- ggplot(a)+
  geom_point(aes(clay,rrSOC, size = SOC.Weight), color="#fdc58f",alpha=0.2) +
  geom_smooth(aes(clay,rrSOC, weight = SOC.Weight), method = "lm", 
              formula = y ~ poly(x, 1, raw = TRUE), color = "black",se=T)+
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") + 
  scale_size_continuous(range=c(1,3.6))+
  facet_wrap(~SOC_categorical,scales = "free_x")+
  #scale_x_continuous(limits = c(0,20))+
  scale_y_continuous(limits = c(-1.5,1.5))+
  stat_poly_eq(aes(clay,rrSOC,
                   label = paste(..eq.label..)),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "left",label.y = "top")+ 
  stat_poly_eq(aes(clay,rrSOC,
                   label = paste(..p.value.label..)),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "right",label.y = "bottom")+ 
  theme_cowplot()+
  xlab("Ln-Clay (%)")+theme(legend.position = "none")+
  ylab(expression(paste(Ln,italic(RR)," of ",SOC[stock],sep="")))+font
p6



ggdraw()+ 
  draw_plot(p1, x=0, y=0.55, width = 0.25, height = 0.45)+
  draw_plot(p2, x=0.25, y=0.55, width = 0.25, height = 0.45)+
  draw_plot(p3, x=0.5, y=0.55, width = 0.25, height = 0.45)+
  draw_plot(p4, x=0.75, y=0.55, width = 0.25, height = 0.45)+
  draw_plot(p5, x=0, y=0, width = 0.5, height = 0.55)+
  draw_plot(p6, x=0.5, y=0, width = 0.5, height = 0.55)+
  draw_plot_label(label = c("a","b","c","d","e","f"), size = 15,
                  x=c(0,0.25,0.5,0.75,0,0.5), y=c(1,1,1,1,0.55,0.55))##11.6*6.2inches
