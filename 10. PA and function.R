library(dplyr)
library(Hmisc)
library(vegan)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggpmisc)

a <- read.csv("Slope BEF data.csv")
var <- c("DOC","MBC","TDN","Olsen.P","MBN","MBP",
         "BG","NAG","LAP","AP","LeafC","LitterC")
pca <- summary(rda(dplyr::select(a,var),scale=T))
sites <- data.frame(pca$sites)
sites$PC1 <- -sites$PC1
PA <- read.csv("EMF-RF.csv")%>%dplyr::select(c("PNA_B_AMF","PPA_B_AMF","PNA_rpoB_AMF","PPA_rpoB_AMF"))
EMF <- read.csv("EMF-RF.csv")%>%dplyr::select(EMF)

data <- rcorr(as.matrix(cbind.data.frame(PA,EMF,sites[,1:3],dplyr::select(a,var))),type = "spearman")
r_value <- data$r
p_value <- data$P
r_related <- r_value[c(5:20),c(1:4)]
p_related <- p_value[c(5:20),c(1:4)]
r_related[p_related>0.05]=NA

PA_EMF <- cbind.data.frame(PA[,3:4],EMF,sites[,1],dplyr::select(a,c("DOC","MBC","NAG","AP")))%>%
  `colnames<-`(value=c("PNA_rpoB_AMF","PPA_rpoB_AMF","EMF","PC1","DOC","MBC","NAG","AP"))
PA_EMF <- data.frame(scale(PA_EMF,center = T,scale = T))
PA_EMF$Position <- c("Top","Top","Top","Top","Top","Top",
                     "Middle","Middle","Middle","Middle","Middle","Middle",
                     "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom")
PA_EMF$Position <- factor(PA_EMF$Position,levels = c("Top","Middle","Bottom"))

p1 <- ggplot(aes(x=PNA_rpoB_AMF, y=EMF),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4,label.x = "right")+
  theme(panel.grid = element_blank(),legend.position = c(0.1,0.2),
        legend.title = element_blank())

p2 <- ggplot(aes(x=PNA_rpoB_AMF, y=PC1),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4,label.x = "right")+ 
  theme(panel.grid = element_blank(),legend.position = "none")

p3 <- ggplot(aes(x=PNA_rpoB_AMF, y=DOC),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4,label.x = "right")+ 
  theme(panel.grid = element_blank(),legend.position = "none")

p4 <- ggplot(aes(x=PNA_rpoB_AMF, y=MBC),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4,label.x = "right")+ 
  theme(panel.grid = element_blank(),legend.position = "none")

p5 <- ggplot(aes(x=PNA_rpoB_AMF, y=NAG),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4,label.x = "right")+ 
  theme(panel.grid = element_blank(),legend.position = "none")

p6 <- ggplot(aes(x=PPA_rpoB_AMF, y=AP),data=PA_EMF) +
  geom_point(size=3,aes(color=Position))+
  #facet_wrap(~variable, scales="free") +
  theme_cowplot()+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_smooth(method="lm", colour="black",size=1.2,fill="lightgray") + 
  stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,size=4)+ 
  theme(panel.grid = element_blank(),legend.position = "none")

plot_grid( p1,p2,p3,p4,p5,p6,
           align = "hv",
           label_size = 15,
           labels = c("a","b","c","d","e","f"),
           nrow = 2,
           ncol=3
)