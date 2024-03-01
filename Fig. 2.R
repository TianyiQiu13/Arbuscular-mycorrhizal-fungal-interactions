library(ggplot2)
library(cowplot)
library(lme4)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)

a <- read.csv("Alpha.csv")
a$Position <- factor(a$Position,levels=c("Top","Mid","Bottom"),
                     labels = c("Top","Middle","Bottom"))

bacteria <- a%>%filter(taxa=="Bacterial 16S")
bacteria$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_B <- lmerTest::lmer(Richness~Position+(1|ID),bacteria)
anova(model_B)
emmeans::emmeans(model_B,pairwise~Position)

rhizobia <- a%>%filter(taxa=="Rhizobial rpoB")
rhizobia$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_rpoB <- lmerTest::lmer(Richness~Position+(1|ID),rhizobia)
anova(model_rpoB)
emmeans::emmeans(model_rpoB,pairwise~Position)

AMF <- a%>%filter(taxa=="AM fungal 18S")
AMF$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_AMF <- lmerTest::lmer(Richness~Position+(1|ID),AMF)
anova(model_AMF)
emmeans::emmeans(model_AMF,pairwise~Position)

p1 <- ggplot(bacteria,aes(x=Position,y=Richness,color=Position))+
  geom_boxplot(outlier.color = NA,width = 0.65)+
  geom_jitter(size=1.8)+
  facet_wrap(.~taxa,scales = "free_y")+
  scale_y_continuous(limits = c(3000,4700))+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  annotate("text",x=1,y=6.4/7*1700+3000,
            label="F = 2.12\nP > 0.05",size=4.2)+
  theme_test()+
  theme(legend.position = "none",axis.title.x = element_blank())
p1

p2 <- ggplot(rhizobia,aes(x=Position,y=Richness,color=Position))+
  geom_boxplot(outlier.color = NA,width = 0.65)+
  geom_jitter(size=1.8)+
  facet_wrap(.~taxa,scales = "free_y")+
  scale_y_continuous(limits = c(400,700))+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  annotate("text",x=1,y=6.4/7*300+400,
           label="F = 0.27\nP > 0.05",size=4.2)+
  theme_test()+
  theme(legend.position = "none",axis.title.x = element_blank())
p2

p3 <- ggplot(AMF,aes(x=Position,y=Richness,color=Position))+
  geom_boxplot(outlier.color = NA,width = 0.65)+
  geom_jitter(size=1.8)+
  facet_wrap(.~taxa,scales = "free_y")+ 
  scale_y_continuous(limits = c(100,800))+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  annotate("text",x=1,y=800*0.925,
           label="F = 0.88\nP > 0.05",size=4.2)+
  theme_test()+
  theme(legend.position = "none",axis.title.x = element_blank())
p3

group <- read.csv("group.csv",row.names = 1)
group$treats <- factor(group$treats,levels = c("Top","Mid","Bottom"))

ASV_B <- read.csv("ASV-B.csv",row.names = 1)
p4 <- amplicon::beta_cpcoa(ASV_B,group,groupID = "treats")+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  theme(legend.position = "none",axis.title.x = element_blank())

ASV_rpoB <- read.csv("ASV-rpoB.csv",row.names = 1)
p5 <- amplicon::beta_cpcoa(ASV_rpoB,group,groupID = "treats")+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  theme(legend.position = "none",axis.title.x = element_blank())

ASV_AMF <- read.csv("ASV-AMF.csv",row.names = 1)
p6 <- amplicon::beta_cpcoa(ASV_AMF,group,groupID = "treats")+
  scale_color_manual(values = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  theme(legend.position = "none",axis.title.x = element_blank())

c <- read.csv("RF1-plot.csv")
c$variable <- factor(c$variable,levels=c("N.P","CAP1_B","SOC","CAP1_rpoB",
                                         "Position","Richness_AMF","Moisture",
                                         "CAP1_AMF","TP","BD","pH","Richness_rpoB","Richness_B"),
                     labels = c("N:P","CPCoA1-Bacteria","SOC","CPCoA1-Rhizobia","Slope position",
                                "Richness-AMF","Moisture","CPCoA1-AMF","TP","BD","pH",
                                "Richness-Rhizobia","Richness-Bacteria"))

p7 <- ggplot()+
  geom_bar(data = c, 
           aes(x= reorder(variable, -IncMSE), 
               y=IncMSE, fill=color),
           stat = "identity")+
  labs(y="%IncMSE\n") +
  scale_fill_manual(values = c("steelblue","gray"))+
  geom_text(data = c,aes(y=IncMSE+0.35,x=variable,label=label),size=4.6)+
  theme_cowplot()+
  annotate("text",x=11,y=8,
           label="Slope multifunctionality\nR2 = 0.67, P < 0.01",size=4.6)+
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1))
p7

p8 <- plot_grid( p1,p2,p3, p4,p5,p6,
           align = 'v', 
           label_size = 15,
           labels = c("a","b","c","d","e","f"),
           nrow = 2,
           ncol=3,
           rel_heights = c(1,1.25)
)

p9 <- plot_grid(p7,
          label_size = 15,
          labels = c("g"),
          nrow = 1,
          ncol=1)

ggdraw()+ 
  draw_plot(p8, x=0, y=1/3, width = 1, height = 2/3)+
  draw_plot(p9, x=0, y=0, width = 1, height = 1/3)
