library(vegan)
library(dplyr)
library(ape)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggrepel)

###PCA###
a <- read.csv("Slope BEF data.csv")
var <- c("DOC","MBC","TDN","Olsen.P","MBN","MBP",
         "BG","NAG","LAP","AP","LeafC","LitterC")
pca <- summary(rda(dplyr::select(a,var),scale=T))
sites <- data.frame(pca$sites)%>%mutate(Position=c("Top","Top","Top","Top","Top","Top",
                                       "Middle","Middle","Middle","Middle","Middle","Middle",
                                       "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"))
sites$Position <- factor(sites$Position,levels=c("Top","Middle","Bottom"))
species <- data.frame(pca$species)%>%mutate(func=rownames(pca$species))
species$func <- factor(species$func,levels = c("DOC","MBC","TDN","Olsen.P","MBN","MBP",
                                               "BG","NAG","LAP","AP","LeafC","LitterC"),
                       labels = c("DOC","MBC","TDN","Olsen P","MBN","MBP",
                                  "BG","NAG","LAP","AP","Leaf C","Litter C"))

adonis2(vegdist(dplyr::select(a,var),method="bray")~Position,data=a)

p1 <- ggplot() +
  geom_point(data=sites,aes(x=-PC1,y=PC2,fill=Position),size=3,color="transparent",shape=21) +
  geom_segment(data = species,x=0,y=0, aes(xend = -1.25*PC1, yend = 1.25*PC2),
               arrow = arrow(angle=22.5,length = unit(0.25,"cm"),
                             type = "closed")) +
  geom_text_repel(data = species,aes(-1.275*PC1,1.275*PC2,label=func),size=3.8)+
  geom_hline(yintercept = 0, linetype="dashed", color="grey") +
  geom_vline(xintercept = 0, linetype="dashed", color="grey") +
  scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  labs(x="PC 1",y="PC 2")+
  theme_classic2()+theme(legend.position = c(0.2,0.15),legend.title = element_blank(),
                              legend.spacing.x = unit(0.01,"lines"),
                              axis.line = element_line(colour = "transparent"),
                              axis.ticks = element_blank())
p1

###Multidimensional EMF###
PC1 <- sites[,c(1,7)]%>%`colnames<-`(value=c("PC","Position"))%>%
  mutate(PC=-PC,comp="Multifunctionality_Dim 1 (40.1%)")
PC2 <- sites[,c(2,7)]%>%`colnames<-`(value=c("PC","Position"))%>%
  mutate(comp="Multifunctionality_Dim 2 (16.3%)")
PC3 <- sites[,c(3,7)]%>%`colnames<-`(value=c("PC","Position"))%>%
  mutate(comp="Multifunctionality_Dim 3 (12.1%)")
multi <- rbind.data.frame(PC1,PC2,PC3)
multi$Position <- factor(multi$Position,levels = c("Top","Middle","Bottom"))

PC1$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_PC1 <- lmerTest::lmer(PC~Position+(1|ID),PC1)
anova(model_PC1)##F=34.2,P<0.001
emmeans::emmeans(model_PC1,pairwise~Position)
PC2$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_PC2 <- lmerTest::lmer(PC~Position+(1|ID),PC2)
anova(model_PC2)##F=2.26,P>0.05
emmeans::emmeans(model_PC2,pairwise~Position)
PC3$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_PC3 <- lmerTest::lmer(PC~Position+(1|ID),PC3)
anova(model_PC3)##F=0.40,P>0.05
emmeans::emmeans(model_PC3,pairwise~Position)

p2 <- ggplot(multi, aes(x=Position, y=PC)) + 
  geom_bar(position = "identity", stat = "summary",aes(fill=Position),width = 0.75) + 
  geom_errorbar(stat="summary", width=0.25)+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  theme_pubclean() + ylab("Multifunctionality") +
  facet_wrap(~comp, ncol = 4)+
  theme(axis.line = element_line(), axis.title.x = element_blank(),legend.position = "none")
p2

###Average EMF###
b <- read.csv("Fig. 1.csv")
b$Position <- factor(b$Position,levels = c("Top","Mid","Bottom"),
                     labels = c("Top","Middle","Bottom"))
b$var <- "Average multifunctionality"

b$ID <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
model_b <- lmerTest::lmer(EMF~Position+(1|ID),b)
anova(model_b)##F=16.6,P<0.001
emmeans::emmeans(model_b,pairwise~Position)

p3 <- ggplot(b, aes(x=Position, y=EMF)) + 
  geom_bar(position = "identity", stat = "summary",aes(fill=Position),width = 0.75) + 
  geom_errorbar(stat="summary", width=0.25)+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  facet_wrap(~var)+
  theme_pubclean() + ylab("Multifunctionality") +
  theme(axis.line = element_line(), axis.title.x = element_blank(),legend.position = "none")
p3

p4 <- plot_grid( p1,p3,
           align = 'h', 
           label_size = 15,
           labels = c("a","b"),
           nrow = 1,
           ncol=2,
           rel_widths = c(1,0.75)
)

ggdraw()+ 
  draw_plot(p4, x=0, y=0.5, width = 1, height = 0.5)+
  draw_plot(p2, x=0, y=0, width = 1, height = 0.5)+
  draw_plot_label(label=c("c"),size = 15,
                  x=c(0),y=c(0.5))+
  draw_label("PERMANOVA\nF = 7.39, P < 0.001",
             x=0.375,y=0.605,size=12,hjust = 0)
