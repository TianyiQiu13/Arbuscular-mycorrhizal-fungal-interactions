library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)

zc_B <- read.csv("B_AMF_microeco_node.csv",row.names = 1)%>%mutate(type="Bacteria + AMF")
zc_rpoB <- read.csv("rpoB_AMF_microeco_node.csv",row.names = 1)%>%mutate(type="Rhizobia + AMF")

zc <- rbind.data.frame(zc_B[,c(6:9,17)],zc_rpoB[,c(6:9,17)])
zc$filed <- factor(zc$filed,levels = c("B","rpoB","AMF"),
                   labels = c("Bacteria","Rhizobia","AMF"))
zc$taxa_roles <- factor(zc$taxa_roles,levels = c("Connectors","Network hubs","Peripheral nodes"),
                        labels = c("Connectors","Network hubs","Peripherals"))

ggplot(zc,aes(x=p,y=z))+
  geom_point(aes(color=filed,shape=taxa_roles))+
  geom_hline(yintercept = 2.5,linetype=2)+
  geom_vline(xintercept = 0.62,linetype=2)+
  scale_color_manual(values = c("#0072B5FF","#BC3C29FF","#FFDC91FF"))+
  labs(x="Among-module connectivity (P)",y="Within-module connectivity (Z)",
       color="Kingdom",shape="Taxa roles")+
  facet_wrap(~type)+
  theme_cowplot()+
  guides(color=guide_legend(order = 1),shape=guide_legend(order = 2))+
  annotate("text",x=0.73,y=-1.1,label = "AMF_ASV186",size=3.8)+
  annotate("text",x=0.63,y=-1.4,label = "AMF_ASV176",size=3.8)