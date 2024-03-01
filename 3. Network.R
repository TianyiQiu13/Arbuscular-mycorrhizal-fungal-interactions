library(ggplot2)
library(igraph)
library(ggnetwork)
library(cowplot)
library(lemon)
library(ggsci)
library(dplyr)
library(ggpubr)

net_B <- read_graph("B_microeco.graphml",format = "graphml")
net_B_AMF <- read_graph("B_AMF_microeco.graphml",format = "graphml")
net_rpoB <- read_graph("rpoB_microeco.graphml",format = "graphml")
net_rpoB_AMF <- read_graph("rpoB_AMF_microeco.graphml",format = "graphml")

net_B <- set_vertex_attr(net_B,"hub.score",V(net_B),hub_score(net_B)$vector)
net_B_AMF <- set_vertex_attr(net_B_AMF,"hub.score",V(net_B_AMF),hub_score(net_B_AMF)$vector)
net_rpoB <- set_vertex_attr(net_rpoB,"hub.score",V(net_rpoB),hub_score(net_rpoB)$vector)
net_rpoB_AMF <- set_vertex_attr(net_rpoB_AMF,"hub.score",V(net_rpoB_AMF),hub_score(net_rpoB_AMF)$vector)

order_B <- read.csv("B_node_new.csv",row.names = 1)
order_B_AMF <- read.csv("B_AMF_node_new.csv",row.names = 1)
genus_rpoB <- read.csv("rpoB_node_new.csv",row.names = 1)
genus_rpoB_AMF <- read.csv("rpoB_AMF_node_new.csv",row.names = 1)

net_B <- set_vertex_attr(net_B,"classification",V(net_B),order_B$Order)
net_B_AMF <- set_vertex_attr(net_B_AMF,"classification",V(net_B_AMF),order_B_AMF$Order)
net_rpoB <- set_vertex_attr(net_rpoB,"classification",V(net_rpoB),genus_rpoB$Genus)
net_rpoB_AMF <- set_vertex_attr(net_rpoB_AMF,"classification",V(net_rpoB_AMF),genus_rpoB_AMF$Genus)

###layout###
ggnet_B <- ggnetwork(net_B,"kamadakawai")%>%mutate(type="Bacteria")
ggnet_B_AMF <- ggnetwork(net_B_AMF,"kamadakawai")%>%mutate(type="Bacteria + AMF")
ggnet_rpoB <- ggnetwork(net_rpoB,"kamadakawai")%>%mutate(type="Rhizobia")
ggnet_rpoB_AMF <- ggnetwork(net_rpoB_AMF,"kamadakawai")%>%mutate(type="Rhizobia + AMF")

networks_B <- rbind.data.frame(ggnet_B,ggnet_B_AMF)
networks_rpoB <- rbind.data.frame(ggnet_rpoB,ggnet_rpoB_AMF)

###plot###
networks_B$classification <- factor(networks_B$classification,
                                    levels = c("Rhizobiales","Chitinophagales","Sphingomonadales",
                                               "Pyrinomonadales","Gemmatimonadales","Others","AMF"))
networks_rpoB$classification <- factor(networks_rpoB$classification,
                                    levels = c("Bradyrhizobium","Afipia","Brucella",
                                               "Bartonella","Mesorhizobium","Others","AMF"))

p1 <- ggplot(networks_B, aes(x = x,  y = y, xend = xend,  yend = yend))+
  geom_edges(size=0.5, alpha = 0.5, curvature = 0.25, aes(colour=label),show.legend = F) +
  geom_nodes(aes(size=hub.score,fill=classification), shape=21,color="#3D405B")+
  geom_nodetext_repel(data=filter(networks_B,hub.score>0.7),
                      label=c("B_ASV45","B_ASV78",
                              "AMF_ASV186","AMF_ASV176"),
                      size=3.8,color="#212d8e")+
  scale_size_continuous(guide = "none")+
  scale_color_manual(values=c("#D57A66","#A3CBB2"))+
  scale_fill_nejm()+
  theme_blank()+
  coord_equal()+
  facet_wrap(~type, ncol = 2)+
  theme(legend.position = "bottom",legend.box.margin = margin(l=0,t=-1.5,b=-0.5,unit="lines"),
        panel.spacing = unit(2,"lines"),
        strip.text = element_text(size=13,face = "bold"),strip.background = element_blank(),
        legend.title = element_blank(),legend.direction = "horizontal",
        legend.text = element_text(size = 12),legend.spacing.x = unit(0.01,"lines"))

p2 <- ggplot(networks_rpoB, aes(x = x,  y = y, xend = xend,  yend = yend))+
  geom_edges(size=0.5, alpha = 0.5, curvature = 0.25, aes(colour=label),show.legend = F) +
  geom_nodes(aes(size=hub.score,fill=classification), shape=21,color="#3D405B")+
  geom_nodetext_repel(data=filter(networks_rpoB,hub.score>0.7),
                      label=c("rpoB_ASV102","AMF_ASV186",
                              "AMF_ASV176"),size=3.8,color="#212d8e")+
  labs(fill="")+
  scale_size_continuous(guide = "none")+
  scale_color_manual(values=c("#D57A66","#A3CBB2"))+
  scale_fill_nejm()+
  theme_blank()+
  coord_equal()+
  facet_wrap(~type, ncol = 2)+
  theme(legend.position = "bottom",legend.box.margin = margin(l=0,t=-1.5,b=-0.5,unit="lines"),
        panel.spacing = unit(2,"lines"),
        strip.text = element_text(size=13,face = "bold"),strip.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),legend.spacing.x = unit(0.01,"lines"))

p3 <- plot_grid( p1,p2,
                 align = 'hv', 
                 label_size = 15,
                 labels = c("a","b"),
                 #hjust = -2, 
                 #vjust= 4,
                 nrow = 2,
                 ncol=1,
                 rel_heights = c(1,1)
)
p3

###biotic association###
PA_B <- read.csv("B_PA2.csv",row.names = 1)%>%
  `colnames<-`(value=c("PNA","PPA"))%>%mutate(type="Bacteria")
PA_B_AMF <- read.csv("B_AMF_PA2.csv",row.names = 1)%>%
  `colnames<-`(value=c("PNA","PPA"))%>%mutate(type="Bacteria + AMF")
PA_rpoB <- read.csv("rpoB_PA2.csv",row.names = 1)%>%
  `colnames<-`(value=c("PNA","PPA"))%>%mutate(type="Rhizobia")
PA_rpoB_AMF <- read.csv("rpoB_AMF_PA2.csv",row.names = 1)%>%
  `colnames<-`(value=c("PNA","PPA"))%>%mutate(type="Rhizobia + AMF")

PAs_B <- rbind.data.frame(PA_B,PA_B_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Bacteria vs. Bacteria + AMF")
PAs_rpoB <- rbind.data.frame(PA_rpoB,PA_rpoB_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Rhizobia vs. Rhizobia + AMF")
PAs_all <- rbind.data.frame(PAs_B,PAs_rpoB)
PAs_all$Position <- factor(PAs_all$Position,levels = c("Top","Middle","Bottom"))

wilcox.test(PNA~type,data = filter(PAs_B,Position=="Top"))

p4 <- ggplot(PAs_all, aes(x=Position, y=abs(PNA))) + 
  geom_boxplot(aes(color=type),outlier.colour = NA) +
  geom_jitter(aes(color=type), alpha=0.2, size=2, shape=16,width = 0.15) +
  labs(color="",y="|Putative negative biotic associations|")+
  theme_pubclean()+
  theme(axis.line = element_line(),axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.key = element_rect(fill="transparent"),legend.box.spacing = unit(-1,"lines"),
        strip.background = element_blank(),strip.text = element_text(color=NA))+
  facet_grid(~var)

p5 <- ggplot(PAs_all, aes(x=Position, y=PPA)) + 
  geom_boxplot(aes(color=type),outlier.colour = NA) +
  geom_jitter(aes(color=type), alpha=0.2, size=2, shape=16,width = 0.15) +
  labs(color="",y="Putative positive biotic associations")+
  theme_pubclean()+
  theme(axis.line = element_line(),axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.key = element_rect(fill="transparent"),legend.box.spacing = unit(-1,"lines"),
        strip.background = element_blank(),strip.text = element_text(color=NA))+
  facet_grid(~var)

###centrality###
cen_B <- read.csv("centrality_sub_B.csv",row.names = 1)%>%
  `colnames<-`(value="centrality")%>%mutate(type="Bacteria")
cen_B_AMF <- read.csv("centrality_sub_B_AMF.csv",row.names = 1)%>%
  `colnames<-`(value="centrality")%>%mutate(type="Bacteria + AMF")
cen_rpoB <- read.csv("centrality_sub_rpoB.csv",row.names = 1)%>%
  `colnames<-`(value="centrality")%>%mutate(type="Rhizobia")
cen_rpoB_AMF <- read.csv("centrality_sub_rpoB_AMF.csv",row.names = 1)%>%
  `colnames<-`(value="centrality")%>%mutate(type="Rhizobia + AMF")

cens_B <- rbind.data.frame(cen_B,cen_B_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Bacteria vs. Bacteria + AMF")
cens_rpoB <- rbind.data.frame(cen_rpoB,cen_rpoB_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Rhizobia vs. Rhizobia + AMF")
cens_all <- rbind.data.frame(cens_B,cens_rpoB)
cens_all$Position <- factor(cens_all$Position,levels = c("Top","Middle","Bottom"))

wilcox.test(centrality~Position,data = filter(cens_B,Position!="Top"))

p6 <- ggplot(cens_all, aes(x=Position, y=centrality)) + 
  geom_boxplot(aes(color=type),outlier.colour = NA) +
  geom_jitter(aes(color=type), alpha=0.2, size=2, shape=16,width = 0.15) +
  labs(color="",y="Kleinberg's hub centrality")+
  theme_pubclean()+
  theme(axis.line = element_line(),axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.key = element_rect(fill="transparent"),legend.box.spacing = unit(-1,"lines"),
        strip.background = element_blank(),strip.text = element_text(color=NA))+
  facet_grid(~var)

###complexity###
comp_B <- read.csv("net_sub_B.csv",row.names = 1)%>%select(density)%>%
  mutate(type="Bacteria")
comp_B_AMF <- read.csv("net_sub_B_AMF.csv",row.names = 1)%>%select(density)%>%
  mutate(type="Bacteria + AMF")
comp_rpoB <- read.csv("net_sub_rpoB.csv",row.names = 1)%>%select(density)%>%
  mutate(type="Rhizobia")
comp_rpoB_AMF <- read.csv("net_sub_rpoB_AMF.csv",row.names = 1)%>%select(density)%>%
  mutate(type="Rhizobia + AMF")

comps_B <- rbind.data.frame(comp_B,comp_B_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Bacteria vs. Bacteria + AMF")
comps_rpoB <- rbind.data.frame(comp_rpoB,comp_rpoB_AMF)%>%
  mutate(Position=c("Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom",
                    "Top","Top","Top","Top","Top","Top",
                    "Middle","Middle","Middle","Middle","Middle","Middle",
                    "Bottom","Bottom","Bottom","Bottom","Bottom","Bottom"),
         var="Rhizobia vs. Rhizobia + AMF")
comps_all <- rbind.data.frame(comps_B,comps_rpoB)
comps_all$Position <- factor(comps_all$Position,levels = c("Top","Middle","Bottom"))

wilcox.test(density~type,data = filter(comps_rpoB,Position=="Top"))

p7 <- ggplot(comps_all, aes(x=Position, y=density*100)) + 
  geom_boxplot(aes(color=type),outlier.colour = NA) +
  geom_jitter(aes(color=type), alpha=0.2, size=2, shape=16,width = 0.15) +
  labs(color="",y="Network density (10-2)")+
  theme_pubclean()+
  #scale_y_continuous(breaks = c(0,1.2,2.4),labels = c(0,1.2,2.4))+
  theme(axis.line = element_line(),axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.key = element_rect(fill="transparent"),legend.box.spacing = unit(-1,"lines"),
        strip.background = element_blank(),strip.text = element_text(color=NA))+
  facet_grid(~var)


p8 <- lemon::grid_arrange_shared_legend(p6,p7,p4,p5,nrow = 2,ncol = 2,position = "top")
p9 <- plot_grid( p8,
                 nrow = 1,
                 ncol=1)
p9

ggdraw()+ 
  draw_plot(p3, x=0, y=0, width = 0.55, height = 1)+
  draw_plot(p9, x=0.54, y=0, width = 0.46, height = 1)+
  draw_plot_label(label=c("c","d","e","f"),size = 15,
                  x=c(0.53,0.765,0.53,0.765),y=c(0.93,0.93,0.46,0.46))+
  draw_label("positive: 1110\nnegative: 559\nmodularity: 0.63",
               x=0.225,y=0.92,size=10,hjust = 0)+
  draw_label("positive: 2586\nnegative: 635\nmodularity: 0.60",
             x=0.305,y=0.625,size=10,hjust = 0)+
  draw_label("positive: 74\nnegative: 29\nmodularity: 0.86",
             x=0.225,y=0.42,size=10,hjust = 0)+
  draw_label("positive: 733\nnegative: 80\nmodularity: 0.77",
             x=0.435,y=0.125,size=10,hjust = 0)+
  draw_label("*",x=0.65,y=0.875,size=12)+
  draw_label("**",x=0.705,y=0.875,size=12)+
  draw_label("***",x=0.7525,y=0.875,size=12)+
  draw_label("***",x=0.93,y=0.875,size=12)+
  draw_label("***",x=0.955,y=0.875,size=12)+
  draw_label("***",x=0.98,y=0.875,size=12)+
  draw_label("***",x=0.62,y=0.4,size=12)+
  draw_label("***",x=0.645,y=0.4,size=12)+
  draw_label("***",x=0.67,y=0.4,size=12)+
  draw_label("***",x=0.845,y=0.4,size=12)+
  draw_label("**",x=0.87,y=0.4,size=12)+
  draw_label("***",x=0.93,y=0.4,size=12)+
  draw_label("***",x=0.955,y=0.4,size=12)+
  draw_label("***",x=0.98,y=0.4,size=12)