library(phyloseq)
library(ape)
library(qcmi)
library(ggtree)
library(Hmisc)
library(gplots)
library(reshape2)
library(ggsci)
library(dplyr)
library(cowplot)
library(vegan)
library(ggvenn)
library(linkET)

env <- read.csv("environment.csv",row.names = 1)
PA <- read.csv("EMF-RF.csv")%>%select(c("PNA_B_AMF","PPA_B_AMF","PNA_rpoB_AMF","PPA_rpoB_AMF"))
AMF.r <- read.csv("EMF-RF.csv")%>%select(value="Richness_AMF")
env2 <- cbind.data.frame(env,PA,AMF.r)%>%`colnames<-`(value=c("Moisture","pH","BD","SOC","N:P","TP",
                                                        "PNA_B_AMF","PPA_B_AMF",
                                                        "PNA_rpoB_AMF","PPA_rpoB_AMF","AMF richness"))

#####Bacteria#####
asv_B <- read.csv("ASV-B_qcmi.csv",row.names = 1)
group_B <- read.csv("group.csv",row.names = 1)
tax_B <- read.csv("B_node_new.csv",row.names = 1)
tax_B$Order <- factor(tax_B$Order,
           levels = c("Rhizobiales","Chitinophagales","Sphingomonadales",
                      "Pyrinomonadales","Gemmatimonadales","Others"),
           labels = c("1Rhizobiales","2Chitinophagales","3Sphingomonadales",
                      "4Pyrinomonadales","5Gemmatimonadales","6Others"))
ps_B <- trans_ps(otu_table=asv_B,taxa_table=tax_B,sample_data = sample_data(group_B))
tree_B <- rtree(ntaxa(ps_B),rooted = T,tip.label = taxa_names(ps_B))
dd<-data.frame(tax_B[tree_B$tip.label,])
bb<-tax_B[tree_B$tip.label,]
groupInfo1 <- split(tree_B$tip.label, as.factor(bb$Order))
aa<- groupOTU(tree_B, groupInfo1)
d1<-data.frame(ASV=dd$name,lab=dd$name)
write.csv(d1,"label_B.csv")
d1 <- read.csv("label_B.csv",row.names = 1)

data <- rcorr(as.matrix(env2[,c(1:8,11)]),t(as.matrix(asv_B)),type="spearman")
r_value <- data$r
p_value <- data$P
r_related <- r_value[c(1:9),c(10:424)]
p_related <- p_value[c(1:9),c(10:424)]
r_related[p_related>0.05]=NA
heat<-t(r_related)

p<-ggtree(aa,aes(color=group),layout ="fan",branch.length='none',open.angle=90,lwd=0.5)
p<-p %<+% d1+geom_tippoint(size=0.0001)+
  geom_tiplab2(aes(label=lab),offset=17.5,align=TRUE,linetype=NA,size=0.8,linesize=0)+
  scale_color_nejm()

p1 <- gheatmap(p,heat,offset=0, width=1, font.size=3.4,hjust=0,color="lightgray",
         colnames_angle=0,colnames = T,colnames_position = "top")+
  scale_fill_gradientn(colours=colorpanel(100,low="#E53935",high="#1E88E5",mid = "white"),
                       na.value="white")+
  scale_size_continuous(range=c(0,1.5))+
  theme(legend.position = c(0.825,0.3),legend.box = "horizontal")+
  guides(colour = guide_legend(title = "", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's r", order = 1,title.position = "top",
                               direction = "horizontal",barheight = unit(0.5,"lines"),
                               barwidth = unit(5,"lines")))
  
#####rpoB#####
asv_rpoB <- read.csv("ASV-rpoB_qcmi.csv",row.names = 1)
group_rpoB <- read.csv("group.csv",row.names = 1)
tax_rpoB <- read.csv("rpoB_node_new.csv",row.names = 1)
tax_rpoB$Genus <- factor(tax_rpoB$Genus,
                      levels = c("Bradyrhizobium","Afipia","Brucella",
                                 "Bartonella","Mesorhizobium","Others"),
                      labels = c("1Bradyrhizobium","2Afipia","3Brucella",
                                 "4Bartonella","5Mesorhizobium","6Others"))
ps_rpoB <- trans_ps(otu_table=asv_rpoB,taxa_table=tax_rpoB,sample_data = sample_data(group_rpoB))
tree_rpoB <- rtree(ntaxa(ps_rpoB),rooted = T,tip.label = taxa_names(ps_rpoB))
dd2<-data.frame(tax_rpoB[tree_rpoB$tip.label,])
bb2<-tax_rpoB[tree_rpoB$tip.label,]
groupInfo2 <- split(tree_rpoB$tip.label, as.factor(bb2$Genus))
aa2<- groupOTU(tree_rpoB, groupInfo2)
d2<-data.frame(ASV=dd2$name,lab=dd2$name)
write.csv(d2,"label_rpoB.csv")
d2 <- read.csv("label_rpoB.csv",row.names = 1)

data2 <- rcorr(as.matrix(env2[,c(1:6,9:11)]),t(as.matrix(asv_rpoB)),type="spearman")
r_value2 <- data2$r
p_value2 <- data2$P
r_related2 <- r_value2[c(1:9),c(10:130)]
p_related2 <- p_value2[c(1:9),c(10:130)]
r_related2[p_related2>0.05]=NA
heat2<-t(r_related2)

p2<-ggtree(aa2,aes(color=group),layout ="fan",branch.length='none',open.angle=90,lwd=0.5)
p2<-p2 %<+% d2+geom_tippoint(size=0.0001)+
  geom_tiplab2(aes(label=lab),offset=14.25,align=TRUE,linetype=NA,size=1.5,linesize=0)+
  scale_color_nejm()

p3 <- gheatmap(p2,heat2,offset=0, width=1, font.size=3.4,hjust=0,color="lightgray",
               colnames_angle=0,colnames = T,colnames_position = "top")+
  scale_fill_gradientn(colours=colorpanel(100,low="#E53935",high="#1E88E5",mid = "white"),
                       na.value="white")+
  scale_size_continuous(range=c(0,1.5))+
  theme(legend.position = c(0.8,0.3),legend.box = "horizontal")+
  guides(colour = guide_legend(title = "", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's r", order = 1,title.position = "top",
                               direction = "horizontal",barheight = unit(0.5,"lines"),
                               barwidth = unit(5,"lines")))

p4 <- plot_grid( p1+theme(plot.background = element_blank()),p3+theme(plot.background = element_blank()),
           align = 'hv', 
           label_size = 15,
           labels = c("a","b"),
           nrow = 1,
           ncol=2
)

###VPA###
vpa_B <- cbind.data.frame(t(asv_B),env2[,c(1:8,11)])
mod_B <- varpart(vegdist(vpa_B[,1:415],method="bray"), ~ Moisture+pH+BD+SOC+`N:P`+TP, 
               ~ PNA_B_AMF+PPA_B_AMF+`AMF richness`, data=vpa_B)
plot(mod_B,bg=c("steelblue","darkred"),col=c("steelblue","darkred"),
     Xnames=c("Abiotic","Biotic"))+
  title(main="Explained variation for bacterial composition")
abiotic_B <- model.matrix(~ Moisture+pH+BD+SOC+`N:P`+TP, data=vpa_B)
aFrac_B <- dbrda(vegdist(vpa_B[,1:415],method="bray")~abiotic_B, data=vpa_B)
anova(aFrac_B)
biotic_B <- model.matrix(~ PNA_B_AMF+PPA_B_AMF+`AMF richness`, data=vpa_B)
bFrac_B <- dbrda(vegdist(vpa_B[,1:415],method="bray")~biotic_B, data=vpa_B)
anova(bFrac_B)

vpa_rpoB <- cbind.data.frame(t(asv_rpoB),env2[,c(1:6,9:11)])
mod_rpoB <- varpart(vegdist(vpa_rpoB[,1:121],method="bray"), ~ Moisture+pH+BD+SOC+`N:P`+TP, 
                 ~ PNA_rpoB_AMF+PPA_rpoB_AMF+`AMF richness`, data=vpa_rpoB)
plot(mod_rpoB)
abiotic_rpoB <- model.matrix(~ Moisture+pH+BD+SOC+`N:P`+TP, data=vpa_rpoB)
aFrac_rpoB <- dbrda(vegdist(vpa_rpoB[,1:121],method="bray")~abiotic_rpoB, data=vpa_rpoB)
anova(aFrac_rpoB)
biotic_rpoB <- model.matrix(~ PNA_rpoB_AMF+PPA_rpoB_AMF+`AMF richness`, data=vpa_rpoB)
bFrac_rpoB <- dbrda(vegdist(vpa_rpoB[,1:121],method="bray")~biotic_rpoB, data=vpa_rpoB)
anova(bFrac_rpoB)

sets <- list(
  Abiotic = sample(1:500,300),
  Biotic = sample(1:500,350)
)

p5 <- ggvenn(sets,fill_color = c("#1E88E5","#E53935"),stroke_color = "transparent",
       set_name_color = c("#1E88E5","#E53935"),set_name_size = 4.6,text_size = 4.2)+
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.background = element_blank())

###Mantel test###
mantel <- cbind.data.frame(t(asv_B),t(asv_rpoB),env2)
mt <- mantel_test(mantel[,1:536],mantel[,537:547],
                  spec_select = list(`Bacterial\ncomposition\n(16S ASVs)`=1:415,
                                     `Rhizobial\ncomposition\n(rpoB ASVs)`=416:536))%>%
  mutate(rd=cut(r,breaks = c(-Inf,0.2,0.4,Inf),labels = c("< 0.2","0.2 - 0.4",">= 0.4")),
         pd=cut(p,breaks = c(-Inf,0.01,0.05,Inf),labels = c("< 0.01","0.01 - 0.05",">= 0.05")))

p7 <- qcorrplot(correlate(mantel[,537:547],method = "spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = mt, 
              curvature = nice_curvature(),
              offset_x = list(`Bacterial\ncomposition\n(16S ASVs)` = -1, 
                              `Rhizobial\ncomposition\n(rpoB ASVs)` = -0.75),
              offset_y = list(`Bacterial\ncomposition\n(16S ASVs)` = 1.25,
                              `Rhizobial\ncomposition\n(rpoB ASVs)` = 0.75)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  scale_fill_gradient2(low = "#E53935",mid="white",high = "#1E88E5")+
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.title = element_blank(),legend.position = "right",plot.background = element_blank(),
        legend.background = element_blank())

ggdraw()+ 
  draw_plot(p4, x=0, y=0.5, width = 1, height = 0.5)+
  draw_plot(p7, x=0, y=0, width = 0.6, height = 0.5)+
  draw_plot(p5, x=0.6, y=0.25, width = 0.4, height = 0.245)+
  draw_plot(p5, x=0.6, y=0, width = 0.4, height = 0.245)+
  draw_plot_label(label=c("c","d"),size = 15,
                  x=c(0,0.6),y=c(0.5,0.5))
