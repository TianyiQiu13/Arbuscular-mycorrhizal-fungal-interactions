library(iCAMP)
library(phyloseq)
library(qcmi)
library(ggpmisc)
library(ggtree)
library(dplyr)
library(cowplot)
library(vegan)
library(tidyverse)
library(ggsci)

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
comm_B <- t(asv_B)
pd_B <- cophenetic(tree_B)

set.seed(123)
qpen.out_B <- qpen(comm = comm_B, pd = pd_B, sig.bNTI = 2, sig.rc = 0.95, rand.time = 1000, nworker = 20)

env <- read.csv("environment.csv",row.names = 1)
PA <- read.csv("EMF-RF.csv")%>%select(c("PNA_B_AMF","PPA_B_AMF","PNA_rpoB_AMF","PPA_rpoB_AMF"))
AMF.r <- read.csv("AMF richness.csv",row.names = 1)%>%`colnames<-`(value="AMF richness")
env2 <- cbind.data.frame(env,PA,AMF.r)%>%`colnames<-`(value=c("Moisture","pH","BD","SOC","N:P","TP",
                                                              "PNA_B_AMF","PPA_B_AMF",
                                                              "PNA_rpoB_AMF","PPA_rpoB_AMF","AMF richness"))

AMF.d <- as.vector(vegdist(env2$`AMF richness`,"bray"))
PPA_B_AMF.d <- as.vector(vegdist(env2$PPA_B_AMF,"bray"))
PNA_B_AMF.d <- as.vector(vegdist(env2$PNA_B_AMF,"bray"))
PPA_rpoB_AMF.d <- as.vector(vegdist(env2$PPA_rpoB_AMF,"bray"))
PNA_rpoB_AMF.d <- as.vector(vegdist(env2$PNA_rpoB_AMF,"bray"))

a <- read.csv("qpen.out.csv")
bNTI_B <- cbind.data.frame(a$result.bNTI,AMF.d,
                           PPA_B_AMF.d,PNA_B_AMF.d,
                           PPA_rpoB_AMF.d,PNA_rpoB_AMF.d)%>%`colnames<-`(value=c("bNTI","AMF.d",
                                                                                 "PPA_B_AMF.d",
                                                                                 "PNA_B_AMF.d",
                                                                                 "PPA_rpoB_AMF.d",
                                                                                 "PNA_rpoB_AMF.d"))

p1 <- ggplot(bNTI_B,aes(x=AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔAMF richness",y="Bacterial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = 0.34\nP < 0.001",x=0,y=3,hjust=0,vjust=0,size=4)

p2 <- ggplot(bNTI_B,aes(x=PPA_B_AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔPPA_B_AMF",y="Bacterial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = 0.25\nP = 0.002",x=0,y=3,hjust=0,vjust=0,size=4)

p3 <- ggplot(bNTI_B,aes(x=PPA_rpoB_AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔPPA_rpoB_AMF",y="Bacterial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = 0.33\nP < 0.001",x=0,y=3,hjust=0,vjust=0,size=4)

p4 <- ggplot(bNTI_B,aes(x=PNA_rpoB_AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔPNA_rpoB_AMF",y="Bacterial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = 0.33\nP < 0.001",x=0,y=3,hjust=0,vjust=0,size=4)

#####Rhizobia#####
b <- read.csv("qpen.out_rpoB.csv")
bNTI_rpoB <- cbind.data.frame(b$result.bNTI,AMF.d,
                              PPA_B_AMF.d,PNA_B_AMF.d,
                              PPA_rpoB_AMF.d,PNA_rpoB_AMF.d)%>%`colnames<-`(value=c("bNTI","AMF.d",
                                                                                    "PPA_B_AMF.d",
                                                                                    "PNA_B_AMF.d",
                                                                                    "PPA_rpoB_AMF.d",
                                                                                    "PNA_rpoB_AMF.d"))

p5 <- ggplot(bNTI_rpoB,aes(x=AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔAMF richness",y="Rhizobial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = -0.20\nP = 0.012",x=0,y=3,hjust=0,vjust=0,size=4)

p6 <- ggplot(bNTI_rpoB,aes(x=PNA_rpoB_AMF.d,y=bNTI))+
  geom_point(size=2.5,alpha=0.75)+
  geom_smooth(method="lm", colour="red",size=1.2,fill="lightgray")+
  geom_hline(yintercept = -2, linetype="dashed", color="grey")+
  geom_hline(yintercept = 2, linetype="dashed", color="grey")+
  labs(x="ΔPNA_rpoB_AMF",y="Rhizobial βNTI")+
  scale_y_continuous(limits = c(-4,4))+
  theme_cowplot()+
  theme(legend.position = "none")+
  annotate("text",label="Mantel's r = -0.29\nP < 0.001",x=0,y=3,hjust=0,vjust=0,size=4)

plot_grid( p1,p2,p3,p4,p5,p6,
           align = "hv",
           label_size = 15,
           labels = c("a","b","c","d","e","f"),
           nrow = 3,
           ncol=2
)

#####Assembly process ratio#####
CA_B <- a[1,1:5]
CA_rpoB <- b[1,1:5]

CA <- rbind.data.frame(CA_B,CA_rpoB)
rownames(CA) <- c("Bacteria","Rhizobia")
CA2 <- as.data.frame(t(CA))
CA2$group <- rownames(CA2)
CA2$group <- factor(CA2$group,levels = c( "ratio.Homogeneous.Selection","ratio.Heterogeneous.Selection",
                                          "ratio.Dispersal.Limitation","ratio.Homogenizing.Dispersal",
                                          "ratio.Undominated"),
                    labels = c("Homogeneous selection","Heterogeneous selection",
                               "Dispersal limitation","Homogenizing dispersal","Drift"))

x1 <- c(1.2,1.2,1.2,1.2,1.2)
x2 <- c(1.8,1.8,1.8,1.8,1.8)
y1 <- c(0,0.05228758,1-0.04575163-0.01307190,1-0.04575163,1)
y2 <- c(0,0.33986928,1-0.04575163,1,1)
segment <- cbind.data.frame(x1,x2,y1,y2)

CA2 %>% 
  pivot_longer(!group) %>% 
  mutate(name=factor(name,levels = c("Bacteria","Rhizobia"))) %>% 
  ggplot(aes(x=name,y=value))+
  geom_bar(aes(fill=group),
           stat="identity",
           position = "fill",
           width = 0.4)+
  geom_segment(data=segment,
               aes(x=x1,xend=x2,y=y1,yend=y2),
               lty="dashed",
               color="black")+
  scale_fill_npg()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100))+
  labs(y="Assembly processes (%)")+
  theme_cowplot()+
  theme(axis.title.x = element_blank(),legend.title = element_blank())