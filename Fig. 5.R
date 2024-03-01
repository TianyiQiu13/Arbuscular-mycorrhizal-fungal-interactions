library(microeco)
library(dplyr)
library(Hmisc)
library(vegan)
library(plspm)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(glmm.hp)

#####Bacteria#####
asv_B <- read.csv("ASV-B_qcmi.csv",row.names = 1)
group_B <- read.csv("group.csv",row.names = 1)
tax_B <- read.csv("taxonomy-B_qcmi.csv",row.names = 1)
data_B <- microtable$new(otu_table=asv_B,sample_table = group_B,tax_table = tax_B)

#####rpoB#####
asv_rpoB <- read.csv("ASV-rpoB_qcmi.csv",row.names = 1)
group_rpoB <- read.csv("group.csv",row.names = 1)
tax_rpoB <- read.csv("taxonomy-rpoB_qcmi.csv",row.names = 1)
data_rpoB <- microtable$new(otu_table=asv_rpoB,sample_table = group_rpoB,tax_table = tax_rpoB)

###Function prediction###
tmp <- clone(data_B)
tmp <- trans_func$new(tmp)
tmp$cal_spe_func(prok_database = "NJC19")
tmp$cal_spe_func_perc(abundance_weighted = T)
tmp$plot_spe_func_perc()
write.csv(tmp$res_spe_func_perc,"FAPROTAX_rpoB.csv")

tmp$show_prok_func(use_func = "ureolysis")

###Cor between func and PA###
func <- read.csv("FAPROTAX_B.csv",row.names = 1)%>%
  dplyr::select(c("aerobic_ammonia_oxidation","aerobic_nitrite_oxidation","chitinolysis",
                  "dark_hydrogen_oxidation","nitrogen_fixation","dark_oxidation_of_sulfur_compounds",
                  "aerobic_chemoheterotrophy","aromatic_compound_degradation","nitrate_respiration",
                  "nitrate_reduction","photoheterotrophy","ureolysis"))%>%
  `colnames<-`(value=c("Aerobic ammonia oxidation","Aerobic nitrite oxidation","Chitinolysis",
                       "Dark hydrogen oxidation","Nitrogen fixation",
                       "Dark oxidation of sulfur compounds","Aerobic chemoheterotrophy",
                       "Aromatic compound degradation","Nitrate respiration","Nitrate reduction",
                       "Photoheterotrophy","Ureolysis"))
PA <- read.csv("EMF-RF.csv")%>%dplyr::select(c("PNA_B_AMF","PPA_B_AMF","PNA_rpoB_AMF","PPA_rpoB_AMF"))
data <- rcorr(as.matrix(cbind.data.frame(func,PA)),type = "spearman")
r_value <- data$r
p_value <- data$P
r_related <- r_value[c(1:12),c(13:16)]
p_related <- p_value[c(1:12),c(13:16)]
r_related[p_related>0.05]=NA

func2 <- read.csv("Function_bar.csv")
func2$Position <- factor(func2$Position,levels = c("Top","Middle","Bottom"))
func2$type <- factor(func2$type,
                     levels=rev(c("aerobic_ammonia_oxidation","aerobic_nitrite_oxidation","chitinolysis",
                              "dark_hydrogen_oxidation","nitrogen_fixation",
                              "dark_oxidation_of_sulfur_compounds","aerobic_chemoheterotrophy",
                              "aromatic_compound_degradation","nitrate_respiration",
                              "nitrate_reduction","photoheterotrophy","ureolysis")),
                     labels = rev(c("Aerobic ammonia oxidation","Aerobic nitrite oxidation",
                                    "Chitinolysis","Dark hydrogen oxidation","Nitrogen fixation",
                                "Dark oxidation of sulfur compounds","Aerobic chemoheterotrophy",
                                "Aromatic compound degradation","Nitrate respiration",
                                "Nitrate reduction","Photoheterotrophy","Ureolysis")))

anova(lm(value~Position,data=filter(func2,type=="Ureolysis")))

p1 <- ggplot(func2, aes(x=type, y=value)) + 
  geom_bar(position = "dodge2", stat = "summary",aes(fill=Position),width = 0.75) + 
  geom_errorbar(position="dodge2",stat="summary",aes(color=Position),width=0.75)+
  scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  ylab("Relative abundance (%)")+
  theme_bw() +
  coord_flip()+
  theme(axis.line = element_line(), axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),legend.background = element_blank(),
        legend.key.size = unit(1,"lines"),
        legend.position = c(-0.875,0.15),legend.title = element_blank())
p1

p2 <- ggballoonplot(r_related,fill="value",size.range = 4,color="white")+
  gradient_fill(colorRampPalette(brewer.pal(n=7,name="RdYlBu"))(100))+
  guides(size=F)+
  theme(axis.text.x = element_text(angle = 20),axis.title.x = element_blank(),
        legend.title = element_blank(),legend.key.size = unit(1,"lines"),
        legend.position = c(-0.925,0.15))
p2

p3 <- plot_grid( p1,p2,
           align = "v",
           label_size = 15,
           labels = c("a","b"),
           nrow = 2,
           ncol=1,
           rel_heights = c(0.95,1)
)
p3

###PLS-PM###
env <- read.csv("environment.csv",row.names = 1)
spe <- read.csv("EMF-RF.csv")%>%dplyr::select(c("Richness_AMF","CAP1_B","CAP1_rpoB"))
PA <- read.csv("EMF-RF.csv")%>%dplyr::select(c("PNA_B_AMF","PPA_B_AMF","PNA_rpoB_AMF","PPA_rpoB_AMF"))
EMF <- read.csv("EMF-RF.csv")%>%dplyr::select(c("Position","EMF"))%>%
  mutate(Position=c(1,1,1,1,1,1,
                    2,2,2,2,2,2,
                    3,3,3,3,3,3),plot=c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6))

aa <- cbind.data.frame(env,spe,PA,EMF)
path <- read.csv("path_EMF.csv",row.names = 1)%>%as.matrix()
blocks = list(14,1:6,7,8:9,10:13,15)
modes = c("A","A","A","A","A","A")
pls <- plspm(aa, path, blocks, modes = modes)
summary(pls)
write.csv(pls$effects,"pls_EMF.csv")

ef <- read.csv("pls_EMF.csv",row.names = 1)
ef <- ef[c(5,9,12,14,15),]
ef$relationships <- factor(ef$relationships,
                           levels = c("Position -> EMF","Env -> EMF","AMF -> EMF",
                                      "Comp -> EMF","PA -> EMF"),
                           labels = c("Position","Environment","AMF richness",
                                      "Microbial composition",
                                      "Putative biotic\nassociations"))
ef_d <- ef[,1:2]%>%`colnames<-`(value=c("var","effect"))%>%mutate(type="Direct effect")
ef_i <- ef[,c(1,3)]%>%`colnames<-`(value=c("var","effect"))%>%mutate(type="Indirect effect")
ef2 <- rbind.data.frame(ef_d,ef_i)
ef2 <- ef2[c(1,2,3,5,6,7,9),]
p4 <- ggplot(ef2, aes(x=var, y=effect)) + 
    geom_bar(position = "dodge2", stat = "summary",aes(fill=var),width = 0.75,alpha=0.85) + 
  scale_fill_manual(values = c("#515366","#226258","#67448D","#957A2A","#8A2A44"))+
  ylab("Standardized effects from PLS-PM")+
  theme_bw()+
  facet_wrap(~type,scales = "free_x")+
  theme(axis.line = element_line(), axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),legend.background = element_blank(),
        legend.key.size = unit(1,"lines"),axis.text.x = element_text(angle = 20,hjust = 1),
        legend.position = "none")
p4

ggdraw()+ 
  draw_plot(p3, x=0, y=0, width = 0.5, height = 1)+
  draw_plot(p4, x=0.5, y=0, width = 0.5, height = 0.5)+
  draw_plot_label(label=c("c","d"),size = 15,
                  x=c(0.5,0.5),y=c(1,0.5))

###glmm.hp###
mod <- lmerTest::lmer(EMF~Moisture+pH+BD+SOC+N.P+TP+Richness_AMF+CAP1_B+CAP1_rpoB+PNA_B_AMF+PPA_B_AMF+
               PNA_rpoB_AMF+PPA_rpoB_AMF+Position+(1|plot),data=aa)
r.squaredGLMM(mod)##0.9184983
result <- glmm.hp(mod)
result$hierarchical.partitioning