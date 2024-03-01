library(qcmi)
library(igraph)
library(ggraph)
library(ggnetwork)
library(sna)
library(geomnet)

env <- read.csv("environment.csv",row.names = 1)
geo <- read.csv("geo.csv",row.names = 1)

#####B+AMF#####
##Network calculation##
asv_B <- read.csv("ASV-B_filtering.csv",row.names = 1)
group_B <- read.csv("group.csv",row.names = 1)
tax_B <- read.csv("taxonomy-B_filtering.csv",row.names = 1)
ps_B <- trans_ps(otu_table=asv_B,taxa_table=tax_B,sample_data = sample_data(group_B))
ps_B <- filter_ps(ps_B,abundance.threshold = 500)
ig_B <- read_graph("B_microeco.graphml",format = "graphml")

asv_AMF <- read.csv("ASV-AMF_filtering.csv",row.names = 1)
group_AMF <- read.csv("group.csv",row.names = 1)
tax_AMF <- read.csv("taxonomy-AMF_filtering.csv",row.names = 1)
ps_AMF <- trans_ps(otu_table=asv_AMF,taxa_table=tax_AMF,sample_data = sample_data(group_AMF))
ps_AMF <- filter_ps(ps_AMF,abundance.threshold = 500)
ig_AMF <- read_graph("AMF_microeco.graphml",format = "graphml")

ps_B_AMF <- merge16S_ITS(ps16s = ps_B,psITS = ps_AMF,N16s = 0,NITS = 0,dat1.lab = "B",dat2.lab = "AMF")
edgelist_B_AMF <- read.csv("B_AMF_microeco_edge.csv",row.names = 1)
ig_B_AMF <- graph_from_edgelist(as.matrix(edgelist_B_AMF[,1:2]),directed = F)
ig_B_AMF <- set_edge_attr(ig_B_AMF,"weight",index=E(ig_B_AMF),as.numeric(edgelist_B_AMF[,3]))

##Within- and cross-kingdom subsets##
edgei_B <- edgelist_B_AMF%>%filter(i=="B")%>%filter(j=="B")
igi_B <- graph_from_edgelist(as.matrix(edgei_B[,1:2]),directed = F)
igi_B <- set_edge_attr(igi_B,"weight",index=E(igi_B),as.numeric(edgei_B[,3]))
write.csv(net_properties(igi_B),"B_i_qcmi_netproperties.csv")

edgei_B_AMF <- edgelist_B_AMF%>%filter(i=="B")%>%filter(j=="A")
igi_B_AMF <- graph_from_edgelist(as.matrix(edgei_B_AMF[,1:2]),directed = F)
igi_B_AMF <- set_edge_attr(igi_B_AMF,"weight",index=E(igi_B_AMF),as.numeric(edgei_B_AMF[,3]))
write.csv(net_properties(igi_B_AMF),"B_AMF_i_qcmi_netproperties.csv")

##Sample subsets##
sub_graph_B <- list()
for (i in colnames(t(vegan_otu(ps_B)))) {
  otu <- t(vegan_otu(ps_B))
  sample_i <- otu[,i]
  select_node <- rownames(data.frame(sample_i))[which(data.frame(sample_i) != 0)]
  sub_graph_B[[i]] <- induced_subgraph(ig_B, select_node)
}
net_sub_B <- data.frame(net_properties(sub_graph_B[[1]]))
for(i in 2:length(sub_graph_B)){
  net_sub_B <- cbind.data.frame(net_sub_B,net_properties(sub_graph_B[[i]]))
}
colnames(net_sub_B) <- names(sub_graph_B)
density_sub_B <- data.frame(edge_density(sub_graph_B[[1]]))
for(i in 2:length(sub_graph_B)){
  density_sub_B <- cbind.data.frame(density_sub_B,edge_density(sub_graph_B[[i]]))
}
colnames(density_sub_B) <- names(sub_graph_B)
rownames(density_sub_B) <- "density"
net_sub_B <- rbind.data.frame(net_sub_B,density_sub_B)
write.csv(t(net_sub_B),"net_sub_B.csv")

centrality_sub_B <- data.frame(mean(hub_score(sub_graph_B[[1]])$vector))
for(i in 2:length(sub_graph_B)){
  centrality_sub_B <- cbind.data.frame(centrality_sub_B,mean(hub_score(sub_graph_B[[i]])$vector))
}
colnames(centrality_sub_B) <- names(sub_graph_B)
write.csv(t(data.frame(centrality_sub_B)),"centrality_sub_B.csv")

sub_graph_B_AMF <- list()
for (i in colnames(t(vegan_otu(ps_B_AMF)))) {
  otu <- t(vegan_otu(ps_B_AMF))
  sample_i <- otu[,i]
  select_node <- rownames(data.frame(sample_i))[which(data.frame(sample_i) != 0)]
  sub_graph_B_AMF[[i]] <- induced_subgraph(ig_B_AMF, select_node)
}
net_sub_B_AMF <- data.frame(net_properties(sub_graph_B_AMF[[1]]))
for(i in 2:length(sub_graph_B_AMF)){
  net_sub_B_AMF <- cbind.data.frame(net_sub_B_AMF,net_properties(sub_graph_B_AMF[[i]]))
}
colnames(net_sub_B_AMF) <- names(sub_graph_B_AMF)
density_sub_B_AMF <- data.frame(edge_density(sub_graph_B_AMF[[1]]))
for(i in 2:length(sub_graph_B_AMF)){
  density_sub_B_AMF <- cbind.data.frame(density_sub_B_AMF,edge_density(sub_graph_B_AMF[[i]]))
}
colnames(density_sub_B_AMF) <- names(sub_graph_B_AMF)
rownames(density_sub_B_AMF) <- "density"
net_sub_B_AMF <- rbind.data.frame(net_sub_B_AMF,density_sub_B_AMF)
write.csv(t(net_sub_B_AMF),"net_sub_B_AMF.csv")

centrality_sub_B_AMF <- data.frame(mean(hub_score(sub_graph_B_AMF[[1]])$vector))
for(i in 2:length(sub_graph_B_AMF)){
  centrality_sub_B_AMF <- cbind.data.frame(centrality_sub_B_AMF,mean(hub_score(sub_graph_B_AMF[[i]])$vector))
}
colnames(centrality_sub_B_AMF) <- names(sub_graph_B_AMF)
write.csv(t(data.frame(centrality_sub_B_AMF)),"centrality_sub_B_AMF.csv")

####qcmi####
otu_table_B_AMF <- as.data.frame(t(vegan_otu(ps_B_AMF)))

results_dl_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                     p=0.05,data=geo,cutoff = 0,method = "dl")
results_moisture_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                           p=0.05,data=env['Moisture'],cutoff = 0,method = "ef")
results_pH_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                     p=0.05,data=env['pH'],cutoff = 0,method = "ef")
results_BD_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                     p=0.05,data=env['BD'],cutoff = 0,method = "ef")
results_SOC_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                      p=0.05,data=env['SOC'],cutoff = 0,method = "ef")
results_N.P_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                      p=0.05,data=env['N.P'],cutoff = 0,method = "ef")
results_TP_B_AMF <- assigned_process(link_table_row = edgelist_B_AMF,OTUabd = otu_table_B_AMF,
                                     p=0.05,data=env['TP'],cutoff = 0,method = "ef")

total_link_B_AMF <- row.names(edgelist_B_AMF)
ef_link_B_AMF <- union(as.character(row.names(results_moisture_B_AMF)),
                       as.character(row.names(results_pH_B_AMF)))
ef_link_B_AMF <- union(ef_link_B_AMF,as.character(row.names(results_BD_B_AMF)))
ef_link_B_AMF <- union(ef_link_B_AMF,as.character(row.names(results_SOC_B_AMF)))
ef_link_B_AMF <- union(ef_link_B_AMF,as.character(row.names(results_N.P_B_AMF)))
ef_link_B_AMF <- union(ef_link_B_AMF,as.character(row.names(results_TP_B_AMF)))
dl_link_B_AMF <- as.character(row.names(results_dl_B_AMF))
bi_link_B_AMF <- setdiff(total_link_B_AMF,ef_link_B_AMF)
bi_link_B_AMF <- setdiff(bi_link_B_AMF,dl_link_B_AMF)
efedge_B_AMF <- edgelist_B_AMF[ef_link_B_AMF,]
dledge_B_AMF <- edgelist_B_AMF[dl_link_B_AMF,]
biedge_B_AMF <- edgelist_B_AMF[bi_link_B_AMF,]

ig.bi_B_AMF <- graph_from_edgelist(as.matrix(biedge_B_AMF[,1:2]),directed = F)
ig.bi_B_AMF <- set_edge_attr(ig.bi_B_AMF,"weight",index=E(ig.bi_B_AMF),as.numeric(biedge_B_AMF[,3]))
results_bi_B_AMF <- qcmi(igraph = ig.bi_B_AMF,OTU = otu_table_B_AMF,pers.cutoff = 0)
bi_B_AMF <- cbind.data.frame(results_bi_B_AMF[[3]],results_bi_B_AMF[[4]])
colnames(bi_B_AMF) <- c("PNA_B_AMF","PPA_B_AMF")
write.csv(bi_B_AMF,"B_AMF_PA2.csv")

#####rpoB+AMF#####
##Network calculation##
asv_rpoB <- read.csv("ASV-rpoB_filtering.csv",row.names = 1)
group_rpoB <- read.csv("group.csv",row.names = 1)
tax_rpoB <- read.csv("taxonomy-rpoB_filtering.csv",row.names = 1)
ps_rpoB <- trans_ps(otu_table=asv_rpoB,taxa_table=tax_rpoB,sample_data = sample_data(group_rpoB))
ps_rpoB <- filter_ps(ps_rpoB,abundance.threshold = 500)
ig_rpoB <- read_graph("rpoB_microeco.graphml",format = "graphml")

edgelist_rpoB_AMF <- read.csv("rpoB_AMF_microeco_edge.csv",row.names = 1)
ig_rpoB_AMF <- graph_from_edgelist(as.matrix(edgelist_rpoB_AMF[,1:2]),directed = F)
ig_rpoB_AMF <- set_edge_attr(ig_rpoB_AMF,"weight",index=E(ig_rpoB_AMF),as.numeric(edgelist_rpoB_AMF[,3]))

##Within- and cross-kingdom subsets##
edgei_rpoB <- edgelist_rpoB_AMF%>%filter(i=="r")%>%filter(j=="r")
igi_rpoB <- graph_from_edgelist(as.matrix(edgei_rpoB[,1:2]),directed = F)
igi_rpoB <- set_edge_attr(igi_rpoB,"weight",index=E(igi_rpoB),as.numeric(edgei_rpoB[,3]))
write.csv(net_properties(igi_rpoB),"rpoB_i_qcmi_netproperties.csv")

ps_rpoB_AMF <- merge16S_ITS(ps16s = ps_rpoB,psITS = ps_AMF,N16s = 0,NITS = 0,dat1.lab = "rpoB",dat2.lab = "AMF")
edgei_rpoB_AMF <- edgelist_rpoB_AMF%>%filter(i=="r")%>%filter(j=="A")
igi_rpoB_AMF <- graph_from_edgelist(as.matrix(edgei_rpoB_AMF[,1:2]),directed = F)
igi_rpoB_AMF <- set_edge_attr(igi_rpoB_AMF,"weight",index=E(igi_rpoB_AMF),as.numeric(edgei_rpoB_AMF[,3]))
write.csv(net_properties(igi_rpoB_AMF),"rpoB_AMF_i_qcmi_netproperties.csv")

##Sample subsets##
sub_graph_rpoB <- list()
for (i in colnames(t(vegan_otu(ps_rpoB)))) {
  otu <- t(vegan_otu(ps_rpoB))
  sample_i <- otu[,i]
  select_node <- rownames(data.frame(sample_i))[which(data.frame(sample_i) != 0)]
  sub_graph_rpoB[[i]] <- induced_subgraph(ig_rpoB, select_node)
}
net_sub_rpoB <- data.frame(net_properties(sub_graph_rpoB[[1]]))
for(i in 2:length(sub_graph_rpoB)){
  net_sub_rpoB <- cbind.data.frame(net_sub_rpoB,net_properties(sub_graph_rpoB[[i]]))
}
colnames(net_sub_rpoB) <- names(sub_graph_rpoB)
density_sub_rpoB <- data.frame(edge_density(sub_graph_rpoB[[1]]))
for(i in 2:length(sub_graph_rpoB)){
  density_sub_rpoB <- cbind.data.frame(density_sub_rpoB,edge_density(sub_graph_rpoB[[i]]))
}
colnames(density_sub_rpoB) <- names(sub_graph_rpoB)
rownames(density_sub_rpoB) <- "density"
net_sub_rpoB <- rbind.data.frame(net_sub_rpoB,density_sub_rpoB)
write.csv(t(net_sub_rpoB),"net_sub_rpoB.csv")

centrality_sub_rpoB <- data.frame(mean(hub_score(sub_graph_rpoB[[1]])$vector))
for(i in 2:length(sub_graph_rpoB)){
  centrality_sub_rpoB <- cbind.data.frame(centrality_sub_rpoB,mean(hub_score(sub_graph_rpoB[[i]])$vector))
}
colnames(centrality_sub_rpoB) <- names(sub_graph_rpoB)
write.csv(t(data.frame(centrality_sub_rpoB)),"centrality_sub_rpoB.csv")

sub_graph_rpoB_AMF <- list()
for (i in colnames(t(vegan_otu(ps_rpoB_AMF)))) {
  otu <- t(vegan_otu(ps_rpoB_AMF))
  sample_i <- otu[,i]
  select_node <- rownames(data.frame(sample_i))[which(data.frame(sample_i) != 0)]
  sub_graph_rpoB_AMF[[i]] <- induced_subgraph(ig_rpoB_AMF, select_node)
}
net_sub_rpoB_AMF <- data.frame(net_properties(sub_graph_rpoB_AMF[[1]]))
for(i in 2:length(sub_graph_rpoB_AMF)){
  net_sub_rpoB_AMF <- cbind.data.frame(net_sub_rpoB_AMF,net_properties(sub_graph_rpoB_AMF[[i]]))
}
colnames(net_sub_rpoB_AMF) <- names(sub_graph_rpoB_AMF)
density_sub_rpoB_AMF <- data.frame(edge_density(sub_graph_rpoB_AMF[[1]]))
for(i in 2:length(sub_graph_rpoB_AMF)){
  density_sub_rpoB_AMF <- cbind.data.frame(density_sub_rpoB_AMF,edge_density(sub_graph_rpoB_AMF[[i]]))
}
colnames(density_sub_rpoB_AMF) <- names(sub_graph_rpoB_AMF)
rownames(density_sub_rpoB_AMF) <- "density"
net_sub_rpoB_AMF <- rbind.data.frame(net_sub_rpoB_AMF,density_sub_rpoB_AMF)
write.csv(t(net_sub_rpoB_AMF),"net_sub_rpoB_AMF.csv")

centrality_sub_rpoB_AMF <- data.frame(mean(hub_score(sub_graph_rpoB_AMF[[1]])$vector))
for(i in 2:length(sub_graph_rpoB_AMF)){
  centrality_sub_rpoB_AMF <- cbind.data.frame(centrality_sub_rpoB_AMF,mean(hub_score(sub_graph_rpoB_AMF[[i]])$vector))
}
colnames(centrality_sub_rpoB_AMF) <- names(sub_graph_rpoB_AMF)
write.csv(t(data.frame(centrality_sub_rpoB_AMF)),"centrality_sub_rpoB_AMF.csv")

####qcmi####
otu_table_rpoB_AMF <- as.data.frame(t(vegan_otu(ps_rpoB_AMF)))

results_dl_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                        p=0.05,data=geo,cutoff = 0,method = "dl")
results_moisture_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                              p=0.05,data=env['Moisture'],cutoff = 0,method = "ef")
results_pH_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                        p=0.05,data=env['pH'],cutoff = 0,method = "ef")
results_rpoBD_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                           p=0.05,data=env['BD'],cutoff = 0,method = "ef")
results_SOC_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                         p=0.05,data=env['SOC'],cutoff = 0,method = "ef")
results_N.P_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                         p=0.05,data=env['N.P'],cutoff = 0,method = "ef")
results_TP_rpoB_AMF <- assigned_process(link_table_row = edgelist_rpoB_AMF,OTUabd = otu_table_rpoB_AMF,
                                        p=0.05,data=env['TP'],cutoff = 0,method = "ef")

total_link_rpoB_AMF <- row.names(edgelist_rpoB_AMF)
ef_link_rpoB_AMF <- union(as.character(row.names(results_moisture_rpoB_AMF)),
                          as.character(row.names(results_pH_rpoB_AMF)))
ef_link_rpoB_AMF <- union(ef_link_rpoB_AMF,as.character(row.names(results_rpoBD_rpoB_AMF)))
ef_link_rpoB_AMF <- union(ef_link_rpoB_AMF,as.character(row.names(results_SOC_rpoB_AMF)))
ef_link_rpoB_AMF <- union(ef_link_rpoB_AMF,as.character(row.names(results_N.P_rpoB_AMF)))
ef_link_rpoB_AMF <- union(ef_link_rpoB_AMF,as.character(row.names(results_TP_rpoB_AMF)))
dl_link_rpoB_AMF <- as.character(row.names(results_dl_rpoB_AMF))
bi_link_rpoB_AMF <- setdiff(total_link_rpoB_AMF,ef_link_rpoB_AMF)
bi_link_rpoB_AMF <- setdiff(bi_link_rpoB_AMF,dl_link_rpoB_AMF)
efedge_rpoB_AMF <- edgelist_rpoB_AMF[ef_link_rpoB_AMF,]
dledge_rpoB_AMF <- edgelist_rpoB_AMF[dl_link_rpoB_AMF,]
biedge_rpoB_AMF <- edgelist_rpoB_AMF[bi_link_rpoB_AMF,]

ig.bi_rpoB_AMF <- graph_from_edgelist(as.matrix(biedge_rpoB_AMF[,1:2]),directed = F)
ig.bi_rpoB_AMF <- set_edge_attr(ig.bi_rpoB_AMF,"weight",index=E(ig.bi_rpoB_AMF),as.numeric(biedge_rpoB_AMF[,3]))
results_rpoBi_rpoB_AMF <- qcmi(igraph = ig.bi_rpoB_AMF,OTU = otu_table_rpoB_AMF,pers.cutoff = 0)
bi_rpoB_AMF <- cbind.data.frame(results_rpoBi_rpoB_AMF[[3]],results_rpoBi_rpoB_AMF[[4]])
colnames(bi_rpoB_AMF) <- c("PNA_rpoB_AMF","PPA_rpoB_AMF")
write.csv(bi_rpoB_AMF,"rpoB_AMF_PA2.csv")