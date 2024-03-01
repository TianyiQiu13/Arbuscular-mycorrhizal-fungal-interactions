library(ggplot2)
library(cowplot)
library(lme4)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)

a <- read.csv("Slope BEF data.csv")%>%dplyr::select(c("DOC","MBC","TDN","Olsen.P","MBN","MBP",
                                                      "BG","NAG","LAP","AP","LeafC","LitterC"))
b <- read.csv("EMF-RF.csv")%>%dplyr::select(c("Richness_AMF","CAP1_B","CAP1_rpoB","EMF"))

data <- rcorr(as.matrix(cbind.data.frame(a,b)),type="spearman")
r_value <- data$r
p_value <- data$P
r_related <- r_value[c(13:15),c(1:12,16)]
p_related <- p_value[c(13:15),c(1:12,16)]
r_related[p_related>0.05]=NA

cols <- c(colorRampPalette(brewer.pal(9, 'RdYlBu'))(12))

p1 <- pheatmap(r_related,fontsize_number = 12,fontsize = 12,cluster_cols = F,cluster_rows = F,
               border= "white",display_numbers = T,number_color="white",color = cols,
               labels_row = c("AMF\nrichness", "Bacteria-\nCPCoA1","Rhizobia-\nCPCoA1"),
               labels_col = c("DOC","MBC","TDN","Olsen P","MBN","MBP",
                              "BG","NAG","LAP","AP","Leaf C","Litter C","Multi-\nfunctionality"),
               angle_col = 45,na_col = "white")

g1 <- as.ggplot(p1)

p2 <- plot_grid( g1 + theme(plot.margin = unit(c(0.05,0,0,0.05), "cm"))
)

ggdraw()+ 
  draw_plot(p2, x=0, y=0, width = 1, height = 1)+
  draw_label("Microbially\ndriven C pools",x=0.025,y=0.055,size=13,fontface = "bold",hjust = 0)+
  draw_label("Nutrient\ncycling",x=0.25,y=0.055,size=13,fontface = "bold",hjust = 0)+
  draw_label("Organic matter\ndecomposition",x=0.475,y=0.055,size=13,fontface = "bold",hjust = 0)+
  draw_label("Plant\nproduction",x=0.685,y=0.055,size=13,fontface = "bold",hjust = 0)