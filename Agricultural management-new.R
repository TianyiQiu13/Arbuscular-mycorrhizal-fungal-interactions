library(ggrepel)
library(rcartocolor)
library(metafor)
library(boot)
library(agricolae)
library(readxl)
library(broom)
library(ggforestplot)
library(funModeling)
library(dplyr)
library(cowplot)

make_pct <- function(x) (exp(x) - 1) * 100
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12))#11.6inches
a <- read.csv("S3-preparation.csv")
a$ID <- as.factor(a$ID)
a$study <- as.factor(a$study)

###Calculation###
####OF####
yield.func1 <- function(a, index) {
  yield <- rma.mv(rryield,varyield,data=a,subset=index,mods=~OF2-1,random=~1|study/ID)
  c(yield$beta, yield$se^2)}
SOC.func1 <- function(a, index) {
  SOC <- rma.mv(rrSOC,varSOC,data=a,subset=index,mods=~OF2-1,random=~1|study/ID)
  c(SOC$beta, SOC$se^2)}
N2O.func1 <- function(a, index) {
  N2O <- rma.mv(rrN2O,varN2O,data=a,subset=index,mods=~OF2-1,random=~1|study/ID)
  c(N2O$beta, N2O$se^2)}
CH4.func1 <- function(a, index) {
  CH4 <- rma.mv(rrCH4,varCH4,data=a,subset=index,mods=~OF2-1,random=~1|study/ID)
  c(CH4$beta, CH4$se^2)}
MWD.func1 <- function(a, index) {
  MWD <- rma.mv(rrMWD,varMWD,data=a,subset=index,mods=~OF2-1,random=~1|study/ID)
  c(MWD$beta, MWD$se^2)}

yield_boot1 <- boot(a, yield.func1, R=49)
SOC_boot1 <- boot(a, SOC.func1, R=49)
N2O_boot1 <- boot(a, N2O.func1, R=49)
CH4_boot1 <- boot(a, CH4.func1, R=49)
MWD_boot1 <- boot(a, MWD.func1, R=49)

yield_ci <- tidy(yield_boot1,conf.int=T,conf.method = "perc")
yield_ci <- data.frame(yield_ci)
yield_ci$mean <- (yield_ci$conf.low+yield_ci$conf.high)/2
yield_ci <- yield_ci[1:2,3:6]
yield_output <- make_pct(yield_ci)
SOC_ci <- tidy(SOC_boot1,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1:2,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot1,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1:2,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot1,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1:2,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot1,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1:2,3:6]
MWD_output <- make_pct(MWD_ci)

yield_p <- rma.mv(rryield,varyield,data=a,mods=~OF2-1,random=~1|study/ID)
SOC_p <- rma.mv(rrSOC,varSOC,data=a,mods=~OF2-1,random=~1|study/ID)
N2O_p <- rma.mv(rrN2O,varN2O,data=a,mods=~OF2-1,random=~1|study/ID)
CH4_p <- rma.mv(rrCH4,varCH4,data=a,mods=~OF2-1,random=~1|study/ID)
MWD_p <- rma.mv(rrMWD,varMWD,data=a,mods=~OF2-1,random=~1|study/ID)

b <- rbind.data.frame(yield_p$pval[1],yield_p$pval[2],SOC_p$pval[1],SOC_p$pval[2],
                      N2O_p$pval[1],N2O_p$pval[2],CH4_p$pval[1],CH4_p$pval[2],
                      MWD_p$pval[1],MWD_p$pval[2])
names(b) <- "pval"
c <- rbind.data.frame(yield_output,SOC_output,N2O_output,CH4_output,MWD_output)
d <- cbind.data.frame(b,c)
d$variable <- c("Yield","Yield","SOC","SOC","N2O","N2O","CH4","CH4","MWD","MWD")
d$type <- c("No","Yes")
#write.csv(d,"S3-OF.csv")

####tillage####
yield.func2 <- function(a, index) {
  yield <- rma.mv(rryield,varyield,data=a,subset=index,mods=~tillage2-1,random=~1|study/ID)
  c(yield$beta, yield$se^2)}
SOC.func2 <- function(a, index) {
  SOC <- rma.mv(rrSOC,varSOC,data=a,subset=index,mods=~tillage2-1,random=~1|study/ID)
  c(SOC$beta, SOC$se^2)}
N2O.func2 <- function(a, index) {
  N2O <- rma.mv(rrN2O,varN2O,data=a,subset=index,mods=~tillage2-1,random=~1|study/ID)
  c(N2O$beta, N2O$se^2)}
CH4.func2 <- function(a, index) {
  CH4 <- rma.mv(rrCH4,varCH4,data=a,subset=index,mods=~tillage2-1,random=~1|study/ID)
  c(CH4$beta, CH4$se^2)}
MWD.func2 <- function(a, index) {
  MWD <- rma.mv(rrMWD,varMWD,data=a,subset=index,mods=~tillage2-1,random=~1|study/ID)
  c(MWD$beta, MWD$se^2)}

yield_boot2 <- boot(a, yield.func2, R=49)
SOC_boot2 <- boot(a, SOC.func2, R=49)
N2O_boot2 <- boot(a, N2O.func2, R=49)
CH4_boot2 <- boot(a, CH4.func2, R=49)
MWD_boot2 <- boot(a, MWD.func2, R=49)

yield_ci <- tidy(yield_boot2,conf.int=T,conf.method = "perc")
yield_ci <- data.frame(yield_ci)
yield_ci$mean <- (yield_ci$conf.low+yield_ci$conf.high)/2
yield_ci <- yield_ci[1:2,3:6]
yield_output <- make_pct(yield_ci)
SOC_ci <- tidy(SOC_boot2,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1:2,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot2,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1:2,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot2,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1:2,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot2,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1:2,3:6]
MWD_output <- make_pct(MWD_ci)

yield_p <- rma.mv(rryield,varyield,data=a,mods=~tillage2-1,random=~1|study/ID)
SOC_p <- rma.mv(rrSOC,varSOC,data=a,mods=~tillage2-1,random=~1|study/ID)
N2O_p <- rma.mv(rrN2O,varN2O,data=a,mods=~tillage2-1,random=~1|study/ID)
CH4_p <- rma.mv(rrCH4,varCH4,data=a,mods=~tillage2-1,random=~1|study/ID)
MWD_p <- rma.mv(rrMWD,varMWD,data=a,mods=~tillage2-1,random=~1|study/ID)

e <- rbind.data.frame(yield_p$pval[1],yield_p$pval[2],SOC_p$pval[1],SOC_p$pval[2],
                      N2O_p$pval[1],N2O_p$pval[2],CH4_p$pval[1],CH4_p$pval[2],
                      MWD_p$pval[1],MWD_p$pval[2])
names(e) <- "pval"
f <- rbind.data.frame(yield_output,SOC_output,N2O_output,CH4_output,MWD_output)
g <- cbind.data.frame(e,f)
g$variable <- c("Yield","Yield","SOC","SOC","N2O","N2O","CH4","CH4","MWD","MWD")
g$type <- c("Conservation","Conventional")
#write.csv(g,"S3-tillage.csv")

####irrigation####
#upland
upland <- a%>%filter(cropland.type=="upland")
yield.func3 <- function(upland, index) {
  yield <- rma.mv(rryield,varyield,data=upland,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(yield$beta, yield$se^2)}
SOC.func3 <- function(upland, index) {
  SOC <- rma.mv(rrSOC,varSOC,data=upland,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(SOC$beta, SOC$se^2)}
N2O.func3 <- function(upland, index) {
  N2O <- rma.mv(rrN2O,varN2O,data=upland,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(N2O$beta, N2O$se^2)}
CH4.func3 <- function(upland, index) {
  CH4 <- rma.mv(rrCH4,varCH4,data=upland,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(CH4$beta, CH4$se^2)}
MWD.func3 <- function(upland, index) {
  MWD <- rma.mv(rrMWD,varMWD,data=upland,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(MWD$beta, MWD$se^2)}

yield_boot3 <- boot(upland, yield.func3, R=49)
SOC_boot3 <- boot(upland, SOC.func3, R=49)
N2O_boot3 <- boot(upland, N2O.func3, R=49)
CH4_boot3 <- boot(upland, CH4.func3, R=49)
MWD_boot3 <- boot(upland, MWD.func3, R=49)

yield_ci <- tidy(yield_boot3,conf.int=T,conf.method = "perc")
yield_ci <- data.frame(yield_ci)
yield_ci$mean <- (yield_ci$conf.low+yield_ci$conf.high)/2
yield_ci <- yield_ci[1:2,3:6]
yield_output <- make_pct(yield_ci)
SOC_ci <- tidy(SOC_boot3,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1:2,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot3,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1:2,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot3,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1:2,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot3,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1:2,3:6]
MWD_output <- make_pct(MWD_ci)

yield_p <- rma.mv(rryield,varyield,data=upland,mods=~irrigation2-1,random=~1|study/ID)
SOC_p <- rma.mv(rrSOC,varSOC,data=upland,mods=~irrigation2-1,random=~1|study/ID)
N2O_p <- rma.mv(rrN2O,varN2O,data=upland,mods=~irrigation2-1,random=~1|study/ID)
CH4_p <- rma.mv(rrCH4,varCH4,data=upland,mods=~irrigation2-1,random=~1|study/ID)
MWD_p <- rma.mv(rrMWD,varMWD,data=upland,mods=~irrigation2-1,random=~1|study/ID)

h <- rbind.data.frame(yield_p$pval[1],yield_p$pval[2],SOC_p$pval[1],SOC_p$pval[2],
                      N2O_p$pval[1],N2O_p$pval[2],CH4_p$pval[1],CH4_p$pval[2],
                      MWD_p$pval[1],MWD_p$pval[2])
names(h) <- "pval"
i <- rbind.data.frame(yield_output,SOC_output,N2O_output,CH4_output,MWD_output)
j <- cbind.data.frame(h,i)
j$variable <- c("Yield","Yield","SOC","SOC","N2O","N2O","CH4","CH4","MWD","MWD")
j$type <- c("Irrigated","Rainfed")
#write.csv(j,"S3-irrigation-upland.csv")
#paddy
paddy <- a%>%filter(cropland.type=="paddy")
yield.func4 <- function(paddy, index) {
  yield <- rma.mv(rryield,varyield,data=paddy,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(yield$beta, yield$se^2)}
SOC.func4 <- function(paddy, index) {
  SOC <- rma.mv(rrSOC,varSOC,data=paddy,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(SOC$beta, SOC$se^2)}
N2O.func4 <- function(paddy, index) {
  N2O <- rma.mv(rrN2O,varN2O,data=paddy,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(N2O$beta, N2O$se^2)}
CH4.func4 <- function(paddy, index) {
  CH4 <- rma.mv(rrCH4,varCH4,data=paddy,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(CH4$beta, CH4$se^2)}
MWD.func4 <- function(paddy, index) {
  MWD <- rma.mv(rrMWD,varMWD,data=paddy,subset=index,mods=~irrigation2-1,random=~1|study/ID)
  c(MWD$beta, MWD$se^2)}

yield_boot4 <- boot(paddy, yield.func4, R=49)
SOC_boot4 <- boot(paddy, SOC.func4, R=49)
N2O_boot4 <- boot(paddy, N2O.func4, R=49)
CH4_boot4 <- boot(paddy, CH4.func4, R=49)
MWD_boot4 <- boot(paddy, MWD.func4, R=49)

yield_ci <- tidy(yield_boot4,conf.int=T,conf.method = "perc")
yield_ci <- data.frame(yield_ci)
yield_ci$mean <- (yield_ci$conf.low+yield_ci$conf.high)/2
yield_ci <- yield_ci[1:3,3:6]
yield_output <- make_pct(yield_ci)
SOC_ci <- tidy(SOC_boot4,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1:3,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot4,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1:3,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot4,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1:3,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot4,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1:3,3:6]
MWD_output <- make_pct(MWD_ci)

yield_p <- rma.mv(rryield,varyield,data=paddy,mods=~irrigation2-1,random=~1|study/ID)
SOC_p <- rma.mv(rrSOC,varSOC,data=paddy,mods=~irrigation2-1,random=~1|study/ID)
N2O_p <- rma.mv(rrN2O,varN2O,data=paddy,mods=~irrigation2-1,random=~1|study/ID)
CH4_p <- rma.mv(rrCH4,varCH4,data=paddy,mods=~irrigation2-1,random=~1|study/ID)
MWD_p <- rma.mv(rrMWD,varMWD,data=paddy,mods=~irrigation2-1,random=~1|study/ID)

k <- rbind.data.frame(yield_p$pval[1],yield_p$pval[2],yield_p$pval[3],
                      SOC_p$pval[1],SOC_p$pval[2],SOC_p$pval[3],
                      N2O_p$pval[1],N2O_p$pval[2],N2O_p$pval[3],
                      CH4_p$pval[1],CH4_p$pval[2],CH4_p$pval[3],
                      MWD_p$pval[1],MWD_p$pval[2],MWD_p$pval[3])
names(k) <- "pval"
l <- rbind.data.frame(yield_output,SOC_output,N2O_output,CH4_output,MWD_output)
m <- cbind.data.frame(k,l)
m$variable <- c("Yield","Yield","Yield","SOC","SOC","SOC","N2O","N2O","N2O",
                "CH4","CH4","CH4","MWD","MWD","MWD")
m$type <- c("Flooding","Intermittent","Rainfed")
#write.csv(m,"S3-irrigation-paddy.csv")

###Interaction###
interaction<- read.csv("Interaction of CC with management.csv")
interaction$study <- as.factor(interaction$study)
#OF
OF <- interaction%>%filter(management=="Organic fertilizer")
Yield1 <- na.omit(OF[OF$variable=="Yield",])
SOC1 <- na.omit(OF[OF$variable=="SOC",])
N2O1 <- na.omit(OF[OF$variable=="N2O",])
CH41 <- na.omit(OF[OF$variable=="CH4",])
MWD1 <- na.omit(OF[OF$variable=="MWD",])

Yield.func5 <- function(Yield1, index) {
  Yield <- rma.mv(rrinter,varinter,data=Yield1,subset=index,random=~1|study)
  c(Yield$beta, Yield$se^2)}
SOC.func5 <- function(SOC1, index) {
  SOC <- rma.mv(rrinter,varinter,data=SOC1,subset=index,random=~1|study)
  c(SOC$beta, SOC$se^2)}
N2O.func5 <- function(N2O1, index) {
  N2O <- rma.mv(rrinter,varinter,data=N2O1,subset=index,random=~1|study)
  c(N2O$beta, N2O$se^2)}
CH4.func5 <- function(CH41, index) {
  CH4 <- rma.mv(rrinter,varinter,data=CH41,subset=index,random=~1|study)
  c(CH4$beta, CH4$se^2)}
MWD.func5 <- function(MWD1, index) {
  MWD <- rma.mv(rrinter,varinter,data=MWD1,subset=index,random=~1|study)
  c(MWD$beta, MWD$se^2)}

Yield_boot5 <- boot(Yield1, Yield.func5, R=49)
SOC_boot5 <- boot(SOC1, SOC.func5, R=49)
N2O_boot5 <- boot(N2O1, N2O.func5, R=49)
CH4_boot5 <- boot(CH41, CH4.func5, R=49)
MWD_boot5 <- boot(MWD1, MWD.func5, R=49)

Yield_ci <- tidy(Yield_boot5,conf.int=T,conf.method = "perc")
Yield_ci <- data.frame(Yield_ci)
Yield_ci$mean <- (Yield_ci$conf.low+Yield_ci$conf.high)/2
Yield_ci <- Yield_ci[1,3:6]
Yield_output <- make_pct(Yield_ci)
SOC_ci <- tidy(SOC_boot5,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot5,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot5,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot5,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1,3:6]
MWD_output <- make_pct(MWD_ci)

OF_output <- rbind.data.frame(Yield_output,SOC_output,N2O_output,CH4_output,MWD_output)

#Tillage
Tillage <- interaction%>%filter(management=="Conservation tillage")
Yield2 <- na.omit(Tillage[Tillage$variable=="Yield",])
SOC2 <- na.omit(Tillage[Tillage$variable=="SOC",])
N2O2 <- na.omit(Tillage[Tillage$variable=="N2O",])
CH42 <- na.omit(Tillage[Tillage$variable=="CH4",])
MWD2 <- na.omit(Tillage[Tillage$variable=="MWD",])

Yield.func6 <- function(Yield2, index) {
  Yield <- rma.mv(rrinter,varinter,data=Yield2,subset=index,random=~1|study)
  c(Yield$beta, Yield$se^2)}
SOC.func6 <- function(SOC2, index) {
  SOC <- rma.mv(rrinter,varinter,data=SOC2,subset=index,random=~1|study)
  c(SOC$beta, SOC$se^2)}
N2O.func6 <- function(N2O2, index) {
  N2O <- rma.mv(rrinter,varinter,data=N2O2,subset=index,random=~1|study)
  c(N2O$beta, N2O$se^2)}
CH4.func6 <- function(CH42, index) {
  CH4 <- rma.mv(rrinter,varinter,data=CH42,subset=index,random=~1|study)
  c(CH4$beta, CH4$se^2)}
MWD.func6 <- function(MWD2, index) {
  MWD <- rma.mv(rrinter,varinter,data=MWD2,subset=index,random=~1|study)
  c(MWD$beta, MWD$se^2)}

Yield_boot6 <- boot(Yield2, Yield.func6, R=49)
SOC_boot6 <- boot(SOC2, SOC.func6, R=49)
N2O_boot6 <- boot(N2O2, N2O.func6, R=49)
CH4_boot6 <- boot(CH42, CH4.func6, R=49)
MWD_boot6 <- boot(MWD2, MWD.func6, R=49)

Yield_ci <- tidy(Yield_boot6,conf.int=T,conf.method = "perc")
Yield_ci <- data.frame(Yield_ci)
Yield_ci$mean <- (Yield_ci$conf.low+Yield_ci$conf.high)/2
Yield_ci <- Yield_ci[1,3:6]
Yield_output <- make_pct(Yield_ci)
SOC_ci <- tidy(SOC_boot6,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot6,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot6,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1,3:6]
CH4_output <- make_pct(CH4_ci)
MWD_ci <- tidy(MWD_boot6,conf.int=T,conf.method = "perc")
MWD_ci <- data.frame(MWD_ci)
MWD_ci$mean <- (MWD_ci$conf.low+MWD_ci$conf.high)/2
MWD_ci <- MWD_ci[1,3:6]
MWD_output <- make_pct(MWD_ci)

Tillage_output <- rbind.data.frame(Yield_output,SOC_output,N2O_output,CH4_output,MWD_output)

#Irrigation
Irrigation <- interaction%>%filter(management=="Intermittent irrigation")
Yield3 <- na.omit(Irrigation[Irrigation$variable=="Yield",])
SOC3 <- na.omit(Irrigation[Irrigation$variable=="SOC",])
N2O3 <- na.omit(Irrigation[Irrigation$variable=="N2O",])
CH43 <- na.omit(Irrigation[Irrigation$variable=="CH4",])

Yield.func7 <- function(Yield3, index) {
  Yield <- rma.mv(rrinter,varinter,data=Yield3,subset=index,random=~1|study)
  c(Yield$beta, Yield$se^2)}
SOC.func7 <- function(SOC3, index) {
  SOC <- rma(rrinter,varinter,data=SOC3,subset=index)
  c(SOC$beta, SOC$se^2)}
N2O.func7 <- function(N2O3, index) {
  N2O <- rma.mv(rrinter,varinter,data=N2O3,subset=index,random=~1|study)
  c(N2O$beta, N2O$se^2)}
CH4.func7 <- function(CH43, index) {
  CH4 <- rma.mv(rrinter,varinter,data=CH43,subset=index,random=~1|study)
  c(CH4$beta, CH4$se^2)}

Yield_boot7 <- boot(Yield3, Yield.func7, R=49)
SOC_boot7 <- boot(SOC3, SOC.func7, R=49)
N2O_boot7 <- boot(N2O3, N2O.func7, R=49)
CH4_boot7 <- boot(CH43, CH4.func7, R=49)

Yield_ci <- tidy(Yield_boot7,conf.int=T,conf.method = "perc")
Yield_ci <- data.frame(Yield_ci)
Yield_ci$mean <- (Yield_ci$conf.low+Yield_ci$conf.high)/2
Yield_ci <- Yield_ci[1,3:6]
Yield_output <- make_pct(Yield_ci)
SOC_ci <- tidy(SOC_boot7,conf.int=T,conf.method = "perc")
SOC_ci <- data.frame(SOC_ci)
SOC_ci$mean <- (SOC_ci$conf.low+SOC_ci$conf.high)/2
SOC_ci <- SOC_ci[1,3:6]
SOC_output <- make_pct(SOC_ci)
N2O_ci <- tidy(N2O_boot7,conf.int=T,conf.method = "perc")
N2O_ci <- data.frame(N2O_ci)
N2O_ci$mean <- (N2O_ci$conf.low+N2O_ci$conf.high)/2
N2O_ci <- N2O_ci[1,3:6]
N2O_output <- make_pct(N2O_ci)
CH4_ci <- tidy(CH4_boot7,conf.int=T,conf.method = "perc")
CH4_ci <- data.frame(CH4_ci)
CH4_ci$mean <- (CH4_ci$conf.low+CH4_ci$conf.high)/2
CH4_ci <- CH4_ci[1,3:6]
CH4_output <- make_pct(CH4_ci)

Irrigation_output <- rbind.data.frame(Yield_output,SOC_output,N2O_output,CH4_output)

OF_output$management <- "Organic fertilizer"
OF_output$variable <- c("Yield","SOC","N2O","CH4","MWD")
Tillage_output$management <- "Conservation tillage"
Tillage_output$variable <- c("Yield","SOC","N2O","CH4","MWD")
Irrigation_output$management <- "Intermittent irrigation"
Irrigation_output$variable <- c("Yield","SOC","N2O","CH4")

Interaction_total <- rbind.data.frame(OF_output,Tillage_output,Irrigation_output)
write.csv(Interaction_total,"Interaction results.csv")


###Plot###
data <- read.csv("Agricultural management.csv")
data$type <- factor(data$type,levels = c("No","Yes","Conservation","Conventional",
                                         "Rainfed","Irrigated","Intermittent","Flooding"))
data$management <- factor(data$management,levels = c("Organic fertilizer","Tillage",
                                                     "Irrigation (upland)","Irrigation (paddy)"))

Yield <- data%>%filter(variable=="Yield")
p1 <- forestplot(df=Yield,name = type,estimate = mean,se=std.error,colour=variable,pvalue = pval)+
  ggforce::facet_col(~management,scales = "free_y",space="free")+
  scale_color_manual(values = "#e4615d")+
  geom_text(aes(label=n,x=13,y=type))+
  xlab("CC effect on yield (%)")+font+
  theme(legend.position = "none",axis.title.y = element_blank())

SOC <- data%>%filter(variable=="SOC")
p2 <- forestplot(df=SOC,name = type,estimate = mean,se=std.error,colour=variable,pvalue = pval)+
  ggforce::facet_col(~management,scales = "free_y",space="free")+
  scale_color_manual(values = "#fdc58f")+
  geom_text(aes(label=n,x=17.2,y=type))+
  xlab(expression(paste("CC effect on ",SOC[stock]," (%)",sep="")))+font+
  theme(legend.position = "none",axis.title.y = element_blank())

N2O <- data%>%filter(variable=="N2O")
p3 <- forestplot(df=N2O,name = type,estimate = mean,se=std.error,colour=variable,pvalue = pval)+
  ggforce::facet_col(~management,scales = "free_y",space="free")+
  scale_color_manual(values = "#99d9d0")+
  geom_text(aes(label=n,x=90,y=type))+
  xlab(expression(paste("CC effect on ",N[2],"O (%)",sep = "")))+font+
  theme(legend.position = "none",axis.title.y = element_blank())

CH4 <- data%>%filter(variable=="CH4")
p4 <- forestplot(df=CH4,name = type,estimate = mean,se=std.error,colour=variable,pvalue = pval)+
  ggforce::facet_col(~management,scales = "free_y",space="free")+
  scale_color_manual(values = "#95b2d6")+
  geom_text(aes(label=n,x=115,y=type))+
  xlab(expression(paste("CC effect on ",CH[4]," (%)",sep = "")))+font+
  theme(legend.position = "none",axis.title.y = element_blank())

MWD <- data%>%filter(variable=="MWD")
p5 <- forestplot(df=MWD,name = type,estimate = mean,se=std.error,colour=variable,pvalue = pval)+
  ggforce::facet_col(~management,scales = "free_y",space="free")+
  scale_color_manual(values = "#ea9c9d")+
  geom_text(aes(label=n,x=26,y=type))+
  xlab("CC effect on MWD (%)")+font+
  theme(legend.position = "none",axis.title.y = element_blank())

inter <- read.csv("Interaction results.csv")
inter$management <- factor(inter$management,levels = c("Organic fertilizer","Conservation tillage",
                                                       "Intermittent irrigation"))
inter$variable <- factor(inter$variable,levels = c("Yield","SOC","N2O","CH4","MWD"))
inter$sig <- factor(inter$sig,labels = c("Negative","Additive","Positive"))

p6 <- ggplot(inter)+ 
  geom_hline(aes(yintercept=0),linetype="dashed",colour="black",size=0.25)+ ##在0画虚线
  geom_point(aes(x=variable, y=mean,color=sig),stat="identity",size=3)+ ##画数据点
  geom_errorbar(aes(x=variable,ymin=conf.low, ymax=conf.high,color=sig), width=0, size=0.8)+ ##画95%CI
  geom_text(aes(label=n, x=variable,y=mean,hjust=-0.5,vjust=0.15),size=3.8)+  ##添加数据量
  facet_wrap(~management,ncol = 1,scales = "free_y")+
  scale_x_discrete(labels = c("Yield",expression(SOC[stock]),
                              expression(paste(N[2],"O",sep="")),
                              expression(CH[4]),"MWD"))+
  labs(color = "") +ylab("CC x agricultural management effect (%)")+
  theme_cowplot()+font+
  theme(legend.position="top",axis.text.x = element_text(angle=45,hjust = 0.2,vjust = 0.51),
        axis.title.x = element_blank(),
        legend.key = element_rect(fill = NA),legend.box.spacing = unit(0.2,"cm"),
        legend.spacing.x = unit(0.18,"cm"),strip.background = element_rect(fill = "#dbdbdb"))+ 
  scale_colour_manual(values=c("#ff0606","#bdbdbd","#525dff"))

ggdraw()+ 
  draw_plot(p3, x=2/3, y=0.5, width = 1/3, height = 0.5)+
  draw_plot(p2, x=1/3, y=0.5, width = 1/3, height = 0.5)+
  draw_plot(p1, x=0, y=0.5, width = 1/3, height = 0.5)+
  draw_plot(p4, x=0, y=0, width = 1/3, height = 0.5)+
  draw_plot(p5, x=1/3, y=0, width = 1/3, height = 0.5)+
  draw_plot(p6, x=2/3+0.01, y=-0.02, width = 1/3-0.01, height = 0.53)+
  draw_plot_label(label = c("a","b","c","d","e","f"), size = 15,
                  x=c(0,1/3,2/3,0,1/3,2/3), y=c(1,1,1,0.5,0.5,0.5))##11.6*8.6

