library(dplR)
library(bootRes)
library(pointRes)

scPDSI <- read.table("F:/zbytka/DendroClim/scPDSI.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
TEMP <- read.table("F:/zbytka/DendroClim/TEMP.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

   # Jednotlive serie
	serieTRW <- read.rwl("F:/zbytka/serie/trw.rwl", format="auto")
	serieVLA <- read.rwl("F:/zbytka/serie/vla_sehrane.rwl", format="auto")

climate <- list (TEMP, scPDSI); timespan.row.trw<-c(178:235); timespan.row.vla<-c(161:218)

# Common interval - pouzivam stromy pouze mezi 60 a 150 lety
common.interval(serieTRW, type="both")
common.interval(serieVLA, type="both")

AGE_COHORT_TRW<-(intersect(which(rwl.stats(serieTRW)
                                         [,4]<150),which(rwl.stats(serieTRW)[,4]>60)))

AGE_COHORT_VLA<-(intersect(which(rwl.stats(serieVLA)
                                         [,4]<150),which(rwl.stats(serieVLA)[,4]>60)))

Zernov_trw <- serieTRW[c(intersect(1:29, AGE_COHORT_TRW))]; Houkvice_trw <- serieTRW[c(intersect(30:51, AGE_COHORT_TRW))]; Zbytka_trw <- serieTRW[c(intersect(52:73, AGE_COHORT_TRW))]; RefOrig_trw <- serieTRW[c(intersect(74:92, AGE_COHORT_TRW))]
Zernov_vla <- serieVLA[c(intersect(1:29, AGE_COHORT_VLA))]; Houkvice_vla <- serieVLA[c(intersect(30:51, AGE_COHORT_VLA))]; Zbytka_vla <- serieVLA[c(intersect(52:73, AGE_COHORT_VLA))]; RefOrig_vla <- serieVLA[c(intersect(74:92, AGE_COHORT_VLA))]
serieTRW<-serieTRW[c(AGE_COHORT_TRW)]
serieVLA<-serieVLA[c(AGE_COHORT_VLA)]


rwi.stats(Zbytka_trw[(178:235),]);rwi.stats(RefOrig_trw[(178:235),]);rwi.stats(Houkvice_trw[(178:235),]);rwi.stats(Zernov_trw[(178:235),])
rwi.stats(Zernov_vla[(161:218),]);rwi.stats(Houkvice_vla[(161:218),]);rwi.stats(RefOrig_vla[(161:218),]);rwi.stats(Zbytka_vla[(161:218),])

###################
##### Vypocet #####
###################

chron.cor(Zbytka_trw, timespan.row.trw)
chron.cor(RefOrig_trw, timespan.row.trw)
chron.cor(Houkvice_trw, timespan.row.trw)
chron.cor(Zernov_trw, timespan.row.trw)

chron.cor(Zbytka_vla, timespan.row.vla)
chron.cor(RefOrig_vla, timespan.row.vla)
chron.cor(Houkvice_vla, timespan.row.vla)
chron.cor(Zernov_vla, timespan.row.vla)

#####################################################
# Moje funkce - vypocet response chronologii + jejich korelaci s klimatem
#####################################################

chron.cor<-function(analysis, timespan) {

PCGA.an<-PCGA(analysis[timespan,]) # Vypocet PCGA

korPCA.an <- korPCA(PCGA.an, axes(analysis,timespan)$PC1, axes(analysis,timespan)$PC2) # Prvni dve komponenty vypocteny pomoci uzivatelske funkce axes

POS_Chron<-chron(analysis[timespan,which(korPCA.an$korelacePC2=="POS")]) # Vytvoreni response chronologii
NON_Chron<-chron(analysis[timespan,which(korPCA.an$korelacePC2=="NON")])
NEG_Chron<-chron(analysis[timespan,which(korPCA.an$korelacePC2=="NEG")])

#Chrons<-cbind(POS_Chron, NON_Chron, NEG_Chron, chron(analysis[timespan,])) # Vytvoreni vystupniho souboru chronologii
#colnames(Chrons)<-c("POS", "POSsd", "NON", "NONsd", "NEG", "NEGsd", "ALL", "ALLsd")

#print(rwi.stats(analysis[timespan,which(korPCA.an$korelacePC2=="POS")]))# EPS/Rbar
#print(rwi.stats(analysis[timespan,which(korPCA.an$korelacePC2=="NEG")]))
#print(rwi.stats(analysis[timespan,which(korPCA.an$korelacePC2=="NON")]))
#print(rwi.stats(analysis[timespan,]))

#POS.rbar<-rwi.stats.running(analysis[timespan,which(korPCA.an$korelacePC2=="POS")], window.length=21, window.overlap=20) # Klouzavy Rbar
#NEG.rbar<-rwi.stats.running(analysis[timespan,which(korPCA.an$korelacePC2=="NEG")], window.length=21, window.overlap=20)
#NON.rbar<-rwi.stats.running(analysis[timespan,which(korPCA.an$korelacePC2=="NON")], window.length=21, window.overlap=20)
#ALL.Rbar<-rwi.stats.running(analysis[timespan,], window.length=21, window.overlap=20)

#Rbar<-cbind(POS.rbar[,c(1:3,14,15)], NON.rbar[,c(14,15)], NEG.rbar[,c(14,15)], ALL.Rbar[,c(14,15)]) # Vystupni soubor s klouzavym Rbar a EPS
#colnames(Rbar)<-c("start.y", "mid.y", "end.y", "Rbar_POS", "EPS_POS", "Rbar_NON", "EPS_NON", "Rbar_NEG", "EPS_NEG", "Rbar_ALL", "EPS_ALL")

#point.POS <- pointer.norm(analysis[timespan,which(korPCA.an$korelacePC2=="POS")], method.thresh="Cropper")$out # Pointer years - Cropper values s defaultnim nastavenim
#point.NEG <- pointer.norm(analysis[timespan,which(korPCA.an$korelacePC2=="NEG")], method.thresh="Cropper")$out
#point.NON <- pointer.norm(analysis[timespan,which(korPCA.an$korelacePC2=="NON")], method.thresh="Cropper")$out
#point.ALL <- pointer.norm(analysis[timespan,], method.thresh="Cropper")$out
#Pointers <- cbind(point.POS[,c(1,5:6)], point.NEG[,c(5:6)], point.NON[,c(5:6)], point.ALL[,c(5:6)])
#colnames(Pointers)<-c("year", "POSnature", "POSCropp", "NEGnature", "NEGCropp", "NONnature", "NONCropp", "ALLnature", "ALLCropp")

#dc.POS<-dcc(POS_Chron, climate, method="correlation") # Vypocet korelacnich koeficientu
#dc.NON<-dcc(NON_Chron, climate, method="correlation")
#dc.NEG<-dcc(NEG_Chron, climate, method="correlation")
#dc.ALL<-dcc(chron(analysis[timespan,]), climate, method="correlation")

#PC1.table <- data.frame(axes(analysis,timespan)$PC1); rownames(PC1.table) <- PCGA.an$period # Vypocet klimatickeho signal PC1 a PC2
#PC2.table <- data.frame(axes(analysis,timespan)$PC2); rownames(PC2.table) <- PCGA.an$period
#dc.PC1 <- dcc(chron(PC1.table), climate, method="correlation")
#dc.PC2 <- dcc(chron(PC2.table), climate, method="correlation")
#ChronsPC <- cbind(chron(PC1.table), chron(PC2.table)); colnames(ChronsPC)<-c("PC1", "PC1sd", "PC2", "PC2sd")

mdc.POS.t<-mdcc(POS_Chron, TEMP, method = "correlation", start=-6, end=9, boot=T, vnames = c("T"), win.size=21, win.offset=1, startlast=TRUE) # Klouzave korelace
mdc.NON.t<-mdcc(NON_Chron, TEMP, method = "correlation", start=-6, end=9, boot=T, vnames = c("T"), win.size=21, win.offset=1, startlast=TRUE)
mdc.NEG.t<-mdcc(NEG_Chron, TEMP, method = "correlation", start=-6, end=9, boot=T, vnames = c("T"), win.size=21, win.offset=1, startlast=TRUE)
mdc.ALL.t<-mdcc(chron(analysis[timespan,]), TEMP, method = "correlation", start=-6, end=9, boot=T, vnames = c("T"), win.size=21, win.offset=1, startlast=TRUE)
mdc.POS.p<-mdcc(POS_Chron, scPDSI, method = "correlation", start=-6, end=9, boot=T, vnames = c("scPDSI"), win.size=21, win.offset=1, startlast=TRUE) # Klouzave korelace
mdc.NON.p<-mdcc(NON_Chron, scPDSI, method = "correlation", start=-6, end=9, boot=T, vnames = c("scPDSI"), win.size=21, win.offset=1, startlast=TRUE)
mdc.NEG.p<-mdcc(NEG_Chron, scPDSI, method = "correlation", start=-6, end=9, boot=T, vnames = c("scPDSI"), win.size=21, win.offset=1, startlast=TRUE)
mdc.ALL.p<-mdcc(chron(analysis[timespan,]), scPDSI, method = "correlation", start=-6, end=9, boot=T, vnames = c("scPDSI"), win.size=21, win.offset=1, startlast=TRUE)

#write.table(Chrons, "F:/zbytka/cluster_oprava/chron_kor/chron.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(Rbar, "F:/zbytka/cluster_oprava/chron_kor/rbar.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.POS, "F:/zbytka/cluster_oprava/chron_kor/dcPOS.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.NON, "F:/zbytka/cluster_oprava/chron_kor/dcNON.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.NEG, "F:/zbytka/cluster_oprava/chron_kor/dcNEG.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.ALL, "F:/zbytka/cluster_oprava/chron_kor/dcALL.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.PC1, "F:/zbytka/cluster_oprava/chron_kor/dcPC1.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(dc.PC2, "F:/zbytka/cluster_oprava/chron_kor/dcPC2.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(ChronsPC, "F:/zbytka/cluster_oprava/chron_kor/chronPC.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#write.table(Pointers, "F:/zbytka/cluster_oprava/chron_kor/point.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")


write.table(mdc.POS.t, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcPOS.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.NON.t, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcNON.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.NEG.t, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcNEG.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.ALL.t, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcALL.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

write.table(mdc.POS.p, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcPOS.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.NON.p, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcNON.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.NEG.p, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcNEG.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(mdc.ALL.p, "F:/zbytka/cluster_oprava/chron_kor/mdc/mdcALL.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
}

######################################
### Pomocna funkce pro vypocet os PCA
######################################


axes <- function(analysis, timespan) {
loadings <- (PCGA(analysis[timespan,])$pca$rot)[,c(1:2)]
PC1 <- 0; PC2 <- 0
	for (i in c(1:nrow(loadings)))
		{PC1 <- PC1 + (loadings[i,1]*analysis[timespan,i])
		PC2 <- PC2 + (loadings[i,2]*analysis[timespan,i])}
return(list(PC1=PC1, PC2=PC2))
}

####################################
# Pomocna funkce pro vypocet korelaci mezi seriemi a osou + jejich roztrideni
####################################

korPCA <- function(vystupPCGA, PC1, PC2) {

korelacePC1<-vector(mode="numeric")
korelacePC2<-vector(mode="numeric")
korelacePC<-vector(mode="numeric")

for(i in 1:ncol(vystupPCGA$pop))
  {
    korPC1<-lm(vystupPCGA$pop[,i]~PC1)
    korPC2<-lm(vystupPCGA$pop[,i]~PC2)
	if(summary(korPC1)$coef[2,4]<0.05&&sign(summary(korPC1)$coef[2,1])<0)
    		{
      			korelacePC1[i]<-"NEG" }
	if(summary(korPC1)$coef[2,4]<0.05&&sign(summary(korPC1)$coef[2,1])>0)
    		{
      			korelacePC1[i]<-"POS" }
	if(summary(korPC1)$coef[2,4]>0.05)
    		{
      			korelacePC1[i]<-"NON" }

	if(summary(korPC2)$coef[2,4]<0.05&&sign(summary(korPC2)$coef[2,1])<0)
    		{
      			korelacePC2[i]<-"NEG" }
	if(summary(korPC2)$coef[2,4]<0.05&&sign(summary(korPC2)$coef[2,1])>0)
    		{
      			korelacePC2[i]<-"POS" }
	if(summary(korPC2)$coef[2,4]>0.05)
    		{
      			korelacePC2[i]<-"NON" }

}
	korelacePC[which(korelacePC2=="NEG")]<-"red" # Jednodussi klasifikace, pripadne nutne vypoznamkovat
	korelacePC[which(korelacePC2=="NON")]<-"grey"
	korelacePC[which(korelacePC2=="POS")]<-"blue"
		
	plot.PCGA(vystupPCGA, col.vec=korelacePC)

return(list(korelacePC1=korelacePC1, korelacePC2=korelacePC2))
}
