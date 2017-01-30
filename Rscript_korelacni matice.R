
trw <- readXL("F:/zbytka/cluster_oprava/cluster_oprava.xlsx",rownames=FALSE, header=TRUE, na="", sheet="trw_chr_R",stringsAsFactors=TRUE)
vla <- readXL("F:/zbytka/cluster_oprava/cluster_oprava.xlsx",rownames=FALSE, header=TRUE, na="", sheet="vla_chr_R",stringsAsFactors=TRUE)

cor.trw <- data.frame(cor(trw[,c("ZB_ALL","MO_ALL","HO_ALL","ZE_ALL","ZB_POS","ZB_NON","ZB_NEG","MO_POS","MO_NON","MO_NEG","HO_POS","HO_NON","HO_NEG","ZE_POS","ZE_NON","ZE_NEG")]))
cor.vla <- data.frame(cor(vla[,c("ZB_ALL","MO_ALL","HO_ALL","ZE_ALL","ZB_POS","ZB_NON","ZB_NEG","MO_POS","MO_NON","MO_NEG","HO_POS","HO_NON","HO_NEG","ZE_POS","ZE_NON","ZE_NEG")]))

write.table(cor.vla, "F:/cor.vla.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

rcorr.adjust(trw[,c("ZB_ALL","MO_ALL","HO_ALL","ZE_ALL","ZB_POS","ZB_NON","ZB_NEG","MO_POS","MO_NON","MO_NEG","HO_POS","HO_NON","HO_NEG","ZE_POS","ZE_NON","ZE_NEG")], type="pearson", use="complete")
rcorr.adjust(vla[,c("ZB_ALL","MO_ALL","HO_ALL","ZE_ALL","ZB_POS","ZB_NON","ZB_NEG","MO_POS","MO_NON","MO_NEG","HO_POS","HO_NON","HO_NEG","ZE_POS","ZE_NON","ZE_NEG")], type="pearson", use="complete")

