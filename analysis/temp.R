# Draft analysis of motifs in topics k = 1 through 10.
homer <- readRDS(file.path("output/Buenrostro_2018",
                           "binarized/filtered_peaks",
                           "homer-buenrostro2018-k=10-noshrink.rds"))

# Compile the motif enrichment p-values into a single table, after
# removing duplicate results.
motifs <- sort(unique(homer$k1[,"Motif Name"]))
n      <- length(motifs)
pvals  <- matrix(0,n,k)
rownames(pvals) <- motifs
colnames(pvals) <- paste0("k",1:k)
for (i in 1:k) {
  dat  <- homer[[i]]
  rows <- which(!duplicated(dat[,"Motif Name"]))
  dat  <- dat[rows,]
  rownames(dat) <- dat[,"Motif Name"]
  pvals[,i] <- dat[motifs,"P-value"]
}
pvals <- as.data.frame(pvals)

# TO DO: Create CSV file containing the motif enrichment results.

# Topics to do still: 10.

# Compute the log p-values.
lpval <- log10(pvals + 1e-256)
i <- c(10,10)
y <- apply(lpval[,i],1,min) - apply(lpval[,-i],1,min)
head(sort(y),n = 20)

motifs <- c("Hoxc9(Homeobox)/Ainv15-Hoxc9-ChIP-Seq(GSE21812)/Homer",
            "Hoxb4(Homeobox)/ES-Hoxb4-ChIP-Seq(GSE34014)/Homer",
            "HOXA1(Homeobox)/mES-Hoxa1-ChIP-Seq(SRP084292)/Homer",
            "ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer",
            "EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Homer",
            "Atf1(bZIP)/K562-ATF1-ChIP-Seq(GSE31477)/Homer",
            "Atf4(bZIP)/MEF-Atf4-ChIP-Seq(GSE35681)/Homer",
            "Atf7(bZIP)/3T3L1-Atf7-ChIP-Seq(GSE56872)/Homer",
            "CEBP(bZIP)/ThioMac-CEBPb-ChIP-Seq(GSE21512)/Homer",
            "CEBP:AP1(bZIP)/ThioMac-CEBPb-ChIP-Seq(GSE21512)/Homer",
            "CEBP:CEBP(bZIP)/MEF-Chop-ChIP-Seq(GSE35681)/Homer",
            "EBF2(EBF)/BrownAdipose-EBF2-ChIP-Seq(GSE97114)/Homer",
            "EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer",
            "EBF1(EBF)/Near-E2A-ChIP-Seq(GSE21512)/Homer",
            "ETS:E-box(ETS,bHLH)/HPC7-Scl-ChIP-Seq(GSE22178)/Homer",
            "GATA:SCL(Zf,bHLH)/Ter119-SCL-ChIP-Seq(GSE18720)/Homer",
            "Gata1(Zf)/K562-GATA1-ChIP-Seq(GSE18829)/Homer",
            "Gata2(Zf)/K562-GATA2-ChIP-Seq(GSE18829)/Homer",
            "GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer",
            "Gata4(Zf)/Heart-Gata4-ChIP-Seq(GSE35151)/Homer",
            "Gata6(Zf)/HUG1N-GATA6-ChIP-Seq(GSE51936)/Homer",
            "Tcf12(bHLH)/GM12878-Tcf12-ChIP-Seq(GSE32465)/Homer",
            "TCF4(bHLH)/SHSY5Y-TCF4-ChIP-Seq(GSE96915)/Homer",
            "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
            "Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer",
            "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
            "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer")
lpval <- lpval[motifs,]

# Plot the p-values in a tile plot.
# Colors from colorbrewer2.org.
colors <- c("#d73027","#f46d43","#fdae61","#fee090",
            "#e0f3f8","#abd9e9","#74add1","#4575b4")
pdat <- NULL
for (i in 1:10) {
  pdat <- rbind(pdat,
                data.frame(motif = motifs,
                           topic = paste0("k",i),
                           lpval = lpval[,i]))
}
pdat <-
  transform(pdat,
            topic = factor(topic,c("k9","k2","k3","k4","k8","k1","k7",
                                   "k5","k6","k10")),
            lpval = cut(lpval,c(-257,-150,-100,-50,-30,-20,-10,-5,0)),
            motif = factor(motif,rev(motifs)))
p <- ggplot(pdat,aes(x = topic,y = motif,fill = lpval)) +
  geom_tile(color = "white",size = 0.5) +
  scale_fill_manual(values = colors) +
  labs(x = "",y = "") +
  theme_cowplot(font_size = 8)
