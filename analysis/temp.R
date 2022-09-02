library(fastTopics)
library(ggplot2)
library(cowplot)
loh_genes <- c("Slc12a1","Umod","Egf","Wfdc15b","Mt2","Ppp1r1a","Mt1",
               "Sostdc1","Slc5a3","Ly6a","Mrps6","Ckb","Ppp1r1b","Cldn10",
               "Cldn19","Ptger3","Gadd45g","Fabp3","Cldn16","Slc5a1","Irx1",
               "Tmem207","Baiap2l2","Shd","Nccrp1","Cyfip2","Clcnka","Gcgr",
               "Cgnl1","Me3","Irx2","Aktip","Fgf9")
dct_genes <- c("Slc12a3","Calb1","Wnk1","Pvalb","Pgam2","Wnk4","Sgms2",
               "Slc16a7","Lhx1","Abca13","Hoxb5os","Emx1","Trpm6","Cwh43",
               "mt-Co1","Oxct1","Trpm7","Papss1","Uroc1","Tfrc","Lhx1os",
               "Tsc22d2","Larp1b","Gm15848")
load(file.path("output/Cusanovich_2018/tissues",
               "de-gene-cusanovich2018-kidney-k=10.RData"))
n             <- nrow(de_gene$postmean)
de_gene$f0    <- rep(0,n)
de_gene$lower <- with(de_gene,postmean - se)
de_gene$upper <- with(de_gene,postmean + se)
genes <- rownames(de_gene$postmean)
genes[!is.element(genes,loh_genes)] <- ""
p1 <- volcano_plot(de_gene,k = 8,labels = genes,ymax = 20,
                   do.label = function (lfc, z) TRUE)
genes <- rownames(de_gene$postmean)
genes[!is.element(genes,dct_genes)] <- ""
p2 <- volcano_plot(de_gene,k = 5,labels = genes,ymax = 20,
                   do.label = function (lfc, z) TRUE)

