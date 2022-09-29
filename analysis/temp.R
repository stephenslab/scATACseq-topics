k <- 9 # 1, 2, 3, 5, 8, 9, 10
prior_genes <- c(pt_genes_S1)
# genes[!(abs((de_gene$z[,k]) > 50 & de_gene$postmean[,k] > 1) |
#         de_gene$postmean[,k] > 3)] <- ""
genes[!(abs((de_gene$z[,k]) > 50 & de_gene$postmean[,k] > 2) |
        de_gene$postmean[,k] > 2.75 |
        is.element(genes,prior_genes))] <- ""
# i <- which(is.element(genes,loh_genes))
# i <- which(is.element(genes,dct_genes))
# i <- which(is.element(genes,podo_genes))
# i <- which(is.element(genes,endo_genes))
# i <- which(is.element(genes,cd_genes))
i <- which(is.element(genes,prior_genes))
genes[i] <- paste0("*",genes[i])
p <- volcano_plot(de_gene,k = k,labels = genes,ymax = 300,
                  do.label = function (lfc, z) TRUE) +
  labs(x = "mean coef",y = "logBF")
