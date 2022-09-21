# Using the Cicero co-accessibility data, generate a data structure
# containing a list of peaks for each gene.
library(tools)

# Get the peaks linked to each gene.
cicero <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
class(cicero) <- "data.frame"
cicero <- cicero[c("Peak1","Peak2","peak1.tss.gene_id",
                   "peak2.tss.gene_id","cluster")]
cicero <- subset(cicero,is.element(cluster,c(11,18,22,25)))
cicero  <- transform(cicero,
                     Peak1             = as.character(Peak1),
                     Peak2             = as.character(Peak2),
                     peak1.tss.gene_id = factor(peak1.tss.gene_id),
                     peak2.tss.gene_id = factor(peak2.tss.gene_id))
# cicero_gene <- tapply(cicero$Peak1,cicero$peak1.tss.gene_id,unique,
#                       simplify = FALSE)
cicero_gene <- tapply(cicero$Peak2,cicero$peak1.tss.gene_id,unique,
                      simplify = FALSE)

# Save these data to an .RData file.
save(list = "cicero_gene",file = "cicero_gene.RData")
resaveRdaFiles("cicero_gene.RData")
