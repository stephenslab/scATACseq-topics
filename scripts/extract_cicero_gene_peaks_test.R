# Using the Cicero co-accessibility data, generate a data structure
# containing a list of peaks for each gene.
library(tools)
library(dplyr)

# Get the peaks linked to each gene.
data.dir <- "/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data"
cicero <- readRDS(file.paxth(data.dir, "master_cicero_conns.rds"))

cicero <- cicero %>%
  filter(cluster %in% c(11,18,22,25)) %>%
  select(Peak1, Peak2, peak1.tss.gene_id, peak2.tss.gene_id, coaccess,
         peak1.isproximal, peak2.isproximal, conn_type, cluster)

cicero  <- transform(cicero,
                     Peak1             = as.character(Peak1),
                     Peak2             = as.character(Peak2),
                     peak1.tss.gene_id = factor(peak1.tss.gene_id),
                     peak2.tss.gene_id = factor(peak2.tss.gene_id))

table(cicero$conn_type)
head(cicero[cicero$conn_type == "proximal_proximal", ])
head(cicero[cicero$conn_type == "distal_proximal", ])
head(cicero[cicero$conn_type == "distal_distal", ])

# the connections should be symmetric
peak1.genesIDs <- unique(cicero$peak1.tss.gene_id)
peak2.genesIDs <- unique(cicero$peak2.tss.gene_id)
setequal(peak1.genesIDs, peak2.genesIDs)

# example gene
gene.of.interest <- "ENSMUSG00000000058"
# gene.of.interest <- "ENSMUSG00000006732"

cicero.test1 <- cicero %>% filter(peak1.tss.gene_id %in% gene.of.interest)

table(cicero.test1$conn_type)

cicero.test1 %>% filter(peak1.tss.gene_id %in% gene.of.interest) %>% head()

cicero_gene_self_peaks <- cicero.test1 %>% filter(peak1.tss.gene_id == gene.of.interest) %>% pull(Peak1) %>% unique()

cicero_gene_connected_peaks <- cicero.test1 %>% filter(peak1.tss.gene_id == gene.of.interest) %>% pull(Peak2) %>% unique()

cicero_gene_all_peaks <- union(cicero_gene_self_peaks, cicero_gene_connected_peaks)

setdiff(cicero_gene_all_peaks, cicero_gene_connected_peaks)

# all genes
cicero_gene_self_peaks <- tapply(cicero.test1$Peak1,cicero.test1$peak1.tss.gene_id,unique,
                                 simplify = FALSE)

cicero_gene_connected_peaks <- tapply(cicero.test1$Peak2,cicero.test1$peak1.tss.gene_id,unique,
                                      simplify = FALSE)

cicero_gene_self_peaks[[gene.of.interest]]

cicero_gene_connected_peaks[[gene.of.interest]]
