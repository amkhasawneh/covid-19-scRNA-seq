#################################Loading########################################
library(shazam)

#Loading the data:
critical119 <- read.table("from_cellranger/critical119/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
critical120 <- read.table("from_cellranger/critical120/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
critical213 <- read.table("from_cellranger/critical213/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
critical238 <- read.table("from_cellranger/critical238/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
critical293 <- read.table("from_cellranger/critical293/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
critical308 <- read.table("from_cellranger/critical308/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
mild186 <- read.table("from_cellranger/mild186/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
mild227 <- read.table("from_cellranger/mild227/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
moderate124 <- read.table("from_cellranger/moderate124/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
moderate138 <- read.table("from_cellranger/moderate138/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
moderate272 <- read.table("from_cellranger/moderate272/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
moderate303 <- read.table("from_cellranger/moderate303/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
severe122 <- read.table("from_cellranger/severe122/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
severe123 <- read.table("from_cellranger/severe123/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
hc1 <- read.table("from_cellranger/healthy1/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
hc2 <- read.table("from_cellranger/healthy2/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)
hc3 <- read.table("from_cellranger/healthy3/vdj_b/filtered_contig_igblast_db-pass.tsv", sep = "\t", header = T)


#Combinging the data frames:
contig.list <- list(critical119 = critical119, critical120 = critical120, critical213 = critical213,
                    critical238 = critical238, critical293 = critical293, critical308 = critical308,
                    mild186 = mild186, mild227 = mild227, 
                    moderate124 = moderate124, moderate138 = moderate138, moderate272 = moderate272, moderate303 = moderate303,
                    severe122 = severe122, severe123 = severe123, 
                    hc1 = hc1, hc2 = hc2, hc3 = hc3)

for (i in names(contig.list)) {
  contig.list[[i]]$sample_id <- i
}


#################################Shazam#########################################
# Group cells in a one-stage process (VJthenLen=FALSE) and using
# both heavy and light chain sequences (onlyHeavy=FALSE)
for (i in 1:length(contig.list)) {
  dist_sc[i] <- distToNearest(contig.list[[i]], cellIdColumn="cell_id", locusColumn="locus", 
                              VJthenLen=FALSE, onlyHeavy=FALSE)
  
}
dist_sc2 <- distToNearest(critical120, cellIdColumn="cell_id", locusColumn="locus", 
                          VJthenLen=FALSE, onlyHeavy=FALSE)





