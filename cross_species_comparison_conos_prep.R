library(conos)
library(pagoda2)
library(Matrix)
library(Seurat)
#source("/home/qiwenhu/hubmap/pip/scripts/active_script/processing.func.R")
setwd("/home/qiwenhu/schizo/data/BICCN_data")
dir <- "/home/qiwenhu/schizo/data/BICCN_data"

human.p2 <- readRDS(file.path(dir, "human.p2.rds"))
mouse.p2 <- readRDS(file.path(dir, "mouse.p2.rds"))
marmoset.p2 <- readRDS(file.path(dir, "marmoset.p2.rds"))

integrated.panel <- list(human.p2, mouse.p2, marmoset.p2)
names(integrated.panel) <- c("human", "mouse", "marmoset")

# conos integration
con <- Conos$new(integrated.panel, n.cores=20)
con$buildGraph(k=30, k.self=5, space='CCA', ncomps=30, 
               n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
con$findCommunities(method=leiden.community, resolution=1)
#con$embedGraph(method="UMAP")
con$embedGraph()
saveRDS(con, file.path(dir, "human.mouse.marmo.conos.cca.rds"))
