Notebook 3: snRNAseq DE analysis, GO analysis
================

<style type="text/css">
.height {
  white-space: pre !important;
  overflow-y: scroll !important;
  height: 50vh !important;
}
</style>

### 1\. Estimating DE genes (WALD test)

``` r
#using same cacoa objects from Notebook2
de.high  <- cao$estimatePerCellTypeDE(test='DESeq2.Wald', n.cores = 20, verbose = TRUE, max.cell.count=80)
de.med  <- caom$estimatePerCellTypeDE(test='DESeq2.Wald', n.cores = 20, verbose = TRUE, max.cell.count=80)
de <- list(high=de.high,med=de.med)
```

### 2\. DE genes bootstrap

``` r
de.high  <- cao$estimatePerCellTypeDE(max.cell.count=80,test='DESeq2.Wald', resampling.method = "bootstrap", n.cores = 30, verbose = TRUE, max.resamplings = 30)
de.med  <- caom$estimatePerCellTypeDE(max.cell.count=80,test='DESeq2.Wald', resampling.method = "bootstrap", n.cores = 30, verbose = TRUE, , max.resamplings = 30)
del.fixedWboot <- list(high=de.high,med=de.med)
```

``` r
#Figure 5.B
boot <- lapply(del.fixedWboot$high, function(x){
  de <- data.table::rbindlist(x$subsamples, idcol = "boot")
  return(de)
})

boot0 <- lapply(boot, function(x){
  boot0 <- x[pvalue < 1e-3, .(nde = .N), by = .(boot)]
  boot0$nde <- boot0$nde/sum(boot0$nde)
  return(boot0)
})

boot1 <- data.table::rbindlist(boot0, idcol = "cell")

#boxplot figure
Fig5.B <- ggplot(boot1,aes(x=reorder(cell,nde,median),y=nde,fill=cell)) + 
    geom_boxplot(notch=T,outlier.shape=NA) + 
    geom_jitter(position=position_jitter(0.1), show.legend=FALSE,alpha=0.1) + 
    coord_cartesian(ylim=c(0, 0.1))+
    #coord_cartesian(ylim=c(0, 1e3))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12),axis.title.x = element_blank()) +ylab("fraction of all DE genes") + scale_fill_manual(values=high.pal) + guides(fill=F)

#layer plot
mv <- tapply(Fig5.B$data$nde, Fig5.B$data$cell, median);
df <- ldf_high; df$ndem <- mv[as.character(df$Subtypes)]; df <- na.omit(df);
Fig5.BL <- ggplot(df, aes(x = reorder(Subtypes,ndem,mean), y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.15, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = high.pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "Helvetica"), axis.title.x = element_blank())
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig5a.jpg" width="60%" style="display: block; margin: auto;" />

### 3\. Calculate GO terms

``` r
library(tibble)
library(dplyr)
library(cowplot)
library(readr)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)
library(pbapply)
require(ComplexHeatmap)
require(circlize)
library(magrittr)
```

#### a)prepare GO categories

``` r
go_datas <- c("BP", "CC", "MF") %>% setNames(., .) %>%
   pblapply(function(n) clusterProfiler:::get_GO_data(org.Hs.eg.db, n, "ENTREZID") %>%
              as.list() %>% as.environment()) # otherwise it pass reference to the environment content
```

#### b)Run enrichment tests

``` r height
enrichGOOpt <- function (gene, OrgDB, goData, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2, minGSSize = 10,
                         maxGSSize = 500, readable = FALSE, pool = FALSE) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))
  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = goData)
  if (is.null(res))
    return(res)
  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(OrgDB)
  if (readable) {
    res <- DOSE::setReadable(res, OrgDB)
  }
  res@ontology <- ont
  return(res)
}

distanceBetweenTerms <- function(go.df) {
  genes.per.go <- sapply(go.df$geneID, strsplit, "/") %>% setNames(go.df$Description)
  all.go.genes <- unique(unlist(genes.per.go))
  all.gos <- unique(go.df$Description)
  genes.per.go.mat <- matrix(0, length(all.go.genes), length(all.gos)) %>%
    `colnames<-`(all.gos) %>% `rownames<-`(all.go.genes)
  for (i in 1:length(genes.per.go)) {
    genes.per.go.mat[genes.per.go[[i]], go.df$Description[[i]]] <- 1
  }
  return(dist(t(genes.per.go.mat), method="binary"))
}

calculate.gos <- function(de, go.datas, n.top.genes=300,n.cores=1) {
  de <- de[unlist(lapply(de,is.list))]

  # add Z scores
  de <- lapply(de,function(d) {
    res.table <- d$res;
    res.table$Z <- -qnorm(res.table$pval/2)
    res.table$Z[is.na(res.table$Z)] <- 0
    res.table$Za <- -qnorm(res.table$padj/2)
    res.table$Za[is.na(res.table$Za)] <- 0
    res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
    res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
    d$res <- res.table;
    d
  })
  
    gns <- list(down=lapply(de,function(x) rownames(x$res)[order(x$res$Z,decreasing=F)[1:n.top.genes]]),
              up=lapply(de,function(x) rownames(x$res)[order(x$res$Z,decreasing=T)[1:n.top.genes]]),
              all=lapply(de,function(x) rownames(x$res)))

    
    return(gns)
   
}
```

#### c)calculate and plot clusters

``` r
goslW <- lapply(de,calculate.gos,n.top.genes=500)
library(circlize)
```

high anno

``` r
g <- lapply(goslW$high, gj)

g$up <- g$up[g$up$pvalue < 0.001,]
g$up$p.adjust2 <- g$up$pvalue %>% p.adjust("bonferroni")
guc <- g$up[g$up$p.adjust2 < 0.05,]


guclusts <- distanceBetweenTerms(guc) %>%  hclust(method='ward.D2') %>% cutree(20) 
gupc <- split(names(guclusts), guclusts)
nguc_per_clust <- sapply(gupc, length)


guc %<>% mutate(GOClust=guclusts[Description])
guname_per_clust <- guc %>% group_by(GOClust, Description) %>% summarise(pvalue=exp(mean(log(pvalue)))) %>% 
    split(.$GOClust) %>% sapply(function(df) df$Description[which.min(df$pvalue)])
guc %<>% mutate(GOClustName=guname_per_clust[as.character(GOClust)])


guc_bp_summ_df <- guc %>% group_by(Type, GOClustName) %>% 
    summarise(p.adjust=min(p.adjust2)) %>% ungroup() %>% mutate(p.adjust=-log10(p.adjust)) %>% 
    tidyr::spread(Type, p.adjust) %>% as.data.frame() %>% set_rownames(.$GOClustName) %>% .[, 2:ncol(.)] #%>% .[, type_order[type_order %in% colnames(.)]]
  guc_bp_summ_df[is.na(guc_bp_summ_df)] <- 0
```

``` r
g$down <- g$down[g$down$pvalue < 0.001,]
g$down$p.adjust2 <- g$down$pvalue %>% p.adjust("bonferroni")
gdc <- g$down[g$down$p.adjust2 < 0.05,]

gdclusts <- distanceBetweenTerms(gdc) %>%  hclust(method='ward.D2') %>% cutree(20) 
gdpc <- split(names(gdclusts), gdclusts)
ngdc_per_clust <- sapply(gdpc, length)



gdc %<>% mutate(GOClust=gdclusts[Description])
gdname_per_clust <- gdc %>% group_by(GOClust, Description) %>% summarise(pvalue=exp(mean(log(pvalue)))) %>% 
    split(.$GOClust) %>% sapply(function(df) df$Description[which.min(df$pvalue)])
gdc %<>% mutate(GOClustName=gdname_per_clust[as.character(GOClust)])


gdc_bp_summ_df <- gdc %>% group_by(Type, GOClustName) %>% 
    summarise(p.adjust=min(p.adjust2)) %>% ungroup() %>% mutate(p.adjust=-log10(p.adjust)) %>% 
    tidyr::spread(Type, p.adjust) %>% as.data.frame() %>% set_rownames(.$GOClustName) %>% .[, 2:ncol(.)] #%>% .[, type_order[type_order %in% colnames(.)]]
  gdc_bp_summ_df[is.na(gdc_bp_summ_df)] <- 0
```

``` r
cols <- list(up=colorRamp2(c(0, 4), c("grey98", "red")),down=colorRamp2(c(0, 4), c("grey98", "blue")))
n.clusters <- 20; max.pval <- 0.05;
```

``` r
hm3 <- Heatmap(as.matrix(guc_bp_summ_df),
             col=cols$up,
              border=T,
              show_row_dend=F,
              show_column_dend=F, 
              heatmap_legend_param = list(title = '-log10(adj.p)'), 
              row_names_max_width = unit(13, "cm"),
              row_names_gp = gpar(fontsize = 10), 
              column_names_max_height = unit(8, "cm"))
```

``` r
hm4 <- Heatmap(as.matrix(gdc_bp_summ_df),
             col=cols$down,
              border=T,
              show_row_dend=F,
              show_column_dend=F, 
              heatmap_legend_param = list(title = '-log10(adj.p)'), 
              row_names_max_width = unit(11, "cm"),
              row_names_gp = gpar(fontsize = 10), 
              column_names_max_height = unit(8, "cm"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/goup.jpg" width="50%" style="display: block; margin: auto;" />

<img src="C:/Users/Katarina/Desktop/scznotebooks/godown.jpg" width="50%" style="display: block; margin: auto;" />

#### d) GO term similarity

``` r
#function used:
library(tidyverse)
library(plyr)

plotCustomSimilarities <- function(gosdt){
gosdt <- data.table(gosdt)
pathway.df <- gosdt[,c("Description", "Type")] 
colnames(pathway.df) <- c("Pathway", "Group")
path.bin <- pathway.df %>%
        dplyr::select(Pathway, Group) %>%
        dplyr::mutate(X=1) %>%
        tidyr::spread(Pathway, X) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$Group) %>%
        .[, 2:ncol(.)] %>%
        as.matrix()
      path.bin[is.na(path.bin)] <- 0

      p.mat <- (1 - (path.bin %>% dist(method="binary") %>% as.matrix)) %>% pmin(0.5)
      cl.tree <- dist(p.mat) %>% hclust()
      clust.order <- cl.tree %$% labels[order]
      clusts <- cutree(cl.tree, h=0.7)[clust.order]
      clusts[clusts %in% names(which(table(clusts) < 5))] <- max(clusts) + 1
      clust.lengths <- rle(clusts)$lengths %>% rev
      diag(p.mat) <- 1

      # Plot
      plotHeatmap(p.mat, color.per.group=NULL, row.order=clust.order, col.order=rev(clust.order),
                  legend.title="Similarity", plot.theme=theme_bw()) +
        scale_fill_distiller(palette="RdYlBu", limits=c(0, 0.5))
}
```

``` r
plotCustomSimilarities(gdc)
plotCustomSimilarities(guc)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/hm1.jpg" width="50%" style="display: block; margin: auto;" />

<img src="C:/Users/Katarina/Desktop/scznotebooks/hm2.jpg" width="50%" style="display: block; margin: auto;" />

#### e) Plotting top GO terms for Up and Down DEs

``` r
gcc <- rbindlist(list("down" = gdc, "up" = guc), idcol = "dir")
guc[order(guc$p.adjust2, decreasing = T),]

res.s9.up <- gcc %>% filter(dir == "up")
res.s9.down <- gcc %>% filter(dir == "down")

upgo <- res.s9.up[order(p.adjust2, decreasing = T),]
downgo <- res.s9.down[order(p.adjust2, decreasing = T),]

upgo$Description %<>% make.unique()
upgo$Description %<>% factor(., levels=.)
upgo$p.adjust2 %<>% log10() %>% {. * -1}
downgo$Description %<>% make.unique()
downgo$Description %<>% factor(., levels=.)
downgo$p.adjust2 %<>% log10() %>% {. * -1}
```

``` r
p2 <- ggplot(downgo[1:135,], aes(p.adjust2, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  theme_bw() +
  labs(x="-log10(adj.p)", y="", fill="Cell types", title="Top GO terms for downregulated DE genes") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = high.pal[downgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(downgo[1:135,]$Description), ".", fixed= TRUE), "[[", 1),downgo[1:135,]$Description)) + xlim(0,12.5) 

p1 <- ggplot(downgo[136:270,], aes(p.adjust2, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x="-log10(adj.p)", y="", fill="Cell types", title="All significant GO terms for downregulated DE genes") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = high.pal[downgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(downgo[136:270,]$Description), ".", fixed= TRUE), "[[", 1),downgo[136:270,]$Description)) + xlim(0,12.5) + guides(fill=guide_legend(ncol=1))

legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

x.grob <- textGrob("-log10(adj.p)")
title <- ggdraw() + draw_label("All significant GO terms for downregulated DE genes", fontface='bold')

p <- plot_grid(plotlist = list(p1 + theme(legend.position = "none", axis.title.x = element_blank(), title = element_blank()),p2 + theme(title = element_blank()), legend), rel_widths = c(0.9,0.8, 0.3), ncol = 3, align = "h", axis = "bt")
grid.arrange(arrangeGrob(plot_grid(title, p, ncol=1, rel_heights=c(0.03, 1)) ,bottom = x.grob))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/godown2.jpg" width="50%" style="display: block; margin: auto;" />

``` r
p2 <- ggplot(upgo[1:96,], aes(p.adjust2, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  theme_bw() +
  labs(x="-log10(adj.p)", y="", fill="Cell types", title="Top GO terms for upregulated DE genes") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = high.pal[upgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(upgo[1:96,]$Description), ".", fixed= TRUE), "[[", 1),upgo[1:96,]$Description)) + xlim(0,8) 

p1 <- ggplot(upgo[97:192,], aes(p.adjust2, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x="-log10(adj.p)", y="", fill="Cell types", title="All significant GO terms for upregulated DE genes") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = high.pal[upgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(upgo[97:192,]$Description), ".", fixed= TRUE), "[[", 1),upgo[97:192,]$Description))  + guides(fill=guide_legend(ncol=1))

library(cowplot)
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
x.grob <- textGrob("-log10(adj.p)")
title <- ggdraw() + draw_label("All significant GO terms for upregulated DE genes", fontface='bold')

p <- plot_grid(plotlist = list(p1 + theme(legend.position = "none", axis.title.x = element_blank(), title = element_blank()),p2 + theme(title = element_blank()), legend), rel_widths = c(0.75,0.8, 0.3), ncol = 3, align = "h", axis = "bt")
grid
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/goup2.jpg" width="50%" style="display: block; margin: auto;" />

#### f) Comparison of DEs with Pardinas, SFARI and DISGENET

``` r
sfari <- read.table("/d0-bayes/home/rasmusr/schizo_plots/SFARI_genes.csv", header = F, sep = ",") %>% .[,1]
disgenet <- read.table("/d0-bayes/home/rasmusr/schizo_plots/PsyGeNET_genes.csv", header = F, sep=",") %>% .[,1]
```

``` r
#function of overrepresentation test

overrep.test <- function(n, gene.set) {
  high.top <- de.high %>%
    .[!sapply(., is.logical)] %>%
    lapply(`[[`, 1) %>%
    lapply(function(x) {
      x %>%
        .[complete.cases(.),] %>%
        .[1:n,]
      }) %>%
    .[sapply(., nrow) > 0] %>%
    lapply(rownames)
  
  high.olap <- high.top %>%
    lapply(function(x) sum(x %in% gene.set))
  
  high.p <- lapply(high.olap, function(x) phyper(x-1, length(gene.set), (3e4-length(gene.set)), n, lower.tail=F)) # 3e4 /home/rasmusr number of human genes
  
  message(paste0("Number of significant overrep. cell types: ", sum(high.p <= 0.05, na.rm = T)))
  message(paste0("Significant cell types (ordered): \n"))
  message(paste0(high.p %>% unlist() %>% .[order(.)] %>% .[. <= 0.05] %>% .[!is.na(.)] %>% names(), collapse="\n"))
  
  return(high.p)
}
```

``` r
top500.sfari <- overrep.test(500, sfari)
top500.disgenet <- overrep.test(500, disgenet)
```

### SFARI

``` r
#Figure 6.D
#bonferroni correction
plot.data <- top500.sfari %>% unlist() %>% as.data.frame() %>% setNames("p") %>% mutate(cell=rownames(.)) %>% .[order(.$p, decreasing = T),]
plot.data$cell %<>% factor(., levels=.)
plot.data$p %<>% log10() %>% {. * -1}
plot.sfari <- ggplot(plot.data, aes(p, cell, fill=cell)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  labs(x="-log10(P)", y="", fill="Cell type", title="SFARI") +
  geom_vline(xintercept = -log10(0.05), linetype="dotted") +
  geom_vline(xintercept = -log10(0.05/length(top500.sfari)), linetype="dotted", colour = "red") +
  theme(legend.position = "none") + scale_fill_manual(values = high.pal) +  coord_cartesian(xlim = c(0, 5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/sfari.jpg" width="50%" style="display: block; margin: auto;" />

### DISGENET

``` r
#Figure 6.C
#bonferroni correction
#for bonferroni we plot p values and cutoffs for p value and p adjusted to show how many cell types with significant p values had also significant p adjusted values
plot.data <- top500.disgenet %>% unlist() %>% as.data.frame() %>% setNames("p") %>% mutate(cell=rownames(.)) %>% .[order(.$p, decreasing = T),]
plot.data$cell %<>% factor(., levels=.)
plot.data$p %<>%log10() %>% {. * -1}
plot.disgenet <- ggplot(plot.data, aes(p, cell, fill=cell)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  labs(x="-log10(P)", y="", fill="Cell type", title="DisGeNET") +
  geom_vline(xintercept = -log10(0.05), linetype="dotted") +
  geom_vline(xintercept = -log10(0.05/length(top500.psygenet)), linetype="dotted", colour = "red") +
  theme(legend.position = "none") + scale_fill_manual(values = high.pal) +  coord_cartesian(xlim = c(0, 5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/disgenet.jpg" width="50%" style="display: block; margin: auto;" />

### Pardinas

  - For this dataset we used a special method called LDSC (please refer
    to <https://github.com/bulik/ldsc>)
  - the method is implemented within CELLECT:
    <https://github.com/perslab/CELLECT>
  - for CELLECT, the result from CELLEX is used for input:
    <https://github.com/perslab/CELLEX>

#### Preparing data for CELLEX:

``` r
cm.merged.raw <- con$getJointCountMatrix(raw=T)
```

Gene SYMBOLS -\> ENSEMBL ID

``` r
gnames <- clusterProfiler::bitr(geneID = colnames(cm.merged.raw), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
gnames %<>% .[match(unique(.$SYMBOL), .$SYMBOL),] %>% .[match(unique(.$ENSEMBL), .$ENSEMBL),]

cm.merged.raw %<>% .[,match(gnames$SYMBOL, colnames(.))]
colnames(cm.merged.raw) <- gnames$ENSEMBL
```

``` r
cell.groups <- tans$high.clean %>%
  as.character() %>%
  setNames(names(tans$high.clean)) %>%
  .[match(rownames(cm.merged.raw),names(.))] %>%
  cbind() %>% 
  data.frame(stringsAsFactors = F) %>% 
  tibble::rownames_to_column() %>% 
  setNames(c("cell_id","cell_type"))
cell.groups$cell_id %<>% gsub(".", "-", ., fixed = T)
```

Isolate SZ cells

``` r
sample.per.cell <- con$getDatasetPerCell()
scz.cells <- sample.per.cell[sample.per.cell %in% names(samplegroups[samplegroups == "Scz"])] %>% factor()
groups.scz <- cell.groups %>% .[.$cell_id %in% names(scz.cells),]
scz.cm <- cm.merged.raw %>% .[rownames(.) %in% groups.scz$cell_id,]
dim(scz.cm); dim(groups.scz) #check if dimensions of matrix and cells are fine
```

#### Save data for Python

``` r
library(HDF5Array)
library(Matrix)

writeHDF5Array(scz.cm, "/d0/home/kdragicevic/CELLEX/scz2.h5", "data", verbose = T)
write.csv(scz.cm %>% colnames(), "/d0/home/kdragicevic/CELLEX/genes2.scz.csv", row.names = F)
write.csv(scz.cm %>% rownames(), "/d0/home/kdragicevic/CELLEX/cells2.scz.csv", row.names = F)
write.csv(groups.scz, "/d0/home/kdragicevic/CELLEX/metadata.scz.csv", row.names = F)
```

### PYTHON:

``` r bg-success
#note: run in SCREEN session because it might take some time to import huge h5 objects - prevents the session from being interrupted

cd /home/kdragicevic/CELLEX
python


import numpy as np
import pandas as pd
import cellex
import h5py

data = pd.DataFrame(np.array(h5py.File("scz2.h5")['data']), dtype="float32")
data.columns = pd.read_csv("cells2.scz.csv").values
data.columns = [v[0] for v in data.columns]
data.index = pd.read_csv("genes2.scz.csv").values
data.index = [v[0] for v in data.index]

metadata = pd.read_csv("metadata.scz.csv", index_col=0)

eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)
eso.results["esmu"].to_csv("Cellex.scz.results.csv.gz")
```

Preparing data for CELLECT

SZ data <https://atlas.ctglab.nl/traitDB/3982>

Munging

``` r bg-success
cd /home/rasmusr/CELLECT
conda activate munge_ldsc

python ldsc/mtag_munge.py \
--sumstats Pardinas2018/Pardinas2018_clean.tsv \
--a1 A1 \
--a2 A2 \
--merge-alleles data/ldsc/w_hm3.snplist \
--keep-pval \
--p PVAL \
--N-cas 40675 \
--N-con 64643 \
--out Pardinas2018/Pardinas2018
```

Run CELLECT LDSC

``` r bg-success
conda activate base
snakemake --use-conda -j 20 -s cellect-ldsc.snakefile --configfile example/config_SCZ.yml
```

``` r
dat.ld <- read.csv("/d0/home/rasmusr/CELLECT/CELLECT-LDSC/prioritization.csv")

dat.ld %<>% dplyr::mutate(type="LD") %>% .[order(.$annotation),]
p.cutoff <- 0.05/nrow(dat.ld)
cellect <- ggplot(dat.ld, aes(annotation, pvalue, col=annotation)) + 
  geom_point() +
  labs(x="", y="P value") +
  geom_hline(yintercept = 0.05, col = "black", linetype = "dotted") +
  geom_hline(yintercept = p.cutoff, col = "red", linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none", axis.text.y = element_text(colour = ifelse(dat.ld$pvalue <= p.cutoff, "red", "black"))) +
  coord_flip() + 
  scale_color_manual(values = dat.ld$pvalue %>% sapply(function(x) {
    if(x <= p.cutoff) return("red") else if(x > p.cutoff & x <= 0.05) return("gray") else return("black")
  }) %>% setNames(dat.ld$annotation)) 
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/pardinas.jpg" width="50%" style="display: block; margin: auto;" />
