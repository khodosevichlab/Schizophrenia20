#' @import dplyr

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

calculateGos <- function(de, go.datas, n.top.genes=300,n.cores=1) {
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
              all=list(all=unique(unlist(lapply(de,function(x) rownames(x$res))))))

  gns.entrez <- lapply(gns,function(x)
    lapply(x, clusterProfiler::bitr, 'SYMBOL', 'ENTREZID', org.Hs.eg.db::org.Hs.eg.db) %>%
      lapply(`[[`, "ENTREZID"))

  gos <- lapply(gns.entrez[c('up','down')],function(gns) {
    lapply(gns, enrichGOOpt, universe=gns.entrez$all$all, ont='BP', goData=go.datas[['BP']],
           readable=T, OrgDB=org.Hs.eg.db::org.Hs.eg.db)# %>% lapply(function(x) x@result)
  })

  return(gos)
}

clusterGos <- function(gos,n.clusters=20,max.pval=0.05) {
  gos_filt <- lapply(gos,function(x) filter(x@result,p.adjust<max.pval))
  gos_joint <- do.call(rbind,gos_filt)

  gos_joint <- gos_filt %>% .[sapply(., nrow) > 0] %>% names() %>% setNames(., .) %>% lapply(function(n) cbind(gos_filt[[n]],Type=n)) %>% Reduce(rbind,.)

  go_bp_df <- gos_joint %>% group_by(Type, Description) %>%
    summarise(p.adjust=min(p.adjust)) %>% ungroup() %>% mutate(p.adjust=-log10(p.adjust)) %>%
    tidyr::spread(Type, p.adjust) %>% as.data.frame() %>% set_rownames(.$Description) %>% .[, 2:ncol(.)] #%>% .[, type_order[type_order %in% colnames(.)]]
  go_bp_df[is.na(go_bp_df)] <- 0
  go_dist <- distanceBetweenTerms(gos_joint)

  if (n.clusters <= ncol(as.matrix(go_dist))) {
    clusts <- hclust(go_dist,method='ward.D2') %>% cutree(min(n.clusters, ncol(go_dist)))
    gos_per_clust <- split(names(clusts), clusts)
    ngos_per_clust <- sapply(gos_per_clust, length)
    #table(clusts)

    gos_per_clust <- split(names(clusts), clusts)
    gos_joint %<>% mutate(GOClust=clusts[Description])
    name_per_clust <- gos_joint %>% group_by(GOClust, Description) %>% summarise(pvalue=exp(mean(log(pvalue)))) %>%
      split(.$GOClust) %>% sapply(function(df) df$Description[which.min(df$pvalue)])
    gos_joint %<>% mutate(GOClustName=name_per_clust[as.character(GOClust)])

    # cluster summary
    go_bp_summ_df <- gos_joint %>% group_by(Type, GOClustName) %>%
      summarise(p.adjust=min(p.adjust)) %>% ungroup() %>% mutate(p.adjust=-log10(p.adjust)) %>%
      tidyr::spread(Type, p.adjust) %>% as.data.frame() %>% set_rownames(.$GOClustName) %>% .[, 2:ncol(.)] #%>% .[, type_order[type_order %in% colnames(.)]]
    go_bp_summ_df[is.na(go_bp_summ_df)] <- 0
  } else {
    clusts <- gos_joint$Description %>% setNames(., .)
    go_bp_summ_df <- go_bp_df
  }

  return(list(joint=gos_joint, clusters=clusts, summary=go_bp_summ_df, go_df=go_bp_df))
}
