Notebook 7: MTG Allen correspondence
================

``` r
library(pagoda2)
  library(conos) #!!!USED 1.2.1 version of conos here
  library(magrittr)
  library(scrattch.io)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(pals)
  library(uwot)
  library(forcats)
  library(ggforce)
  library(Hmisc)
  library(viridis)
  library(rhdf5) #from bioconductor, for hdf5 matrices
  library(readxl) #import xlsx files
  library(dplyr)
  library(writexl)
  library(cowplot) #plot_grid function
  library(patchwork)
  library(ggalluvial)
  library("ggrepel")
  library("car") #leveneTest
  library("ggsignif") #geom_signif stat_signif
```

### 1\. Create conos object of list of count matrices that contain our data + TC/MTG data

``` r
#pagoda2
p2.list.mtg.schizo <- lapply(cms_mtg_schizo, basicP2proc, n.cores = 30, n.odgenes = 3000, min.cells.per.gene = 0,
                  min.transcripts.per.cell = 0, get.largevis = F, get.tsne = T, make.geneknn = F
)
#conos
con_mtg_schizo <- Conos$new(p2.list.mtg.schizo,n.cores = 30)
con_mtg_schizo$buildGraph(k = 30, k.self = 5, 
               space = "PCA",
               ncomps = 30, n.odgenes = 3000, matching.method = "mNN", metric = "angular",
               score.component.variance = T, verbose = T, 
               
               
               balancing.factor.per.cell = setNames(
                        con_mtg_schizo$getDatasetPerCell() %>% gsub("MB.*", "SCHIZO", .) %>% gsub("^H200.*", "ALLEN", .),
                        names(con_mtg_schizo$getDatasetPerCell())
                        ),
               same.factor.downweight = 0.1,               
               
               alignment.strength = 0.2
)
#create UMAP
con_mtg_schizo$embedGraph(method = 'UMAP', min.dist = 0.1, spread = 10, n.cores = 30, min.prob.lower = 1e-7)
```

### 2\. propagate conos labels from our scz data to allen data

``` r
#annotations_scz2 is the vector containing annotation for our scz data - cluster labels of cells
#propage label, high resolution
annot_tranfer_mtg_highres <-
  con_mtg_schizo$propagateLabels(labels = annotations_scz2$high, verbose = T)

#propage label, mid resolution
annot_tranfer_mtg_midres <-
  con_mtg_schizo$propagateLabels(labels = annotations_scz2$med, verbose = T)
```

### 3\. clean data

``` r
#REMOVE LABELS OF SEVERAL OPCS AND MICROGLIA CELLS THAT APPEAR AS OTHER CELL TYPES (VIP)
#REMOVE ABBERANTLY ASSIGNED CELLS

#make vector with names of cells to be removed
mtg_remove_cells_hres <-  annot_tranfer_mtg_highres$labels[con_allen_mtg$clusters$leiden$groups[con_allen_mtg$clusters$leiden$groups == "16" | con_allen_mtg$clusters$leiden$groups == "17"] %>% 
                                                             names] %>%
    grepl("Glia|Other", .) %>% 
  {!.} %>% 
  (con_allen_mtg$clusters$leiden$groups[con_allen_mtg$clusters$leiden$groups == "16" | con_allen_mtg$clusters$leiden$groups == "17"] %>%
     names)[.]




mtg_remove_cells_midres <- annot_tranfer_mtg_midres$labels[con_allen_mtg$clusters$leiden$groups[con_allen_mtg$clusters$leiden$groups == "16" | con_allen_mtg$clusters$leiden$groups == "17"] %>% 
                                                             names] %>%grepl("Glia|Other", .) %>% 
  {!.} %>% 
  (con_allen_mtg$clusters$leiden$groups[con_allen_mtg$clusters$leiden$groups == "16" | con_allen_mtg$clusters$leiden$groups == "17"] %>%
     names)[.]


#remove cells from hres annotation transfer
annot_tranfer_mtg_highres$labels <- 
  annot_tranfer_mtg_highres$labels[!names(annot_tranfer_mtg_highres$labels) %in% mtg_remove_cells_hres]

annot_tranfer_mtg_highres$label.distribution <- 
  annot_tranfer_mtg_highres$label.distribution[!rownames(annot_tranfer_mtg_highres$label.distribution) %in% mtg_remove_cells_hres, ]

annot_tranfer_mtg_highres$uncertainty <- 
  annot_tranfer_mtg_highres$uncertainty[! names(annot_tranfer_mtg_highres$uncertainty) %in% mtg_remove_cells_hres]



#remove cells from midres annotation transfer
annot_tranfer_mtg_midres$labels <- 
  annot_tranfer_mtg_midres$labels[! names(annot_tranfer_mtg_midres$labels) %in% mtg_remove_cells_midres]

annot_tranfer_mtg_midres$label.distribution <- 
  annot_tranfer_mtg_midres$label.distribution[! rownames(annot_tranfer_mtg_midres$label.distribution) %in% mtg_remove_cells_midres, ]

annot_tranfer_mtg_midres$uncertainty <- 
  annot_tranfer_mtg_midres$uncertainty[! names(annot_tranfer_mtg_midres$uncertainty) %in% mtg_remove_cells_midres]
```

### 4\. clean labels in annot transfer

``` r
annot_tranfer_mtg_midres$labels %<>% 
  gsub("L2_CUX2_FREM3", "L2_3_CUX2_FREM3", .) %>% 
  gsub("L2_CUX2_PRSS12", "L3_CUX2_PRSS12", .)

colnames(annot_tranfer_mtg_midres$label.distribution) %<>% 
  gsub("L2_CUX2_FREM3", "L2_3_CUX2_FREM3", .) %>% 
  gsub("L2_CUX2_PRSS12", "L3_CUX2_PRSS12", .)
```

### 5\. Subtype probability plot, MTG/TC, high resolution

``` r
#Plot probability HIGH RES of each nuclei to belong to a specific cluster
#select Allen cells
probab.mtg.hres.df <- annot_tranfer_mtg_highres$label.distribution %>% rownames() %>% 
  grepl("^MB", .) %>% 
  {!.} %>% 
  annot_tranfer_mtg_highres$label.distribution[.,]  %>% 
  as.data.frame()


#add cluster names
probab.mtg.hres.df$cluster <- annot_tranfer_mtg_highres$labels %>% names() %>% 
  grepl("^MB", .) %>% 
  {!.} %>% 
  annot_tranfer_mtg_highres$labels[.] %>% 
  factor(levels = c("L2_CUX2_LAMP5_MARCH1", "L2_CUX2_LAMP5_PDGFD", "L2_3_CUX2_FREM3_UNC5D", "L2_3_CUX2_FREM3_SV2C", "L3_CUX2_PRSS12", "L4_RORB_SCHLAP1_MET", "L4_RORB_SCHLAP1_MME", "L4_RORB_SCHLAP1_ARHGAP15", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A", "L5_6_FEZF2_TLE4_ABO", "L5_6_FEZF2_TLE4_SCUBE1","L5_6_FEZF2_TLE4_HTR2C", "L5_6_THEMIS_SEMA3A", "L5_6_THEMIS_NTNG2", "ID2_LAMP5_NMBR", "ID2_LAMP5_CRH", "ID2_LAMP5_NOS1", "ID2_PAX6", "ID2_NCKAP5", "VIP_ABI3BP", "VIP_TYR", "VIP_RELN", "VIP_SEMA3C", "VIP_SSTR1", "VIP_CRH", "PVALB_MEPE", "PVALB_SST", "PVALB_CRH", "SST_TH", "SST_TAC3", "SST_CALB1","SST_NPY", "SST_NOS1", "SST_STK32A", "Glia", "Other"))
  
  
#get probability of nucleus belonging to a certain cluster
probab.mtg.hres.df$cluster_probability <-  probab.mtg.hres.df %>%
  apply(., 1, function(x) {x[colnames(.) == x["cluster"]]}) %>% 
  as.numeric()

#plot probability of nuclei belonging to clusters
plot.probab.mtg.hres <-
  ggplot(probab.mtg.hres.df, aes(x = cluster, y = cluster_probability, color = cluster)) + 
    geom_sina(show.legend = F, alpha = 0.2, size = 1, scale = "width", maxwidth = 0.6) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", size = 0.2,
                   show.legend = FALSE, color = "black", shape = 23, fill = "grey") +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_color_manual(values = high.pal[levels(annotations_scz2$high)]) +
  labs(title = "High resolution subtype probability, MTG estimate",  y = "Probability") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/tch1.jpg" width="60%" style="display: block; margin: auto;" />

### 6\. Subtype probability plot, MTG/TC, medium resolution

``` r
#Plot probability MIDDLE RES of each nuclei to belong to a specific cluster
#select Allen cells
probab.mtg.mres.df <-
  annot_tranfer_mtg_midres$label.distribution %>% 
  rownames() %>% 
  grepl("^MB", .) %>% 
  {!.} %>% 
  annot_tranfer_mtg_midres$label.distribution[.,]  %>% 
  as.data.frame()


#add cluster names
probab.mtg.mres.df$cluster <-
  annot_tranfer_mtg_midres$labels %>% 
  names() %>% 
  grepl("^MB", .) %>% 
  {!.} %>%
  annot_tranfer_mtg_midres$labels[.] %>% 
  factor(levels = c("L2_CUX2_LAMP5", "L2_3_CUX2_FREM3", "L3_CUX2_PRSS12", "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A", "L5_6_FEZF2_TLE4", "L5_6_THEMIS","ID2_LAMP5", "ID2_PAX6", "ID2_NCKAP5", "VIP", "PVALB", "SST", "Glia", "Other"))
                                    
  
  
#get probability of nucleus belonging to a certain cluster
probab.mtg.mres.df$cluster_probability <- probab.mtg.mres.df %>%
  apply(., 1, function(x) {x[colnames(.) == x["cluster"]]}) %>%
  as.numeric()

#plot probability of nuclei belonging to clusters
plot.probab.mtg.mres <-
  ggplot(probab.mtg.mres.df, aes(x = cluster, y = cluster_probability, color = cluster)) + 
    geom_sina(show.legend = F, alpha = 0.2, size = 1, scale = "width", maxwidth = 0.6) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", size = 0.2,
                   show.legend = FALSE, color = "black", shape = 23, fill = "grey") +
  coord_cartesian(ylim = c(0, 1.05)) +
 scale_color_manual(values = med.pal[levels(annotations_scz2$med)]) +
  labs(title = "Medium resolution subtype probability, MTG estimate",  y = "Probability") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/tcm1.jpg" width="40%" style="display: block; margin: auto;" />

### 7\. Position plot, MTG/TC estimate, high resolution

``` r
layers.hres.mtg <-
  ggplot(data.frame(
    Subtypes = annot_tranfer_mtg_highres$labels[annotation_allen$sample_name[annotation_allen$sample_name %in% names(annot_tranfer_mtg_highres$labels)]] %>% 
         factor(levels = c("L2_CUX2_LAMP5_MARCH1", "L2_CUX2_LAMP5_PDGFD", "L2_3_CUX2_FREM3_UNC5D", "L2_3_CUX2_FREM3_SV2C", "L3_CUX2_PRSS12", "L4_RORB_SCHLAP1_MET", "L4_RORB_SCHLAP1_MME", "L4_RORB_SCHLAP1_ARHGAP15", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A", "L5_6_FEZF2_TLE4_ABO", "L5_6_FEZF2_TLE4_SCUBE1","L5_6_FEZF2_TLE4_HTR2C", "L5_6_THEMIS_SEMA3A", "L5_6_THEMIS_NTNG2", "ID2_LAMP5_NMBR", "ID2_LAMP5_CRH", "ID2_LAMP5_NOS1", "ID2_PAX6", "ID2_NCKAP5", "VIP_ABI3BP", "VIP_TYR", "VIP_RELN", "VIP_SEMA3C", "VIP_SSTR1", "VIP_CRH", "PVALB_MEPE", "PVALB_SST", "PVALB_CRH", "SST_TH", "SST_TAC3", "SST_CALB1","SST_NPY", "SST_NOS1", "SST_STK32A", "Glia", "Other")),
    Layers = annotation_allen$cortical_layer_label[annotation_allen$sample_name %in% names(annot_tranfer_mtg_highres$labels)] %>%
      factor(levels = c("L6", "L5", "L4", "L3", "L2", "L1"))), 
    aes(x = Subtypes, y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.7, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = high.pal[levels(annotations_scz2$high)]) +
  labs(title = "Layer distribution of high resolution subtypes, MTG estimate") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        text = element_text(family = "Helvetica"), axis.title.x = element_blank())
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/tch2.jpg" width="60%" style="display: block; margin: auto;" />

### 8\. Position plot, MTG/TC estimate, medium resolution

``` r
layers.mres.mtg <-
 ggplot(data.frame(
   Subtypes = annot_tranfer_mtg_midres$labels[annotation_allen$sample_name[annotation_allen$sample_name %in% names(annot_tranfer_mtg_midres$labels)]] %>% 
     factor(levels = c("L2_CUX2_LAMP5", "L2_3_CUX2_FREM3", "L3_CUX2_PRSS12", "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A", "L5_6_FEZF2_TLE4", "L5_6_THEMIS","ID2_LAMP5", "ID2_PAX6", "ID2_NCKAP5", "VIP", "PVALB", "SST", "Glia", "Other")),
   Layers = annotation_allen$cortical_layer_label[annotation_allen$sample_name %in% names(annot_tranfer_mtg_midres$labels)] %>%
     factor(levels = c("L6", "L5", "L4", "L3", "L2", "L1"))), 
   aes(x = Subtypes, y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.7, show.legend = F) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        text = element_text(family = "Helvetica"), axis.title.x = element_blank()) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = med.pal[levels(annotations_scz2$med)]) +
  labs(title = "Layer distribution of medium resolution subtypes, MTG estimate")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/tcm2.jpg" width="40%" style="display: block; margin: auto;" />

``` r
annotation_allen <- fread(file = "~/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/sample_annotations.csv")
```

``` r
#High Res
df.correspond.hres.mtg <- data.frame(KU = factor(annot_tranfer_mtg_highres$labels[!names(annot_tranfer_mtg_highres$labels) %>% 
                                                            grepl("MB", .)],levels = hsub_order),
            Allen = setNames(annotation_allen$cluster_label,annotation_allen$sample_name)[names(annot_tranfer_mtg_highres$labels[!names(annot_tranfer_mtg_highres$labels) %>% 
                                                            grepl("MB", .)])])


#sumamrize data and filter - remove correspondence <=7% of corresponding Allen subtypes - some kind of trash
df.correspond.hres.mtg <- df.correspond.hres.mtg %>% 
  group_by(., KU, Allen) %>% 
  summarize(., n =  n()) %>% 
  mutate(freq = n / sum(n)) %>%
  ungroup() %>% 
  filter(freq > 0.07)

df.correspond.hres.mtg$Allen <- factor(df.correspond.hres.mtg$Allen, levels = df.correspond.hres.mtg$Allen %>% unique)
hsub_order_allen_mtg <- df.correspond.hres.mtg %>% 
  split(.$Allen) %>% 
  lapply(function(x) {x %$% setNames(n, KU) %>% 
      `/`(sum(.))}) %>% 
  sapply(function(y) as.numeric(hsub_ranks[names(y)[which.max(y)]]) +  sum(hsub_ranks[names(y)] * y) / 10) %>% 
  sort() %>% 
  names()

df.correspond.hres.mtg %<>% mutate(KU = factor(KU, levels = hsub_order), Allen = factor(Allen, levels = hsub_order_allen_mtg))
```

``` r
#Medium Res
df.correspond.mres.mtg <-  data.frame(KU = factor(annot_tranfer_mtg_midres$labels[!names(annot_tranfer_mtg_midres$labels) %>% 
                                                                                    grepl("MB", .)],levels = msub_order),
       Allen = setNames(annotation_allen$cluster_label, annotation_allen$sample_name)[names(annot_tranfer_mtg_midres$labels[!names(annot_tranfer_mtg_midres$labels) %>% 
                                                                                    grepl("MB", .)])])


#sumamrize data and filter - remove correspondence <=7% of corresponding Allen subtypes - some kind of trash
df.correspond.mres.mtg <- df.correspond.mres.mtg %>% group_by(., KU, Allen) %>% summarize(., n =  n()) %>% mutate(freq = n / sum(n)) %>%
  ungroup() %>% filter(freq > 0.07)

df.correspond.mres.mtg$Allen <- factor(df.correspond.mres.mtg$Allen, levels = df.correspond.mres.mtg$Allen %>% unique)


msub_order_allen_mtg <- df.correspond.mres.mtg %>% 
  split(.$Allen) %>% 
  lapply(function(x) {x %$% setNames(n, KU) %>% 
      `/`(sum(.))}) %>% 
  sapply(function(y) as.numeric(msub_ranks[names(y)[which.max(y)]])  +  sum(msub_ranks[names(y)] * y) / 10) %>% 
  sort() %>% 
  names()

df.correspond.mres.mtg %<>% mutate(KU = factor(KU, levels = msub_order), Allen = factor(Allen, levels = msub_order_allen_mtg))
```

Plot them

``` r
plot.corresp.hres.mtg <- plotAlluviumHigh(df.correspond.hres.mtg, 0.65, 0.58)
plot.corresp.mres.mtg <- plotAlluviumMed(df.correspond.mres.mtg, 0.4, 0.53)
```

``` r
fig.s6_correspond <-
  
  plot.corresp.hres.mtg + plot.corresp.mres.mtg +
  plot_layout(ncol = 2) &
  plot_annotation(title = "KU - Allen MTG correnspondence") & theme(plot.title = element_text(family = "Liberation Sans", size = 18, face = "bold",
                                                                                                             hjust = 0.5))

fig.s6_correspond
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/kuallen.jpg" width="60%" style="display: block; margin: auto;" />
