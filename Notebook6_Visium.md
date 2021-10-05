Notebook 5: Visium data
================

### 1\. Integrate data: MB14, MB15, MB17 + allen(glia)

``` r
cms_scz <- lapply(setNames(c("MB14", "MB15", "MB17"), c("MB14", "MB15", "MB17")), function(smplname){
             cm <- Seurat::Read10X_h5(
                paste0("/home/mbatiuk/projects/human_schizo/MB6-MB23/counts_Allen_Genome/MB14_MB15_MB17_proper_gtf/",
                       smplname, "/outs/filtered_feature_bc_matrix.h5"))
             colnames(cm) <- paste0(smplname, "_", colnames(cm))
             cm})
```

``` r
annotations_scz2 <- read_rds("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/annotations_scz2.rds")
cms_scz <- lapply(cms_scz, function(smpl) {smpl[, colnames(smpl) %in% names(annotations_scz2$high.clean)]}) 
```

Filter Allen dataset, retain non-neural cells

``` r
#load Allen annotation

annotation_allen <- fread("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/sample_annotations.csv")

#schizo annotations and subtype orders
annotations_scz2 <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/annotations_scz2.rds")

hsub_order <- read_rds("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/hsub_order.rds")
msub_order <- read_rds("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/msub_order.rds")


#load Allen annotation
annotation_allen <- fread(file = "/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/sample_annotations.csv")

#Allen 2019 human data, smart-seq4 
#https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/transcriptome.zip


allen_exon <- read_tome_dgCMatrix(tome = "/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/transcrip.tome", 
                                  target = "data/t_exon")
allen_intron <- read_tome_dgCMatrix(tome = "/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/transcrip.tome", 
                                    target = "data/t_intron")


#load sample and gene names
allen_smpl_names <- read_tome_sample_names(tome = "/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/transcrip.tome")
allen_gene_names <- read_tome_gene_names(tome = "/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/transcrip.tome")

#sum up introns + exons:
allen_cm <- allen_exon + allen_intron

#add colnames and rownames to matrix
rownames(x = allen_cm) <- allen_gene_names
colnames(x = allen_cm) <- allen_smpl_names
```

``` r
allen_list_1 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  annotation_allen$region_label == "CgG" &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("ACC_", annotation_allen$external_donor_name_label %>% unique()))


allen_list_2 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  annotation_allen$region_label == "MTG"  &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("MTG_", annotation_allen$external_donor_name_label %>% unique()))


allen_list_3 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  (annotation_allen$region_label == "M1ul" | annotation_allen$region_label == "M1lm") &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("M1_", annotation_allen$external_donor_name_label %>% unique()))


allen_list_4 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  annotation_allen$region_label == "V1C"  &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("V1C_", annotation_allen$external_donor_name_label %>% unique()))




allen_list_5 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  (annotation_allen$region_label == "S1ul" | annotation_allen$region_label == "S1lm") &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("S1_", annotation_allen$external_donor_name_label %>% unique()))


allen_list_6 <- lapply(annotation_allen$external_donor_name_label %>% unique(), 
                     function(x) 
                       allen_cm[ ,annotation_allen$external_donor_name_label == x &
                                  annotation_allen$region_label == "A1C"  &
                                   annotation_allen$class_label == "Non-neuronal"]) %>% 
  setNames(paste0("A1C_", annotation_allen$external_donor_name_label %>% unique()))

allen_list <- c(allen_list_1, allen_list_2, allen_list_3, allen_list_4, allen_list_5, allen_list_6)
```

Make Seurat Visium object

``` r
#merge allen and scz countmatrics lists
scz_all_cms <- append(cms_scz, allen_list)


#create Seurat objects
seu_scz_all <- lapply(scz_all_cms, CreateSeuratObject)


#run SCTransform on all memebers
seu_scz_all <- lapply(seu_scz_all, SCTransform) 


#select features for downstream integration
scz_all_features <- SelectIntegrationFeatures(object.list = seu_scz_all, nfeatures = 3000)

#run PrepSCTIntegration
seu_scz_all <- PrepSCTIntegration(object.list = seu_scz_all, anchor.features = scz_all_features)


#find anchors
scz_all_anchors <- FindIntegrationAnchors(object.list = seu_scz_all, normalization.method = "SCT", 
                                          anchor.features = scz_all_features,
                                          k.filter = 130 #needed to decrease, error when integr small and big datasets
                                          )


#integrate schizo and allen data
scz_all_anchors <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/rds_objects/scz_all_anchors.rds")
scz_all_integrated <- IntegrateData(anchorset = scz_all_anchors, normalization.method = "SCT")


#Run PCA
scz_all_integrated <- RunPCA(scz_all_integrated)

#Run UMAP
scz_all_integrated <- RunUMAP(scz_all_integrated, dims = 1:30)
saveRDS(scz_all_integrated, "scz_all_integrated.RDS")
```

#### load data

``` r
cm_glia <- read_rds("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/rds_objects/allen_10x_m1_non_neur.rds")
annot_glia <- read_rds("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/rds_objects/anno_allen_10x.rds")
so <- readRDS("scz_all_integrated.RDS")
annotation <- read_rds("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/rds_objects/annotation_scz_allen.rds")

so$subtypes_med <- annotation$med[names(Idents(so))]
so$subtypes_high <- annotation$high[names(Idents(so))]
so$origin <- setNames(
  gsub("^[^(MB)].*", "Allen", names(annotation$high)) %>% gsub("^MB.*", "Scz", .),
  names(annotation$high))[names(Idents(so))]

so_subs <- subset(so, (origin == "Scz") & (subtypes_high != "Glia"))

#embeddingPlot(so_subs@reductions$umap@cell.embeddings, groups=so$subtypes_high, size=0.1, font.size=c(2,3))
```

``` r
visium_allen_old <- read_rds("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/rds_objects/visium_allen.rds")
```

Estimate variable genes for stereoscope deconvolotion

``` r
var_genes_3k <- visium_allen_old$MB11 %>% FindVariableFeatures(assay="Spatial", nfeatures=3000) %>% 
  .@assays %>% .$Spatial %>% .@var.features
var_genes_1k <- visium_allen_old$MB11 %>% FindVariableFeatures(assay="Spatial", nfeatures=1000) %>% 
  .@assays %>% .$Spatial %>% .@var.features
```

#### Save the data for stereoscope

``` r
asChn <- function(x) setNames(as.character(x), names(x))
```

``` r
cm_merged <- so_subs@assays$RNA@counts %>% list(cm_glia) %>% 
  conos:::mergeCountMatrices(transposed=F)
annotation_merged <- asChn(so_subs$subtypes_med) %>% c(asChn(annot_glia$med))
annotation_merged_high <- asChn(so_subs$subtypes_high) %>% c(asChn(annot_glia$high))
p_vals <- colSums(cm_merged > 0) %>% split(annotation_merged[names(.)]) %>% 
  sapply(mean)
ggplot(tibble(y=p_vals, n=names(p_vals))) +
  geom_bar(aes(x=n, y=y), stat="identity") +
  theme(axis.text=element_text(angle=45, hjust=1))
```

``` r
cm_merged[var_genes_3k,] %>% Matrix::t() %>% as.matrix() %>% 
  as.data.frame() %>% data.table::fwrite("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/cm_with_glia2.tsv", sep="\t", row.names=T)
cm_norm <- Pagoda2$new(cm_merged)$counts
cm_norm[, var_genes_3k] %>% as.matrix() %>% as.data.frame() %>% 
  data.table::fwrite("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/cm_with_glia_norm2.tsv", sep="\t", row.names=T)
```

``` r
tr <- names(visium_allen) %>% 
    pblapply(function(n) {
      cm <- visium_allen[[n]]@assays$Spatial@counts
      Matrix::t(cm[var_genes_3k,]) %>% as.matrix() %>% as.data.frame() %>% 
        data.table::fwrite(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/", n, "_cm.tsv"), sep="\t", row.names=T)
      cm.norm <- Pagoda2$new(cm, verbose=FALSE)$counts
      cm.norm[, var_genes_3k] %>% as.matrix() %>% as.data.frame() %>% 
        data.table::fwrite(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/", n, "_cm_norm.tsv"), sep="\t", row.names=T)
      
      cm.norm[, var_genes_1k] %>% as.matrix() %>% as.data.frame() %>% 
        data.table::fwrite(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/", n, "_cm_norm_1k.tsv"), sep="\t", row.names=T)
})
```

``` r
as_tibble(annotation_merged, rownames="cell") %>% 
  set_colnames(c("cell", "bio_celltype")) %>% 
  write_delim("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/annotation_med2.tsv", delim="\t")
as_tibble(annotation_merged_high, rownames="cell") %>% 
  set_colnames(c("cell", "bio_celltype")) %>% 
  write_delim("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/annotation_med_high.tsv", delim="\t")
```

STEREOSCOPE

``` r bg-success
taskset --cpu-list 0-14 stereoscope run --sc_cnt cm_with_glia2.tsv --sc_labels annotation_med_high.tsv -sce 20000  -o deconv -n 5000 --st_cnt MB*cm.tsv -ste 30000 -stb 1000 -scb 1000
```

### 2\. Read results from stereoscope and plot visium panel on MB11

``` r
library(pagoda2)
library(magrittr)
library(tidyverse)
library(pbapply)
library(dataorganizer)
library(Seurat)
library(scrattch.io)
```

Code below is from revision, therefore some chunks might contatin just
functions which load new data (since cohort 1 was created in the exactly
same way)

``` r
#find the ensemble count folder of spatial data cohort 2
dataPath <- function(...) 
  file.path("/d0/home/mbatiuk/projects/human_schizo/primary/spatial/counts/ENSEMBLE_97_genome_cohort2/", ...)
annotationPath <- function(...)
  file.path("/d0/home/kdragicevic/RProjects/schizo/spatial/layeranno", ...)
kDatasetNames <- c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52") %>% 
  setNames(., .)
#cohort 2:
kConditionPerDataset <- c(MB6="SCZ", MB50="SCZ", MB52="SCZ", 
                          MB9="CNT", MB19="CNT", MB21="CNT")
```

``` r
samples <- kDatasetNames %>% pblapply(function(n) {
 sp <- dataPath(n, "outs/spatial/tissue_positions_list.csv") %>% 
     read_csv(col_names=F) %>% as.data.frame() %>% set_rownames(.$X1) %>% .[,5:6] %>% 
     set_colnames(c("x", "y"))
   cm <- dataPath(n, "outs/filtered_feature_bc_matrix.h5") %>% Seurat::Read10X_h5()
   ann <- paste0(n, "_MB.csv") %>% annotationPath() %>% read_csv() %$% 
     setNames(.[[2]], Barcode) %>% .[. != "WM"]
   cm <- cm[, names(ann)]
   sp <- sp[names(ann),]
   colnames(cm) <- rownames(sp) <- names(ann) <- paste0(n, "_", names(ann))
   p2 <- vpscutils::GetPagoda(cm, n.cores=10, embeding.type="UMAP", verbose=FALSE)
   list(cm=cm, annotation=ann, position=sp, p2=p2)})
 
 write_rds(samples, "/d0/home/kdragicevic/RProjects/schizo/spatial/spatial_samples2.rds")
```

#### samples1 (cohort 1) and samples 2 (cohort 2)

``` r
#samples_sp1 <-  read_rds("/d0-mendel/home/viktor_petukhov/Copenhagen/Schizophrenia20/cache/spatial_samples.rds")
#samples_sp2 <- read_rds("/d0/home/kdragicevic/RProjects/schizo/spatial/spatial_samples2.rds")
```

#### subtype position in tissue plot

``` r
visium_allen_plot <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")),
                       function(smpl) {
                         Load10X_Spatial(paste0(
                           "/home/mbatiuk/projects/human_schizo/MB6-MB23/primary/spatial/counts/Allen_genome/", smpl, "/outs/"), 
                           slice = smpl)})
names(visium_allen_plot) <- names(visium_allen_plot) %>% gsub("MB18-2", "MB18.2", .)

#Import layer labels with filtered out iffy areas
layers.filter.list <- lapply(
  setNames(c("MB11", "MB12", "MB14", "MB15", "MB18.2", "MB22", "MB23", "MB7"),
  c("MB11", "MB12", "MB14", "MB15", "MB18.2", "MB22", "MB23", "MB7")),
  function(x) {
    df <- read.csv(paste0("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/layers_filtered/", x, ".csv"), 
             header = T)
    setNames(df[, 2], df[, 1])})
```

``` r
## add layers
visium_allen_plot <- lapply(setNames(c("MB11", "MB12", "MB14", "MB15", "MB18.2", "MB22", "MB23", "MB7"),
  c("MB11", "MB12", "MB14", "MB15", "MB18.2", "MB22", "MB23", "MB7")),
    function(x) {
   #Also set white matter as exclude   
   visium_allen_plot[[x]]$Layers <- layers.filter.list[[x]][names(Idents(visium_allen_plot[[x]]))] %>% as.vector %>%
     replace_na("Exclude") %>% gsub("^WM$", "Exclude", .) %>% setNames(names(Idents(visium_allen_plot[[x]])))
   visium_allen_plot[[x]]})

## Subset visium seurat objects to remove NAs - no layer annotation, and WM
visium_allen_plot <- lapply(visium_allen_plot, function(x) {subset(x, subset = Layers == "Exclude", invert = T)})
```

``` r
## Add stereoscope predictions
visium_allen_plot <- lapply(setNames(names(visium_allen_plot), names(visium_allen_plot)), function(x) {
  anno.med <- stereo.med
  anno.med$MB18.2 <- anno.med$`MB18-2`
  anno.med$`MB18-2` <- NULL
  anno.high <- stereo.high
  anno.high$MB18.2 <- anno.high$`MB18-2`
  anno.high$`MB18-2` <- NULL
  visium_allen_plot[[x]][["Stereo_MED"]] <- anno.med[[x]][names(Idents(visium_allen_plot[[x]])), ] %>% t %>% CreateAssayObject
  visium_allen_plot[[x]][["Stereo_HIGH"]] <- anno.high[[x]][names(Idents(visium_allen_plot[[x]])), ] %>% t %>% CreateAssayObject
  visium_allen_plot[[x]]})
```

Plot Medium stereoscope res with seurat

``` r
visium_allen_plot <- lapply(visium_allen_plot, function(x) {
  DefaultAssay(x) <- "Stereo_MED"
  x})

lots.seur.stereo.MED <- lapply(setNames(names(visium_allen_plot), names(visium_allen_plot)),
                                      function(x) {
                        SpatialFeaturePlot(visium_allen_plot[[x]], 
                                           features = rownames(visium_allen_plot[[x]][["Stereo_MED"]]), 
                                           pt.size.factor = 1, 
                                           ncol = 2, 
                                           crop = F)})
```

same for HIGH annotation

``` r
visium_allen_plot <- lapply(visium_allen_plot, function(x) {
  DefaultAssay(x) <- "Stereo_HIGH"
  x})

plots.seur.stereo.HIGH <- lapply(setNames(names(visium_allen_plot), names(visium_allen_plot)),
                                      function(x) {
                                        SpatialFeaturePlot(visium_allen_plot[[x]], 
                                                           features = rownames(visium_allen_plot[[x]][["Stereo_HIGH"]]),
                                                           pt.size.factor = 1,
                                                           ncol = 2,
                                                           crop = F)})
```

Plot layers based on MB11 sample

``` r
#set default assay
DefaultAssay(visium_allen_plot$MB11) <- "Stereo_HIGH"
MB11.stereo.plot.list <- 
  append(list(Layers = SpatialPlot(visium_allen_plot$MB11, 
                                   group.by = "Layers", 
                                   pt.size.factor = 1.6, 
                                   crop = T, 
                                   cols = palette_45[1:6],
                                   stroke = 0) + 
  labs(subtitle = "Layers") + theme(plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold"), 
                            legend.text = element_text(size = 11), legend.position = "bottom", 
                            legend.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "mm"),
          legend.spacing.x = unit(0, "mm")) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1))),
  lapply(setNames(hsub_order_allen_10x %>% gsub("_", "-", .), hsub_order_allen_10x), function(subt) {
                      SpatialPlot(visium_allen_plot[["MB11"]], features = subt, pt.size.factor = 1.6, crop = T, alpha = c(0.5, 4), stroke = 0) +
  labs(subtitle = subt %>% gsub("-", "_", .)) +
                        theme(plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold"), 
                                legend.text = element_text(size = 11), legend.position = "bottom", 
                                legend.title = element_blank(), legend.key.width = unit(10, "mm"),
                              legend.key.height = unit(3, "mm"),
                              plot.margin = unit(c(0, 0, 0, 0), "mm"))})) %>%
  cowplot::plot_grid(plotlist=., ncol=5)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/mb11.jpg" width="60%" style="display: block; margin: auto;" />

### 3\. Cell proportions

code is duplicated because of existance of 2 cohorts

``` r
layers.filter.list0 <- lapply(
  setNames(c("MB11", "MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"),
  c("MB11", "MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")),
  function(x) {
    df <- read.csv(paste0("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/schizo_visium/layers_filtered/", x, ".csv"), 
             header = T)
    setNames(df[, 2], df[, 1])})

layers.filter.list1 <- lapply(
  setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"),
  c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52")),
  function(x) {
    df <- read.csv(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/layeranno/", x, "_MB.csv"), 
             header = T)
    setNames(df[, 2], df[, 1])})

layers.filter.list2 <- c(layers.filter.list1,layers.filter.list0)
```

MEDIUM

``` r
stereo.med1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                       y <-        read.delim(paste0("/home/viktor_petukhov/Data/Schizophrenia/deconv_v2/stereoscope_res/",
                                                    x, "_cm.tsv"), row.names = 1)
                       #filter low quality spots
                        y[rownames(y) %in% names(layers.filter.list0[[x]]), ]})

stereo.med1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                       y <- stereo.med1[[x]]
                       #change rownames
                       colnames(y) <- colnames(y) %>% gsub("\\.", "_", x = .)
                       y})

stereo.med2 <- lapply(setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"), 
                                c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52")), function(x) {
                       y <- read.delim(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/deconv2/",
                                                    x, ".tsv"), row.names = 1)
                       #filter low quality spots
                        y[rownames(y) %in% names(layers.filter.list1[[x]]), ]})

stereo.med2 <- lapply(setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"), 
                                c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52")), function(x) {
                       y <- stereo.med2[[x]]
                       #change rownames
                       colnames(y) <- colnames(y) %>% gsub("\\.", "_", x = .)
                       y})

stereo.med <- c(stereo.med1, stereo.med2)
```

HIGH

``` r
stereo.high1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                       y <-        read.delim(paste0("/d0/home/viktor_petukhov/Data/Schizophrenia/deconv_v2/stereoscope_res_high/",
                                                    x, "_cm.tsv"), row.names = 1)
                       #filter low quality spots
                         y[rownames(y) %in% names(layers.filter.list0[[x]]), ]})

stereo.high1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                       y <- stereo.high1[[x]]
                       #change rownames
                        colnames(y) <- colnames(y) %>% gsub("\\.", "_", x = .)
                        y})

stereo.high2 <- lapply(setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"), 
                                c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52")), function(x) {
                       y <- read.delim(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/deconv/",
                                                    x, ".tsv"), row.names = 1)
                       #filter low quality spots
                         y[rownames(y) %in% names(layers.filter.list1[[x]]), ]})

stereo.high2 <- lapply(setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"), 
                                c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52")), function(x) {
                       y <- stereo.high2[[x]]
                       #change rownames
                        colnames(y) <- colnames(y) %>% gsub("\\.", "_", x = .)
                        y})

stereo.high <- c(stereo.high1, stereo.high2)
```

load order

``` r
hsub_order <- readRDS("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/hsub_order.rds")
msub_order <- readRDS("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/msub_order.rds")
hsub_order_allen_10x <- factor(c(hsub_order[-c(length(hsub_order), length(hsub_order) -1)], 
                            "Astro_L1_FGFR3_SERPINI2", "Astro_L1_6_FGFR3_AQP1", "Astro_L1_6_FGFR3_PLCG1", 
                            "Oligo_L2_6_OPALIN_FTH1P3", "Oligo_L3_6_OPALIN_ENPP6", "Oligo_L2_6_OPALIN_MAP6D1",
                            "Oligo_L5_6_OPALIN_LDLRAP1", "OPC_L1_6_PDGFRA_COL20A1", "Micro_L1_6_TYROBP_CD74", 
                            "Endo_L2_5_NOSTRIN_SRGN", "VLMC_L1_5_PDGFRA_COLEC12"), 
                           
                           levels = c(hsub_order[-c(length(hsub_order), length(hsub_order) -1)], 
                             "Astro_L1_FGFR3_SERPINI2", "Astro_L1_6_FGFR3_AQP1", "Astro_L1_6_FGFR3_PLCG1", 
                            "Oligo_L2_6_OPALIN_FTH1P3", "Oligo_L3_6_OPALIN_ENPP6", "Oligo_L2_6_OPALIN_MAP6D1",
                            "Oligo_L5_6_OPALIN_LDLRAP1", "OPC_L1_6_PDGFRA_COL20A1", "Micro_L1_6_TYROBP_CD74", 
                            "Endo_L2_5_NOSTRIN_SRGN", "VLMC_L1_5_PDGFRA_COLEC12"),
                           ordered = T)
msub_order_allen_10x <- factor(c(msub_order[-c(length(msub_order), length(msub_order) -1)], 
                            "Astro", "Oligo", "OPC", "Micro_PVM", "Endo", "VLMC"), 
                           levels = c(msub_order[-c(length(msub_order), length(msub_order) -1)], 
                             "Astro", "Oligo", "OPC", "Micro_PVM", "Endo", "VLMC"),ordered = T)
```

Filter WM - not well represented in Scz and has variable size - skew in
bulk analysis. Make data frames

``` r
coh <- setNames(names(stereo.med),names(stereo.med))

library(data.table)
coh2 <- setNames(c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"),
  c("MB6", "MB9", "MB19", "MB21", "MB50", "MB52"))


c(MB6="SCZ", MB50="SCZ", MB52="SCZ", 
                          MB9="CNT", MB19="CNT", MB21="CNT")

samplegroups2 <- list(
  control = c('MB7', 'MB9', 'MB11', 'MB13', 'MB15', 'MB16', 'MB17', 'MB19', 'MB20', 'MB21', 'MB18-2'),
  schizo = c('MB6', 'MB10', 'MB12', 'MB14', 'MB22', 'MB23','MB8-2', 'MB50', 'MB52')
)

diseasef <- as.factor(setNames(rep(names(samplegroups2),unlist(lapply(samplegroups2,length))),unlist(samplegroups2)))

library(Matrix)
#MED
predict.stereo.med <- 
  lapply(coh, function(x) {
  t(stereo.med[[x]])[, names(layers.filter.list2[[x]][!layers.filter.list2[[x]] == "WM"])] %>% 
   apply(1, mean)
  }
) %>% as.data.frame %>% data.frame(subtypes = factor(rownames(.), levels = msub_order_allen_10x), .)
predict.stereo.med.melt <- melt(predict.stereo.med, id.vars = "subtypes", variable.name = "samples", 
                          value.name = "Mean probability")
#add condition
predict.stereo.med.melt$condition <- diseasef[predict.stereo.med.melt$samples %>% as.vector]
#HIGH
predict.stereo.high <- 
  lapply(coh, function(x) {
  t(stereo.high[[x]])[, names(layers.filter.list2[[x]][!layers.filter.list2[[x]] == "WM"])] %>% 
   apply(1, mean)
  }
) %>% as.data.frame %>% data.frame(subtypes = factor(rownames(.), levels = hsub_order_allen_10x), .)
predict.stereo.high.melt <- melt(predict.stereo.high, id.vars = "subtypes", variable.name = "samples", 
                          value.name = "Mean probability")
#add condition
predict.stereo.high.melt$condition <- diseasef[predict.stereo.high.melt$samples %>% as.vector]
```

LAYERS stereo data frames

``` r
#MED LAYERS
stereo.med.layers1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                   y <- read.delim(paste0("/home/viktor_petukhov/Data/Schizophrenia/deconv_v2/stereoscope_res/",
                                                    x, "_cm.tsv"), row.names = 1) %>% t
                   #correct annotations                 
                   rownames(y) <- rownames(y) %>% gsub("\\.", "_", x = .)                       
                   #filter low quality spots
                   y[, colnames(y) %in% names(layers.filter.list0[[x]])] %>%
                    data.frame(subtypes = factor(rownames(.), levels = msub_order_allen_10x), .) %>%
                    melt(id.vars = "subtypes", variable.name = "barcodes", value.name = "probability") %>%
                    cbind(., sample = rep(x, dim(.)[1])) %>% cbind(., condition = diseasef[rep(x, dim(.)[1])]) %>%
                    cbind(., layer = `[`(layers.filter.list2[[x]], .[["barcodes"]])) %>%
                     group_by(., subtypes, sample, condition, layer) %>% 
                     summarize_at(., "probability", list(mean_probability = mean, spot_number = length))})


stereo.med.layers2 <- lapply(coh2 , function(x) {
                       y <- read.delim(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/deconv2/",
                                                    x, ".tsv"), row.names = 1) %>% t
                       #correct annotations   
                       rownames(y) <- rownames(y) %>% gsub("\\.", "_", x = .)   
                       #filter low quality spots
                       y[, colnames(y) %in% names(layers.filter.list1[[x]])] %>%
                        data.frame(subtypes = factor(rownames(.), levels = msub_order_allen_10x), .) %>%
                        melt(id.vars = "subtypes", variable.name = "barcodes", value.name = "probability") %>%
                        cbind(., sample = rep(x, dim(.)[1])) %>% cbind(., condition = diseasef[rep(x, dim(.)[1])]) %>%
                        cbind(., layer = `[`(layers.filter.list2[[x]], .[["barcodes"]])) %>%
                         group_by(., subtypes, sample, condition, layer) %>% 
                         summarize_at(., "probability", list(mean_probability = mean, spot_number = length))})

stereo.med.layers<- c(stereo.med.layers1, stereo.med.layers2)
stereo.med.layers.comb <- bind_rows(stereo.med.layers)

#Filter out areas/layers with less than 5 spots
stereo.med.layers.comb.filt <- stereo.med.layers.comb[stereo.med.layers.comb$spot_number > 5, ]

#Filter out WM
stereo.med.layers.comb.filt <- stereo.med.layers.comb.filt[
  !stereo.med.layers.comb.filt$layer == "WM",]
```

``` r
#HIGH LAYERS
stereo.high.layers1 <- lapply(setNames(c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7"), 
                                c("MB11","MB12", "MB14", "MB15", "MB18-2", "MB22", "MB23", "MB7")), function(x) {
                   y <- read.delim(paste0("/d0/home/viktor_petukhov/Data/Schizophrenia/deconv_v2/stereoscope_res_high/",
                                                    x, "_cm.tsv"), row.names = 1) %>% t
                   #correct annotations    
                   rownames(y) <- rownames(y) %>% gsub("\\.", "_", x = .)    
                   #filter low quality spots
                   y[, colnames(y) %in% names(layers.filter.list0[[x]])] %>%
                    data.frame(subtypes = factor(rownames(.), levels = hsub_order_allen_10x), .) %>%
                    melt(id.vars = "subtypes", variable.name = "barcodes", value.name = "probability") %>%
                    cbind(., sample = rep(x, dim(.)[1])) %>% 
                     cbind(., condition = diseasef[rep(x, dim(.)[1])]) %>%
                     cbind(., layer = `[`(layers.filter.list0[[x]], .[["barcodes"]])) %>%
                     group_by(., subtypes, sample, layer, condition) %>% 
                     summarize_at(., "probability", list(mean_probability = mean, spot_number = length))})

stereo.high.layers2 <- lapply(coh2, function(x) {
                   y <- read.delim(paste0("/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/deconv/",
                                                    x, ".tsv"), row.names = 1) %>% t
                   #correct annotations      
                   rownames(y) <- rownames(y) %>% gsub("\\.", "_", x = .)                 
                   #filter low quality spots
                   y[, colnames(y) %in% names(layers.filter.list1[[x]])] %>%
                      data.frame(subtypes = factor(rownames(.), levels = hsub_order_allen_10x), .) %>%
                      melt(id.vars = "subtypes", variable.name = "barcodes", value.name = "probability") %>%
                      cbind(., sample = rep(x, dim(.)[1])) %>% cbind(., condition = diseasef[rep(x, dim(.)[1])]) %>%
                      cbind(., layer = `[`(layers.filter.list2[[x]], .[["barcodes"]])) %>%
                     group_by(., subtypes, sample, layer, condition) %>% 
                     summarize_at(., "probability", list(mean_probability = mean, spot_number = length))})

stereo.high.layers <- c(stereo.high.layers1,stereo.high.layers2)
stereo.high.layers.comb <- bind_rows(stereo.high.layers)

#Filter out areas/layers with less than 5 spots
stereo.high.layers.comb.filt <- stereo.high.layers.comb[stereo.high.layers.comb$spot_number > 5, ]

#Filter out WM
stereo.high.layers.comb.filt <- stereo.high.layers.comb.filt[
  !stereo.high.layers.comb.filt$layer == "WM",]
```

helper funcs

``` r
null.lenth.signif <- function(x, y) {
  if (signif(x)[pastes(y)][x[pastes(y)] <= 0.05 & !is.na(x[pastes(y)])] %>% length > 0) {
    signif(x)[pastes(y)][x[pastes(y)] <= 0.05 & !is.na(x[pastes(y)])]
} else {
  NULL}}

pastes <- function(x) {
  c(paste0(x, ".L1"), paste0(x, ".L2"), paste0(x, ".L3"), paste0(x, ".L4"), paste0(x, ".L5"),
     paste0(x, ".L6")
    )}

null.length <- function(vect, x, y) {
  if (x[pastes(y)][x[pastes(y)] <= 0.05 & !is.na(x[pastes(y)])] %>% length > 0) {
  vect[x[pastes(y)] <= 0.05 & !is.na(x[pastes(y)])]
} else {
  NULL}}
```

``` r
library(ggplot2)
library(ggsignif)
palette_45_2 <- readRDS("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/palette_45_2.rds")
palette_45 <- readRDS("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/palette_45.rds")

list.stereo.high.plot <-
    lapply(setNames(hsub_order_allen_10x, hsub_order_allen_10x), function(subt) {
    ggplot(data = 
             stereo.high.layers.comb.filt[stereo.high.layers.comb.filt[["subtypes"]] == as.character(subt), ], 
            aes(x = layer, y = mean_probability, dodge = condition, fill = condition)) +
      geom_point(aes(y = mean_probability, color = condition), 
                 position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.5), 
                   size = 1.1, alpha = 0.3, show.legend = F) +
      geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.3, show.legend = F, lwd = 0.4) +
      geom_signif(annotations = null.lenth.signif(wilcox.stereo.layer.high.vec.adj, subt),
                 xmin = null.length(c(0.75:5.75), wilcox.stereo.layer.high.vec.adj, subt), 
                 xmax = null.length(c(1.25:6.25), wilcox.stereo.layer.high.vec.adj, subt), 
                 y_position = null.length(c(rep(max(stereo.high.layers.comb.filt[stereo.high.layers.comb.filt[["subtypes"]] == 
                                                as.character(subt), ]$mean_probability)*1.1, 6)), wilcox.stereo.layer.high.vec.adj, subt), textsize = 7, vjust = 0.5) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 11),
              axis.title.y = element_text(size = 11),
              legend.title = element_blank(),
              legend.text = element_text(size = 11),
              plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold"),
              legend.position = "bottom",
              text = element_text(),
              plot.margin = unit(c(2.5, 13, 2.5, 1), "mm")) +
              labs(subtitle = subt, x = NULL, y = "Mean fraction") +
        scale_color_manual(values = palette_45_2[c(8, 15)]) +
        scale_fill_manual(values = palette_45_2[c(8, 15)]) +
        ylim(min = -0.001, max = max(
          stereo.high.layers.comb.filt[stereo.high.layers.comb.filt[["subtypes"]] == 
                                                as.character(subt), ]$mean_probability)*1.25)})
figS.stereo.high.plot <- cowplot::plot_grid(plotlist = list.stereo.high.plot, ncol = 5)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/meanfracvisium.jpg" width="60%" style="display: block; margin: auto;" />

### 4\. DE and GO analysis

``` r
#Read in Visium matrices
samples1 <-  read_rds("/d0-mendel/home/viktor_petukhov/Copenhagen/Schizophrenia20/cache/spatial_samples.rds")
samples2 <- readRDS("/d0/home/kdragicevic/RProjects/schizo/spatial/spatial_samples2.rds")

deconv_path1 <- "/home/viktor_petukhov/Data/Schizophrenia/deconv_v2/stereoscope_res/"
for(n in names(samples1)) {
  samples1[[n]]$deconv_probs <- read.table(paste0(deconv_path1, n, "_cm.tsv"))
  cm <- samples1[[n]]$cm
  samples1[[n]]$mit_frac <- Matrix::colSums(cm[grep("^MT-", rownames(cm)),]) / Matrix::colSums(cm)
}
deconv_path2 <- "/d0/home/kdragicevic/RProjects/schizo/spatial/prepdata/deconv2/"
for(n in names(samples2)) {
  samples2[[n]]$deconv_probs <- read.table(paste0(deconv_path2, n, ".tsv"))
  cm <- samples2[[n]]$cm
  samples2[[n]]$mit_frac <- Matrix::colSums(cm[grep("^MT-", rownames(cm)),]) / Matrix::colSums(cm)
}

samples <- c(samples1, samples2)
annot_spat <- lapply(samples, `[[`, "annotation") %>% Reduce(c, .)
```

``` r
#Create conos and cacoa object for DE
con_spat <- lapply(samples, `[[`, "p2") %>% Conos$new()
cao <- Cacoa$new(con_spat, 
                 ref.level="Ctr", 
                 target.level="Scz", 
                 sample.groups=setNames(c(rep("Ctr", 7),rep("Scz", 7)),unlist(samp_per_cond_spat, use.names = FALSE)), 
                 cell.groups=(annot_spat), 
                 n.cores=5)
de_spat_annot <- cao$estimatePerCellTypeDE(n.cores = 20)
```

``` r
#GO
library(org.Hs.eg.db)
org <- org.Hs.eg.db

library(pbapply)
go_datas <- c("BP", "CC", "MF") %>% setNames(., .) %>%
  pblapply(function(n) clusterProfiler:::get_GO_data(org.Hs.eg.db::org.Hs.eg.db, n, "ENTREZID") %>%
           as.list() %>% as.environment())

gos_spat <- list(annot = calculateGos(de_spat_annot, go_datas, n.top.genes=500))
gos_spat_dfs <- lapply(gos_spat, lapply, extractAllGODf)
gos_spat_clust <- lapply(gos_spat, lapply, clusterGos, n.clusters=20, max.pval=0.05)
```

``` r
cols <- list(up=circlize::colorRamp2(c(0, 1.3), c("grey98", "red")),
             down=circlize::colorRamp2(c(0, 1.3), c("grey98", "blue")))
```

\#UP GO

``` r
plotPValueHeatmap(gos_spat_clust$annot$up$summary, cols$up)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/goupvisium.jpg" width="60%" style="display: block; margin: auto;" />

``` r
pal <- setNames(RColorBrewer::brewer.pal(6, name = "Set2"), paste0("L", seq(1,6)))
upgo <- split.data.frame(gos_spat_dfs$annot$up, gos_spat_dfs$annot$up$Type)
upgo <- lapply(upgo, head,10)
upgo <- data.table::rbindlist(upgo)
upgo <- upgo[order(-p.adjust),]
upgo$Description %<>% make.unique()
upgo$Description %<>% factor(., levels=.)
upgo$p.adjust %<>% log10() %>% {. * -1}

ggplot(upgo, aes(p.adjust, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x="-log10(adj. P)", y="", fill="Layers - spatial data", title="Top 10 GO terms for up-regulated DE genes (spatial)") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = pal[upgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(upgo$Description), ".", fixed= TRUE), "[[", 1),upgo$Description))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/goup10.jpg" width="60%" style="display: block; margin: auto;" />

\#DOWN GO

``` r
plotPValueHeatmap(gos_spat_clust$annot$down$summary, cols$down)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/godownvisium.jpg" width="60%" style="display: block; margin: auto;" />

``` r
downgo <- split.data.frame(gos_spat_dfs$annot$down, gos_spat_dfs$annot$down$Type)
downgo <- lapply(downgo, head,10)
downgo <- data.table::rbindlist(downgo)
downgo <- downgo[order(-p.adjust),]
downgo$Description %<>% make.unique()
downgo$Description %<>% factor(., levels=.)
downgo$p.adjust %<>% log10() %>% {. * -1}

ggplot(downgo, aes(p.adjust, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x="-log10(adj. P)", y="", fill="Layers - spatial data", title="Top 10 GO terms for down-regulated DE genes (spatial)") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = pal[downgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(downgo$Description), ".", fixed= TRUE), "[[", 1),downgo$Description))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/godown10.jpg" width="60%" style="display: block; margin: auto;" />

### snRNAseq and Visium overlap

``` r
gos_sc <- goslW
go_df <- names(gos_sc$high$up) %>% setNames(., .) %>% 
  lapply(function(n) as_tibble(gos_sc$high$up[[n]]@result) %>% mutate(CellType=n)) %>% 
  plyr::rbind.fill()

gos_sc_filt <- gos_sc %>% lapply(lapply, function(x) x) %>% .["high"] %>% .[[1]]
gos_sc_aggr <- gos_sc_filt %>% names() %>% lapply(function(dir) {
  gos_sc_filt[[dir]] %>% names() %>%
    lapply(function(type) {
      gos_sc_filt[[dir]][[type]]@result %>% 
        filter(pvalue < 0.05) %>% 
        mutate(Type=type, Direction=dir) %>% 
        {if(nrow(.) < 4) . else .[1:3,]}
    })
}) %>% 
  do.call(c, .) %>% .[sapply(., nrow) > 0] %>% do.call(rbind, .) %>%
  .[order(.$p.adjust, decreasing=T),] %>% split(.$Direction) %>% lapply(function(df) {
  df$Description %<>% make.unique() %>% factor(., levels=.)
  df$p.adjust %<>% log10() %>% {. * -1}
  df
})

gos_spat_dfs_filt <- lapply(gos_spat_dfs, function(dfs) {
  names(dfs) %>% setNames(., .) %>% 
    lapply(function(n) {
      filter(dfs[[n]], Description %in% gos_sc_aggr[[n]]$Description) %>% 
        mutate(p.adjust=p.adjust(pvalue, method="BH"))
    })
})
```

``` r
plotPValueDfHeatmap(gos_spat_dfs_filt$annot$up, cols$up, p.threshold = 0.05) + ggtitle("Upregulated GO terms; Visium:snRNAseq overlap")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/overlapup.jpg" width="60%" style="display: block; margin: auto;" />

``` r
upgo <- split.data.frame(gos_spat_dfs_filt$annot$up, gos_spat_dfs_filt$annot$up$Type)

upgo <- lapply(upgo, head,3)
upgo <- data.table::rbindlist(upgo)
upgo <- upgo[order(-p.adjust),]
upgo$Description %<>% make.unique()
upgo$Description %<>% factor(., levels=.)
upgo$p.adjust %<>% log10() %>% {. * -1}
library(ggplot2)
upgobar <- ggplot(upgo, aes(p.adjust, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw()+
  theme(plot.title.position = "plot")+
  labs(x="-log10(adj. P)", y="", fill="Layers - spatial data", title="Top 3 upregulated GO terms; Visium:snRNA-seq overlap") + 
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = pal[upgo$Type %>% rev %>% unique]) +
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(upgo$Description), ".", fixed= TRUE), "[[", 1),upgo$Description))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/goupv.jpg" width="60%" style="display: block; margin: auto;" />

``` r
downgo <- split.data.frame(gos_spat_dfs_filt$annot$down, gos_spat_dfs_filt$annot$down$Type)

downgo <- lapply(downgo, head,3)
downgo <- data.table::rbindlist(downgo)
downgo <- downgo[order(-p.adjust),]
downgo$Description %<>% make.unique()
downgo$Description %<>% factor(., levels=.)
downgo$p.adjust %<>% log10() %>% {. * -1}

downgobar <- ggplot(downgo, aes(p.adjust, Description, fill=Type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw()+
  theme(plot.title.position = "plot")+
  labs(x="-log10(adj. P)", y="", fill="Layers - spatial data", title="Top 3 downregulated GO terms; Visium:snRNA-seq overlap") +  
  geom_vline(xintercept = 1.3, linetype="dotted") + scale_fill_manual(values = pal[downgo$Type %>% rev %>% unique]) + 
    scale_y_discrete(labels=setNames(sapply(strsplit(as.character(downgo$Description), ".", fixed= TRUE), "[[", 1),downgo$Description))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/godownv.jpg" width="60%" style="display: block; margin: auto;" />
