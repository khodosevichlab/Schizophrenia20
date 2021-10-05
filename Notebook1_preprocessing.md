Notebook \#1 - loading objects and preprocessing data, checking metadata
================

<style type="text/css">
.height {
  white-space: pre !important;
  overflow-y: scroll !important;
  height: 50vh !important;
}
</style>

### CONOS AND PAGODA PLOTS

``` r
library(pagoda2)
library(parallel)
library(conos)
library(Matrix)
library(ggplot2)
library(cowplot)
library(magrittr)
library(tidyverse)
library(tidyr)


samplegroups <- list(
  control = c('MB7', 'MB9', 'MB11', 'MB13', 'MB15', 'MB16', 'MB17', 'MB19', 'MB21', 'MB18-2', "MB51", "MB53", "MB55", "MB57"),
  schizo = c('MB6', 'MB10', 'MB12', 'MB14', 'MB22', 'MB23','MB8-2', "MB54", "MB56")
)
diseasef <- as.factor(setNames(rep(names(samplegroups),unlist(lapply(samplegroups,length))),unlist(samplegroups)))
```

### 1\. Loading count matrices from h5 files.

``` r
path <- '/home/mbatiuk/projects/human_schizo/primary/cohort2/counts/' #path to counts
x <- list.files(path,pattern='MB.*')

sn <- function(x) {setNames(x,x)} #helper function

cdl <- mclapply(sn(x),function(pat) {
  cd <- Seurat::Read10X_h5(paste(path, pat,'/outs/filtered_feature_bc_matrix.h5', sep = ''))
  colnames(cd) <- paste(pat, colnames(cd), sep = '_')
  cd
},mc.cores=30)
```

### 2\. Preprocessing and integration of data.

``` r
cdlf <- lapply(cdl,pagoda2::gene.vs.molecule.cell.filter,min.cell.size=500,plot=T)
```

count matrices were previously analyzed by “scrublet” in order to
identify potential doublets

### 3\. Applying scrublet filter.

``` r
scrubletf <- readRDS("scrubletf.RDS") #reading in result from scrublet which was saved as RDS previously
cdlf <- lapply(cdlf2, function(x) x[, scrubletf[colnames(x)] < 0.3])  #applying the filter; treshold is 0.3
```

### 4\. Apply pagoda2 processing.

``` r
cdlf.p2 <- mclapply(cdlf,function(m) basicP2proc(m,n.cores=1,
                                                 get.largevis = F,
                                                 get.tsne = T,
                                                 min.transcripts.per.cell = 500,
                                                 min.cells.per.gene = 0,
                                                 n.odgenes = 2e3,
                                                 make.geneknn = F),
                    mc.cores=30)
```

### 5\. Creating conos object from pagoda objects.

``` r
con <- Conos$new(cdlf.p2, n.cores = 30)
```

### 6\. Building UMAP embedding.

``` r
con$buildGraph(k=15, k.self=5, space='PCA', ncomps=30, n.odgenes=3000, matching.method='mNN', metric='angular', verbose=TRUE)
con$findCommunities(r=0.5)
con$embedGraph(method='UMAP')
con$plotGraph(alpha=0.1, size=0.1)
```

### 7\. Read in annotation and plot graph with it.

``` r
tans <- readRDS("tans.RDS") #contains both high level and medium level annotation; contains both clean (without clusters "Other" and "Glia") and full annotation
high.pal <- readRDS("high.pal.RDS") #palette for celltypes
med.pal <- readRDS("med.pal.RDS") #palette for celltypes

Figure1.B <- con$plotGraph(alpha=0.1,size=0.1,
              groups=tans$high,
              font.size=c(3,5),
              raster=T,
              plot.na=T,
              palette=high.pal)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig1b.jpg" width="40%" style="display: block; margin: auto;" />

### 8\. Estimate marker genes and plot marker genes of interest

\!\!\! since the number if cells is huge we downsampled the cells for
each cell type and estimated afterwards so that the session does not
crash

both made for medium and high annotation

``` r
dpal <- c('control'='dodgerblue1','schizo'='indianred1')
dpalf <- function(n) dpal
celldiseasef <- con$getDatasetPerCell()
celldiseasef <- setNames(diseasef[celldiseasef],names(celldiseasef))


#medium
inh.cells <- names(tans$med.clean)[tans$med.clean %in% grep("^Id|^Pval|^Sst|^Vip",levels(tans$med.clean),val=T)]
max.cells <- 1e3; # maximum number of cells per population
set.seed(0)
inh.fac <- tans$high.clean
inh.fac <- inh.fac[names(inh.fac) %in% inh.cells]; inh.fac <- droplevels(inh.fac);
inh.sampled.cells <- unlist(tapply(names(inh.fac),inh.fac,function(x) sample(x,min(length(x),max.cells))))
inh.de <- con$getDifferentialGenes(groups=inh.fac[inh.sampled.cells],n.cores=30,append.auc=TRUE,z.threshold=0,upregulated.only=T)

#Figure 3.D
pm <- pp <- plotDEheatmap(con,
                          inh.fac[inh.sampled.cells],
                          inh.de,
                          n.genes.per.cluster = 10,
                          show.gene.clusters=T,
                          column.metadata=list(samples=con$getDatasetPerCell(),
                                               disease=celldiseasef),
                          order.clusters = T,
                          use_raster=T,
                          raster_device = "CairoPNG", 
                          column.metadata.colors = list(clusters=high.pal,
                                                        disease=dpal),
                          row.label.font.size = 8)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig3d.jpg" width="40%" style="display: block; margin: auto;" />

``` r
sampled.cells <- unlist(tapply(names(tans$high),tans$high,function(x) sample(x,min(length(x),max.cells))))
high.de <- con$getDifferentialGenes(groups=tans$high[sampled.cells],
                                   n.cores=30,
                                   append.auc=TRUE,
                                   z.threshold=0,
                                   upregulated.only=T)

#Figure 1.D
ph <- plotDEheatmap(con,tans$high,
                    med.de,
                    n.genes.per.cluster = 10 ,
                    show.gene.clusters=T,
                    column.metadata = list(samples = con$getDatasetPerCell(),
                                           disease = celldiseasef),
                    order.clusters = T,
                    use_raster=T,
                    raster_device = "CairoPNG", 
                    column.metadata.colors = list(clusters = high.pal, 
                                                  disease = dpal),
                    row.label.font.size = 8)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig1d.jpg" width="40%" style="display: block; margin: auto;" />

### 9\. Plot chosen marker genes on UMAP panel

``` r
#Figure 1.C
plot.list.markers <- lapply(X = c("GAD1", "SATB2", "CUX2", "RORB", "FEZF2", "THEMIS", "VIP",
                                  "ID2", "SST", "PVALB", "SLC1A3", "MBP"),
                               function(x) {
  con$plotGraph(gene = x, 
                    title = x, 
                    alpha= 1, 
                    size = 0.3, 
                    plot.na = F, 
                    raster = T, 
                    raster.dpi = T) +
  theme_void() +
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(size = 15, 
                                  hjust = 0.5, 
                                  face = "bold"), 
        legend.position="none")})

plot.markers <- grid.arrange(grobs = plot.list.markers, ncol = 4)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig1c.jpg" width="60%" style="display: block; margin: auto;" />

### 10\. plot estimated cell position (layers)

``` r
#read transfered annotations from Allen data
annot_tranfer_mtg_highres <-  readRDS("/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/mtg_high_transf_annotations.rds")
#load Allen annotation
annotation_allen <- fread(file = "~/projects/human_schizo/MB6-MB23/2ndary/Allen_31.01.2020/sample_annotations.csv")

#Figure 1.E
paper.layers.mtg.high <-
 ggplot(data.frame(Subtypes = annot_tranfer_mtg_highres$labels[annotation_allen$sample_name[annotation_allen$sample_name %in% names(annot_tranfer_mtg_highres$labels)]] %>%
        `[`(., .!= "Other") %>%
         factor(levels = c("L2_CUX2_LAMP5_MARCH1", "L2_CUX2_LAMP5_PDGFD", "L2_3_CUX2_FREM3_UNC5D", "L2_3_CUX2_FREM3_SV2C",
            "L3_CUX2_PRSS12", "L4_RORB_SCHLAP1_MET", "L4_RORB_SCHLAP1_MME", "L4_RORB_SCHLAP1_ARHGAP15", "L4_5_FEZF2_LRRK1",
            "L5_FEZF2_ADRA1A", "L5_6_FEZF2_TLE4_ABO", "L5_6_FEZF2_TLE4_SCUBE1", "L5_6_FEZF2_TLE4_HTR2C", "L5_6_THEMIS_SEMA3A",
            "L5_6_THEMIS_NTNG2", "ID2_LAMP5_NMBR", "ID2_LAMP5_CRH", "ID2_LAMP5_NOS1", "ID2_PAX6", "ID2_NCKAP5", 
            "VIP_ABI3BP", "VIP_TYR", "VIP_RELN", "VIP_SEMA3C", "VIP_SSTR1", "VIP_CRH", "PVALB_MEPE", "PVALB_SST", 
            "PVALB_CRH", "SST_TH", "SST_TAC3", "SST_CALB1", "SST_NPY", "SST_NOS1", "SST_STK32A", "Glia")),
        Layers = setNames(annotation_allen$cortical_layer_label, annotation_allen$sample_name)[annot_tranfer_mtg_highres$labels[annotation_allen$sample_name[annotation_allen$sample_name %in% names(annot_tranfer_mtg_highres$labels)]] %>%
        `[`(., .!= "Other") %>% names] %>%
         factor(levels = c("L6", "L5", "L4", "L3", "L2", "L1"))), aes(x = Subtypes, 
                                                                      y = Layers,  
                                                                      color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.4, alpha = 0.5, size = 0.1, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), 
             linetype = 2, 
             size = 1, 
             colour = "grey") +
  scale_color_manual(values = high.pal[levels(tans$high)]) +
  labs(title = "Position, MTG estimate, high resolution") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5, 
                                   size = 7), 
        plot.title = element_text(hjust = 0.5, 
                                  size = 15, 
                                  face = "bold"),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_blank())
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig1e.jpg" width="60%" style="display: block; margin: auto;" />

### METADATA TESTS

### 1\. Categorical covariates

``` r
metadata <- read.csv("~/RProjects/schizo/revision/meta.csv")
```

Now follows a very long list of different factors from metadata

``` r height
annotations_scz2 <- tans

phfactor <- setNames(metadata$pH_liquor,
                     metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

source("http://pklab.med.harvard.edu/rasmus/scRNA_helper.R")
mitpercfactor <- mitoFraction(con, species="human")

mitpercfactor <- mitpercfactor*100
#Sex
sexfactor <- setNames(metadata$Gender,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

#PMI
pmifactor <- setNames(metadata$PMI_hrs,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#colors PMI
color.PMI <-
    colorRamps::matlab.like2(n = (length((min(metadata$PMI_hrs) * 100) : (max(metadata$PMI_hrs) * 100)))) %>%
      setNames((min(metadata$PMI_hrs) * 100) : (max(metadata$PMI_hrs) * 100)) %>%
      `[`(., (metadata$PMI_hrs * 100) %>% unique %>% sort %>% as.character) %>% as.vector
#Source
sourcefactor <- setNames(metadata$Source,
                         metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
# Cause of death
deathfactor <- 
  setNames(metadata$Cause_of_death %>% factor(., levels = c(`[`(unique(.), unique(.) != "n.a."), "n.a.")),
           metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

#FACS
facsfactor <- setNames(metadata$FACS_10x_1st_day_date,
                       metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#% NeuN+
neunfactor <- setNames(metadata$`%_NeuN_positive`,
                       metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#colors NeuN
color.neun <-
      colorRamps::matlab.like2(n = (length((min(metadata$`%_NeuN_positive`) * 10) : 
                                             (max(metadata$`%_NeuN_positive`) * 10)))) %>%
      setNames((min(metadata$`%_NeuN_positive`) * 10) : (max(metadata$`%_NeuN_positive`) * 10)) %>%
      `[`(., (metadata$`%_NeuN_positive` * 10) %>% unique %>% sort %>% as.character) %>% as.vector %>% rev
#GEM kit lot
gemfactor <- setNames(metadata$`10x_GEM_kit_lot`,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#cDNA_clean_pre-amp date
preampfactor <- setNames(metadata$cDNA_cleanup_preAmp_date,
                         metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#library prep date
libdatefactor <- setNames(metadata$library_prep_date,
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Lib prep kit lot
  libfactor <- setNames(metadata$`10x_libr_prep_kit_batch`,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

#Color volume susp loaded

#cDNA concentration
cdnafactor <- setNames(metadata$`cDNA_concentration_ng/ul`,
                       metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color cDNA concentration
color.cdna <-
      colorRamps::matlab.like2(n = (length((min(metadata$`cDNA_concentration_ng/ul`) * 100) : 
                                             (max(metadata$`cDNA_concentration_ng/ul`) * 100)))) %>%
      setNames((min(metadata$`cDNA_concentration_ng/ul`) * 100) : (max(metadata$`cDNA_concentration_ng/ul`) * 100)) %>%
      `[`(., (metadata$cDNA_concentration * 100) %>% unique %>% sort %>% as.character) %>% as.vector %>% rev
#ng per library prep
amountlibfactor <- setNames(metadata$`Amount_(ng)_per_reaction`,
                            metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Index primer
indexfactor <- setNames(metadata$Index_primer,
                        metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Autopsy brain hypoxia signs
hypoxiafactor <- setNames(factor(metadata$Postmortem_signs_of_brain_hipoxia %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .),
                                 levels = c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Colors Y/N + n.a.
color.yn <- c('dodgerblue1', 'indianred1', "grey92")
#Tangles
tanglesfactor <- setNames(factor(metadata$`Dementia_tangles_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                 
                                 levels = c("No", "Low", "Med", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#color tangles
color.tangles <- c(colorRamps::matlab.like2(n = (metadata$`Dementia_tangles_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")
#Amyloid
amyloidfactor <- setNames(factor(metadata$`Dementia_amyloid_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                 levels = c("No", "Low", "Med", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color Amyloid
color.amyloid <- c(colorRamps::matlab.like2(n = (metadata$`Dementia_amyloid_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")
#Atherosclerosis
atherofactor <- setNames(factor(metadata$`Brain_vasculature_atherosclerosis_N(No)_Low_Mid_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                levels = 
                                  c("No", "Low", "Med", "High", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color Amyloid
color.athero <- c(colorRamps::matlab.like2(n = (metadata$`Brain_vasculature_atherosclerosis_N(No)_Low_Mid_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")
#Lewy bodies
lewyfactor <- setNames(factor(metadata$`Parkinsonism_Lewy_bodies_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .), 
                              levels = 
                                  c("No", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color lewy
color.lewy <- c("dodgerblue1", "grey92")
#Lithium treatment
lithiumfactor <- setNames(factor(metadata$Lithium %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#neuroleptic duration
neuroleptfactor <- setNames(factor(metadata$Neuroleptic_treatment_duration_years, levels = 
                                     c(0, 0.003, 0.16, 2, 10, 30, 40, "n.a.")),
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color neuroleptic duration years
color.neurolept <-
       c(colorRamps::matlab.like2(n = (metadata$Neuroleptic_treatment_duration_years %>% 
                                                   unique %>% length) -1), "grey92")
#typical neuroleptics
typicalfactor <- setNames(factor(metadata$Typical_neuroleptics %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#atypical neuroleptics
atypicalfactor <- setNames(factor(metadata$Atypical_neuroleptics %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
 
#tranquilizers benzodiazepines
benzofactor <- setNames(factor(metadata$Tranquilizers_benzodiazepines %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#tranquilizers barbiturates
barbifactor <- setNames(factor(metadata$Tranquilizers_barbiturates %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Antidepressants tricyclic
tricycfactor <- setNames(factor(metadata$Antidepressant_tricyclic %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Antidepressants SIRT/NE inhib
sirtfactor <- setNames(factor(metadata$`Antidepressant_serotonin_norepinephrine-transporter_inhibitors` %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

#Chemo tamoxifen
tamoxifactor <- setNames(factor(metadata$Chemotherapy_tamoxifen %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Chemo antiprolif
proliffactor <- setNames(factor(metadata$Chemotherapy_anti_proliferative %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Radiotherapy
radiofactor <- setNames(factor(metadata$Radiotherapy %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Cancer
cancerfactor <- setNames(factor(metadata$Cancer %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Onset of psychiatry
onsetfactor <- setNames(factor(metadata$First_schizophrenia_symptoms_age, 
                               levels = metadata$First_schizophrenia_symptoms_age %>% 
                                 as.vector %>% `[`(. != "n.a.") %>% `[`(. != "No") %>% as.numeric %>%
                                           unique %>% sort %>% c(., "No", "n.a.")),
                               metadata$Identifier)[
                                 as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
color.onset <-
      colorRamps::matlab.like2(n = 
                                  (length((min(metadata$First_schizophrenia_symptoms_age %>% as.vector %>% 
                                                 `[`(. != "n.a.")  %>% `[`(. != "No") %>% as.numeric)) : 
                                    (max(metadata$First_schizophrenia_symptoms_age %>% as.vector %>% 
                                           `[`(. != "n.a.")  %>% `[`(. != "No") %>% as.numeric))))) %>%
  
  
      setNames((min(metadata$First_schizophrenia_symptoms_age %>% as.vector %>% `[`(. != "n.a.") %>% `[`(. != "No") %>%
                      as.numeric)) : 
                 (max(metadata$First_schizophrenia_symptoms_age %>% as.vector %>% `[`(. != "n.a.") %>% `[`(. != "No")  %>%
                        as.numeric))) %>%
      `[`(., ((metadata$First_schizophrenia_symptoms_age %>% as.vector %>% `[`(. != "n.a.") %>% `[`(. != "No") %>%
                 as.numeric)) %>% 
            unique %>% sort %>% as.character) %>% as.vector %>% c(., "thistle2", "grey92")
#Hyperthension
hyperfactor <- setNames(factor(metadata$Hypertension %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Depression
deprfactor <- setNames(factor(metadata$Depression %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Malnutrition/weight loss
malnutrfactor <- setNames(factor(metadata$Malnutrition_weight_loss %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#breathing problems
breathfactor <- setNames(factor(metadata$Breathing_problems %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Opioid
opioidfactor <- setNames(factor(metadata$Opioid_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Dopamine
dopaminefactor <- setNames(factor(metadata$Dopaminergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Serotonine
serotoninefactor <- setNames(factor(metadata$Serotoninergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Glu-ergic
glufactor <- setNames(factor(metadata$Glutamatergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#GABA
gabafactor <- setNames(factor(metadata$GABAergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Histamine
histfactor <- setNames(factor(metadata$Histaminergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Muscarine
muscarfactor <- setNames(factor(metadata$Muscarinic_acetycholinergic__modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Nicotinic
nicofactor <- setNames(factor(metadata$Nicotinic_acetylcholine_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Adrenergic
adrenfactor <- setNames(factor(metadata$Adrenergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Adenosin
adenfactor <- setNames(factor(metadata$Adenosinergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Glycine
glycinefactor <- setNames(factor(metadata$Glycinergic_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Na
nafactor <- setNames(factor(metadata$Na_channels_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#K
kfactor <- setNames(factor(metadata$K_channels_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Ca
cafactor <- setNames(factor(metadata$Ca_channels_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Proton pump inhibit
protonfactor <- setNames(factor(metadata$Proton_pump_inhibitors_PPI %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#MAO inhibitors
maofactor <- setNames(factor(metadata$MAO_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Mitoch electron transport
mitofactor <- setNames(factor(metadata$Mitochondrial_electron_transport_modulators %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Smoking
smokefactor <- setNames(factor(metadata$Nicotine_abuse %>%
                                   gsub("N", "No", .) %>% gsub("Y", "Yes", .), 
                                 levels = 
                                  c("No", "Yes", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
```

``` r
#editing the sex factor additionaly
Toch=function(conosCluster){
  cluster2=as.character(conosCluster)
  names(cluster2)=names(conosCluster)
  return(cluster2)
}
Female = Toch(sexfactor)
Female[Female=='F']='Yes'
Female[Female=='M']='No'
Female=as.factor(Female)
```

``` r
#list of factors to plot:

vec=list('Signs of postmortem brain hypoxia'=hypoxiafactor,
         'Female'=Female,
         'Typical neuroleptics'=typicalfactor,
         'Atypical neuroleptics'=atypicalfactor,
         'Benzodiazepines'=benzofactor,
         'Barbiturates'=barbifactor,
         'Tricyclic'=tricycfactor,
         'SIRT NE transporter inhibitors antidepressants'=sirtfactor,
         'Tamoxifen chemotherapy'=tamoxifactor,
         'Anti-proliferative chemotherapy'=proliffactor,
         'Radiotherapy'=radiofactor,
         'Cancer'=cancerfactor,
         'Hypertension'=hyperfactor,
         'Depression'=deprfactor,
         'Malnutrition'=malnutrfactor,
         'Breathing problems'=breathfactor,
         'Opioid modulators'=opioidfactor,
         'Dopaminergic modulators'=dopaminefactor,
         'Serotoninergic modulators'=serotoninefactor,
         'Glutamatergic modulators'=glufactor,
         'GABAergic modulators'=gabafactor,
         'Histaminergic modulators'=histfactor,
         'Muscarine receptors modulators'=muscarfactor,
         'Nicotine receptors modulators'=nicofactor,
         'Adrenergic modulators'=adrenfactor,
         'Purinergic modulators'=adenfactor,
         'Glycine receptors modulators'=glycinefactor,
         'Sodium channels modulators'=nafactor,
         'Potassium channels modulators'=kfactor,
         'Calcium channels modulators'=cafactor,
         'Proton pump modulators'=protonfactor,
         'MAO modulators'=maofactor,
         'Mitochondrial e transport modulators'=mitofactor,
         'Smoking'=smokefactor)
```

remove category with only one sample

``` r
nn=names(vec)
nn=nn[!nn %in% c("Barbiturates",'Tamoxifen chemotherapy')]
anoSample=con$getDatasetPerCell()
cname=names(annotations_scz2$high.clean)


sn=function(x) { names(x) <- x; return(x); }
res=lapply(sn(names(vec)),function(x){
  temp=vec[[x]]
  temp=Toch(temp[cname]) %>% .[.!='n.a.'] %>% as.factor()
  tsample.groups <- tapply(anoSample[names(temp)],temp,unique)
  tsample.groups <- as.factor(setNames(rep(names(tsample.groups),unlist(lapply(tsample.groups,length))),unlist(tsample.groups)))
  tsample.groups
})
```

``` r
res=res[nn]

#get samples which are Scz
schizo <- names(cao$sample.groups[which(cao$sample.groups == "Scz")])

fraction=Toch(anoSample)
index=fraction %in% schizo

fraction[index]='schizo'
fraction=as.factor(fraction)
sample.groups <- tapply(anoSample,fraction,unique)
sample.groups <- as.factor(setNames(rep(names(sample.groups),unlist(lapply(sample.groups,length))),unlist(sample.groups)))
sample.groups
n1=sample.groups

tabl=lapply(sn(names(res)),function(x) {
  n2=res[[x]]
  tab=table(n2,n1[names(n2)])
  tab
})

datNo=do.call(rbind,lapply(tabl,function(x) x['No',]))
datYes=do.call(rbind,lapply(tabl,function(x) x['Yes',]))
```

``` r
library(reshape2)
d1=setNames(melt(datYes), c('rows', 'vars', 'values'))
d1$type='No'
d2=setNames(melt(datNo), c('rows', 'vars', 'values'))
d2$type='Yes'
dd=rbind(d1,d2)
#dd
dd$name=paste(dd$rows,dd$type)


library(ggplot2)
p <- ggplot(data=dd, aes(x=name, y=values, fill=vars)) +
geom_bar(stat="identity", color="black")+ coord_flip()+xlab('')+ylab('')
fraction.pal <- c("control" = "#980e5c", "schizo" = "#435790")
p=p+scale_fill_manual(values=fraction.pal)
p <- p + ylim(0, 21)
```

``` r
pvalue.list = lapply(tabl,function(x) fisher.test(data.matrix(x))$p.value) 
pvalue.l=do.call(c,pvalue.list)
dd=data.frame('condition'=names(pvalue.l),'p.value'=pvalue.l)

any(pvalue.l < 0.01) # all NS

#Extended data figure 1.A.
p <- p + geom_signif(y_position = c(20.5), xmin = c(seq(0.5,62.5,2)), 
              xmax = c(seq(0.5,62.5,2)), annotation = c(rep("NS",length(pvalue.l))),
              tip_length = 0, lwd = 0) + guides(fill = guide_legend(title = "")) + theme(legend.position="bottom") + ggtitle("Categorical covariates") + theme(plot.title = element_text(face = "bold", size = 15, hjust = -0.2))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig1a.jpg" width="40%" style="display: block; margin: auto;" />

### 2\. Numeric covariates

``` r
metadata$Diagnosis %<>% gsub("Ctr", "control", .) %>% gsub("Scz", "schizo", .)
library(magrittr)
meta.list <- 
    list(age = metadata[, c("Identifier", "Age", "Diagnosis")],
         pmi = metadata[, c("Identifier", "PMI_hrs", "Diagnosis")],
         ph = metadata[, c("Identifier", "pH_liquor", "Diagnosis")],
         duration = metadata[, c("Identifier", "Neuroleptic_treatment_duration_years", "Diagnosis")]) %>% lapply(., function(x) {
      split(x, f = x$Diagnosis)
}) %>% unlist(recursive = F)

#remove MB8

meta.list <-lapply(meta.list, function(x) {
 x <- x[x$Identifier != "MB8", ]
 x  
})

#clean up
#remove n.a
meta.list$ph.control <- meta.list$ph.control[meta.list$ph.control$pH_liquor != "n.a.", ]
meta.list$ph.schizo <- meta.list$ph.schizo[meta.list$ph.schizo$pH_liquor != "n.a.", ]

meta.list$ph.control$pH_liquor %<>% as.numeric()
meta.list$ph.schizo$pH_liquor %<>% as.numeric()

meta.list$duration.control <- meta.list$duration.control[meta.list$duration.control$Neuroleptic_treatment_duration_years != "n.a.",]
meta.list$duration.schizo <- meta.list$duration.schizo[meta.list$duration.schizo$Neuroleptic_treatment_duration_years != "n.a.",]
meta.list$duration.schizo <-na.omit(meta.list$duration.schizo)
meta.list$duration.control$Neuroleptic_treatment_duration_years %<>% as.numeric()
meta.list$duration.schizo$Neuroleptic_treatment_duration_years %<>% as.numeric()
meta.list$duration.schizo <- na.omit(meta.list$duration.schizo)
meta.list <- lapply(meta.list, unique) #to ensure unique rows in case something got duplicated


#shapiro test
shapiro.meta <- lapply(meta.list, function(x) {tryCatch(shapiro.test(unlist(x[,2])), error = function(e) NULL)})
#p<0.05 - non-normal distribution; p>0.05 - normal distribution
```

``` r
meta.list2 <- 
    list(age = metadata[, c("Identifier", "Age", "Diagnosis")],
         pmi = metadata[, c("Identifier", "PMI_hrs", "Diagnosis")],
         ph = metadata[, c("Identifier", "pH_liquor", "Diagnosis")],
         duration = metadata[, c("Identifier", "Neuroleptic_treatment_duration_years", "Diagnosis")])

#remove MB8

meta.list2 <-
lapply(meta.list2, function(x) {
 x <- x[x$Identifier != "MB8", ]
 x})

#RENAME schizo to SZ
meta.list2 <-lapply(meta.list2, function(x) {
  x$Diagnosis <- x$Diagnosis %>% gsub("schizo", "SZ", .)
  x})


meta.list2$ph <- meta.list2$ph[meta.list2$ph$pH_liquor != "n.a.", ]
meta.list2$ph$pH_liquor %<>% as.numeric()
meta.list2$duration <- meta.list2$duration[meta.list2$duration$Neuroleptic_treatment_duration_years != "n.a.", ]
meta.list2$duration$Neuroleptic_treatment_duration_years %<>% as.numeric()
meta.list2 <- lapply(meta.list2, unique)

#stat test
library(car)
levene.meta <- lapply(meta.list2, function(x) leveneTest(x[[colnames(x[,2])]] ~ x$Diagnosis))

levene.meta.p <-
lapply(levene.meta, function(x) {
  y <- x$`Pr(>F)`[1]
}) 
#unlist(levene.meta.p) <= 0.05
 
wilcox.meta <-lapply(meta.list2, function(x) {wilcox.test(x[[colnames(x[,2])]] ~ x$Diagnosis)})
ttest.meta <-lapply(meta.list2, function(x) {t.test(x[[colnames(x[,2])]] ~ x$Diagnosis, var.equal = T)})
```

make plots

age plot

``` r
library(ggplot2)
library(ggsignif)
library(dplyr)

palette_45_2 <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/palette_45_2.rds")
age.plot <-
  ggplot(data = meta.list2$age, aes(x = 1, y = Age, dodge = Diagnosis, fill = Diagnosis)) +
  geom_point(aes(y = Age, color = Diagnosis), position = position_jitterdodge(dodge.width = 2.3, jitter.width = 0.2), 
               size = 1.1, alpha = 0.5, show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.3) +
  geom_signif(annotations = signif(ttest.meta$age$p.value),
             xmin = 0.8, 
             xmax = 1.2, 
             y_position = max(meta.list2$age$Age)*1.07,
             textsize = 5,annotation = c("NS")) +
  stat_summary(fun = "mean", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.75), 
                 mapping = aes(group = Diagnosis), fill = "grey",
                 shape = 23) + 
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()
          ) +
          labs(x = NULL, y = "Age, years") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
    ylim(0, max(meta.list2$age$Age)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/ageplot.jpg" width="10%" style="display: block; margin: auto;" />

pmiplot

``` r
pmi.plot <-
  ggplot(data = meta.list2$pmi, aes(x = 1, y = PMI_hrs, dodge = Diagnosis, fill = Diagnosis)) +
  geom_point(aes(y = PMI_hrs, color = Diagnosis), 
             position = position_jitterdodge(dodge.width = 2.3, jitter.width = 0.2), 
               size = 1.1, alpha = 0.5, show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.3) +
  geom_signif(annotations = signif(wilcox.meta$pmi$p.value),
             xmin = 0.8, 
             xmax = 1.2, 
             y_position = max(meta.list2$pmi$PMI_hrs)*1.07,
             textsize = 5,annotation = c("NS")) +
  stat_summary(fun = "mean", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.75), 
                 mapping = aes(group = Diagnosis), fill = "grey",
                 shape = 23) + 
  theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) + 
  labs(x = NULL, y = "PMI, hrs") +
  scale_color_manual(values = palette_45_2[c(8, 15)]) +
  scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(0, max(meta.list2$pmi$PMI_hrs)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/pmiplot.jpg" width="10%" style="display: block; margin: auto;" />

phplot

``` r
ph.plot <-
  ggplot(data = meta.list2$ph, aes(x = 1, y = pH_liquor, dodge = Diagnosis, fill = Diagnosis)) +
  geom_point(aes(y = pH_liquor, color = Diagnosis), 
             position = position_jitterdodge(dodge.width = 2.3, jitter.width = 0.2), 
               size = 1.1, alpha = 0.5, show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.3) +
  geom_signif(annotations = signif(wilcox.meta$ph$p.value),
             xmin = 0.8, 
             xmax = 1.2, 
             y_position = max(meta.list2$ph$pH_liquor)*1.07,
             textsize = 5,annotation = c("NS")) +
   stat_summary(fun = "mean", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.75), 
                 mapping = aes(group = Diagnosis), fill = "grey",
                 shape = 23) + 
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) +
          labs(x = NULL, y = "pH") +
  scale_color_manual(values = palette_45_2[c(8, 15)]) +
  scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(0, max(meta.list2$ph$pH_liquor)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/phplot.jpg" width="10%" style="display: block; margin: auto;" />

durationplot

``` r
meta.list2$duration <- na.omit(meta.list2$duration)
durat.plot <-
  ggplot(data = meta.list2$duration, aes(x = 1, 
          y = Neuroleptic_treatment_duration_years, dodge = Diagnosis, fill = Diagnosis)) +
  geom_point(aes(y = Neuroleptic_treatment_duration_years, color = Diagnosis), 
             position = position_jitterdodge(dodge.width = 2.3, jitter.width = 0.2), 
               size = 1.1, alpha = 0.5, show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.3) +
  geom_signif(annotations = signif(wilcox.meta$duration$p.value),
             xmin = 0.8, 
             xmax = 1.2, 
             y_position = max(meta.list2$duration$Neuroleptic_treatment_duration_years)*1.07,
             textsize = 5, annotation = c("*")) +
  stat_summary(fun = "mean", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.75), 
                 mapping = aes(group = Diagnosis), fill = "grey",
                 shape = 23) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) +
          labs(x = NULL, y = "Neuroleptics, years") +
  scale_color_manual(values = palette_45_2[c(8, 15)]) +
  scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(0, max(meta.list2$duration$Neuroleptic_treatment_duration_years)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/duratplot.jpg" width="10%" style="display: block; margin: auto;" />

``` r
#Extended data figure 1.B
library(patchwork)
covar.plot <-
  (age.plot + theme(plot.margin = unit(c(0, 0.8, 0, 0), "cm"))) + 
  (pmi.plot + theme(plot.margin = unit(c(0, 0.8, 0, 0), "cm"))) +
  (ph.plot + theme(plot.margin = unit(c(0, 0.8, 0, 0), "cm"))) +
  (durat.plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) +
  plot_layout(guides = 'collect', ncol = 4) &
  plot_annotation(subtitle = "Numeric covariates") & 
  theme(legend.position = "top", plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig1b.jpg" width="60%" style="display: block; margin: auto;" />

### 3\. Estimate number of genes per nucleus in each samples and number of nuclei in each sample

``` r
cms_scz_list <- cdl #raw matrices
cms_scz_list <- lapply(setNames(unlist(samplegroups), unlist(samplegroups)),
     function(sn) {
       paste0(sn, "_", colnames(cms_scz_list[[sn]])) %>%
         set_colnames(cms_scz_list[[sn]], .)})
```

``` r
#Filter out cells in scz samples and not present sannotation
cms_scz_list %<>% lapply(function(smpl) {
  smpl[, colnames(smpl) %in% names(tans$high)]})


#Count amount of expressed genes/cell >0
library(Matrix)
genes_per_nuc <- lapply(cms_scz_list, function(y){
  y@x <- rep(1, length(y@x))
  colSums(y)})

genes_per_nuc_vec <-
    setNames(unlist(genes_per_nuc),
    gsub("^MB[[:digit:]][[:digit:]]?-?[[:digit:]]?.", "", names(unlist(genes_per_nuc))))
```

``` r
#now for umis
umi_per_nuc <- lapply(cms_scz_list, function(y) colSums(y))

umi_per_nuc_vec <-
    setNames(unlist(umi_per_nuc),
    gsub("^MB[[:digit:]][[:digit:]]?-?[[:digit:]]?.", "", names(unlist(umi_per_nuc))))
```

``` r
genes.p.nuc.df <- data.frame(genes_expressed = genes_per_nuc_vec,
                             sample = gsub("(^MB[[:digit:]][[:digit:]]?-?[[:digit:]]?)(.*)", "\\1",
                                           names(genes_per_nuc_vec)),
                                           condition = diseasef[gsub("(^MB[[:digit:]][[:digit:]]?-?[[:digit:]]?)(.*)", "\\1",
                                           names(genes_per_nuc_vec))],
                             subtype = tans$high[names(genes_per_nuc_vec)],
                             umis_detected = umi_per_nuc_vec[names(genes_per_nuc_vec)]
)
genes.p.nuc.df$sample %<>% factor(levels = names(sort(diseasef)))
```

``` r
#load ordering vectors
hsub_order <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/hsub_order.rds")
msub_order <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/msub_order.rds")
genes.p.nuc.df$subtype %<>% factor(levels = hsub_order)
```

``` r
#get samples per subtype
samp.p.subt.df <-  genes.p.nuc.df %>% group_by(sample, subtype) %>% summarize(n_of_cells = n())

#filter subtypes with less than 5 or 10 cells
samp.p.subt.df <- samp.p.subt.df[samp.p.subt.df$n_of_cells >= 5,]
 
#summarize
samp.p.subt.df %<>% ungroup %>% group_by(subtype) %>% summarize(freq = n()/length(diseasef))
```

``` r
#get of subtypes per sample
subs.p.sample.df <-  genes.p.nuc.df %>% group_by(sample, subtype) %>% summarize(n_of_cells = n())

#filter subtypes with less than 5 or 10 cells
subs.p.sample.df <- subs.p.sample.df[subs.p.sample.df$n_of_cells >= 5, ]
 
#summarize
subs.p.sample.df %<>% ungroup %>% group_by(sample) %>% summarize(freq = n()/length(tans$high %>% unique))
subs.p.sample.df$condition <- diseasef[subs.p.sample.df$sample]

#nuclei per sample
nuc.p.sample.df <-  genes.p.nuc.df %>% group_by(sample, subtype) %>% summarize(count = n()) %>% 
      group_by(subtype) %>% mutate(fraction = count/sum(count))

nuc.p.condition.df <- genes.p.nuc.df %>% group_by(condition, subtype) %>% summarize(count = n()) %>% 
      group_by(subtype) %>% mutate(fraction = count/sum(count))
```

``` r
nucs.p.sample.plot <-
  ggplot(data = genes.p.nuc.df, aes(x = sample, fill = condition)) +
  geom_bar(width = 0.7) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 15),
          
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "# of nuclei") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2c.jpg" width="60%" style="display: block; margin: auto;" />

``` r
levels(genes.p.nuc.df$sample) <- names(cao$sample.groups %>% sort)

genes.p.nuc.plot <-  ggplot(data = genes.p.nuc.df, aes(x = sample, y = genes_expressed, fill = condition)) +
  rasterise(geom_sina(size = 0.1, alpha = 0.02, aes(color = condition), show.legend = T, scale = F), dpi = 600) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1, width = 0.2, color = "black", fill = "white", coef = 0) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Genes/nucleus") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2b.jpg" width="60%" style="display: block; margin: auto;" />

``` r
#Extended data figure 2. B and C
library(cowplot)
p <- plot_grid(plotlist = list(genes.p.nuc.plot,nucs.p.sample.plot), labels = c("B", "C"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2bc.jpg" width="80%" style="display: block; margin: auto;" />

### 4\. Plot doublets on UMAP embedding

``` r
#Import data
scrubletf <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/scrubletf.rds")

#Filter out cells not present in scrubletf
cms_scz_list %<>% lapply(function(smpl) {
  smpl[, colnames(smpl) %in% names(scrubletf)]
    })

p2.list.scrub <- lapply(cms_scz_list, basicP2proc, n.cores = 10, n.odgenes = 3000, min.cells.per.gene = 0,
                  min.transcripts.per.cell = 0, get.largevis = F, get.tsne = T, make.geneknn = F)
con.scrub <- Conos$new(p2.list.scrub,n.cores = 10)
con.scrub$buildGraph(k = 15, k.self = 5, 
               space = "PCA",
               ncomps = 30, n.odgenes = 3000, matching.method = "mNN", metric = "angular",
               score.component.variance = T, verbose = T)
con.scrub$embedGraph(method = "UMAP",n.cores = 70, min.dist = 0.1, spread = 15)
```

``` r
#create gg object for getting cluster annotations
gg_annot <- con.scrub$plotGraph(groups = annotations_scz2$med)

doublet.plot <-con.scrub$plotGraph(alpha = 0.5, 
                                   size = 0.05, 
                                   colors = scrubletf, 
                                   show.legend = T, 
                                   font.size = 4, 
                                   raster = T, 
                                   raster.dpi = 300) +
  #get cluster annotations 
  gg_annot$layers[[2]] +
  scale_size_continuous(range = c(4, 4), guide = F, trans='identity') +
  theme_void() +
  labs(subtitle = "Doublet score", color = "Score") +
  theme(text = element_text(size = 15),
        plot.subtitle = element_text(size = 15, 
                                     hjust = 0.5, 
                                     face = "bold"), 
        legend.position = c(1, 0.9))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2d.jpg" width="40%" style="display: block; margin: auto;" />

``` r
doublet.plot.excluded <- con.scrub$plotGraph(alpha = 0.5, 
                                             size = 0.05, 
                                             groups = setNames(cut(scrubletf, 
                                                                   br = c(0, 0.25, 1), 
                                                                   include.lowest = T, 
                                                                   labels = c("Keep", "Exclude")), 
                                                               names(scrubletf)), 
                                             show.legend = F, 
                                             mark.groups = F, 
                                             font.size = 4, 
                                             raster = T, raster.dpi = 300) +
  #get cluster annotations 
  gg_annot$layers[[2]] +
  scale_size_continuous(range = c(4, 4), guide = F, trans='identity') +
  theme_void() +
  labs(subtitle = "Doublet exclusion") +
  theme(text = element_text( size = 15),
        plot.subtitle = element_text(size = 15, 
                                     hjust = 0.5, 
                                     face = "bold"), 
        legend.position = "top",
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, 
                                                   size = 3))) +
  scale_color_manual(values = c("#779876", "magenta"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2e.jpg" width="40%" style="display: block; margin: auto;" />

``` r
#Extended data figure 2. D and E
p <-
      (doublet.plot + theme(plot.margin = unit(c(0.5, 3, 0 ,0), "cm")) + 
         doublet.plot.excluded + theme(plot.margin = unit(c(0.5, 0, 0 ,1), "cm"))) + 
  plot_layout(heights = c(1, 3.5))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig2de.jpg" width="60%" style="display: block; margin: auto;" />

### 5\. plot rest of raw sample statistics (extended figure 3)

Genes per nucleus

``` r
genes.p.nuc.subt.plot <-
  ggplot(data = genes.p.nuc.df, aes(x = subtype, y = genes_expressed)) +
  rasterise(geom_sina(size = 0.1, alpha = 0.2, aes(color = subtype), show.legend = F, scale = F), 
            dpi = 600) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1, width = 0.3, color = "black", fill = "white", coef = 0) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Genes per\nnucleus") +
    scale_color_manual(values = high.pal[levels(genes.p.nuc.df$subtype)])
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3a.jpg" width="80%" style="display: block; margin: auto;" />

UMI per nucleus

``` r
umi.p.nuc.subt.plot <-
  ggplot(data = genes.p.nuc.df, aes(x = subtype, y = umis_detected)) +
  rasterise(geom_sina(size = 0.1, alpha = 0.2, aes(color = subtype), show.legend = F, scale = F), 
            dpi = 600) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1, width = 0.3, color = "black", fill = "white", coef = 0) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),,
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "UMIs per\nnucleus") +
    scale_color_manual(values = high.pal[levels(genes.p.nuc.df$subtype)])
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3b.jpg" width="80%" style="display: block; margin: auto;" />

Samples with more than 5 nuclei

``` r
percent <- c("0%", "25%", "50%", "75%", "100%")
sampl.p.subt.plot <-
  ggplot(data = samp.p.subt.df, aes(x = subtype, y = freq, fill = subtype)) +
  geom_col(show.legend = F, width = 0.9) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Samples\nwith ≥ 5 nuclei") +
    scale_color_manual(values = high.pal[levels(genes.p.nuc.df$subtype)]) +
    scale_fill_manual(values = high.pal[levels(genes.p.nuc.df$subtype)]) +
  scale_y_continuous(labels = percent, expand = c(0, 0))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3c.jpg" width="80%" style="display: block; margin: auto;" />

Subtypes with more than 5 nuclei

``` r
subs.p.sample.df$condition <- recode_factor((subs.p.sample.df$sample), !!!sort(setNames(as.character(samplegroups), names(samplegroups))))

sub.p.sample.plot <-
  ggplot(data = subs.p.sample.df, aes(x = sample, y = freq, fill = condition)) +
  geom_col(show.legend = T, width = 0.9) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
         legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Subtypes\nwith ≥ 5 nuclei") +
    scale_fill_manual(values = palette_45_2[c(8, 15)], labels  = c("control", "schizophrenia")) +
    scale_color_manual(values = palette_45_2[c(8, 15)],  labels  = c("control", "schizophrenia")) +
  scale_y_continuous(labels = percent, expand = c(0, 0))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3d.jpg" width="80%" style="display: block; margin: auto;" />

Fraction of nuclei of each sample per cluster

``` r
nuc.p.sample.df$condition <- recode_factor((nuc.p.sample.df$sample), !!!sort(setNames(as.character(samplegroups), names(samplegroups))))

library(pals)
library(RColorBrewer)
cols <- sample(rainbow(n = 27, s = 0.6, v = 0.8))
nuc.p.sample.df
nuc.p.sample.plot <-
  ggplot(data = nuc.p.sample.df, aes(x = subtype, y = fraction, fill = sample)) +
  geom_col(show.legend = T, width = 0.9, position = "stack") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
         legend.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Fraction of nuclei") +
    scale_fill_manual(values =cols) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 3,override.aes = list(size = 0.2)), shape = guide_legend(override.aes = list(size = 0.2)))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3e.jpg" width="80%" style="display: block; margin: auto;" />

Fraction of nuclei (ctr vs scz) per cluster

``` r
nuc.p.cond.plot <-
  ggplot(data = nuc.p.condition.df, aes(x = subtype, y = fraction, fill = condition)) +
  geom_col(show.legend = T, width = 0.9, position = "stack") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
         legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Fraction of nuclei") +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
  scale_y_continuous(labels = percent, expand = c(0, 0))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3f.jpg" width="80%" style="display: block; margin: auto;" />

``` r
#Extended data figure 3 (all plots)
p <-
      (genes.p.nuc.subt.plot /
      umi.p.nuc.subt.plot /
      sampl.p.subt.plot / 
      (sub.p.sample.plot + nuc.p.cond.plot + plot_layout(width = c(19, 37))) /
      nuc.p.sample.plot) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig3.jpg" width="60%" style="display: block; margin: auto;" />

### 6\. Panel of many UMAPS showing certain factors

``` r height
col27 <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
"#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
"#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
"#8A7C64", "#599861", "#F9998B")


plot.anno.high <- 
      con$plotGraph(groups = annotations_scz2$high,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, 
        raster = T, raster.dpi = 600) +
 
  theme_void() +
  labs(subtitle = "Subtypes, high resolution") +
  theme(plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.position="none") +
  scale_color_manual(values = high.pal)

plot.batch.samples <- 
      con$plotGraph(groups = con$getDatasetPerCell()[
        names(con$getDatasetPerCell()) %in% names(annotations_scz2$high)],
                    alpha = 0.03, size = 0.1, palette = col27,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Samples") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=2))

phfactor <- setNames(metadata$pH_liquor,
                     metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#colors pH
color.ph <-
       colorRamps::matlab.like2(n = 
                                  (length((min(metadata$pH_liquor %>% as.vector %>% `[`(. != "n.a.") %>% as.numeric) * 100) : 
                                    (max(metadata$pH_liquor %>% as.vector %>% `[`(. != "n.a.") %>% as.numeric) * 100)))) %>%
  
  
      setNames((min(metadata$pH_liquor %>% as.vector %>% `[`(. != "n.a.") %>% as.numeric) * 100) : 
                 (max(metadata$pH_liquor %>% as.vector %>% `[`(. != "n.a.") %>% as.numeric) * 100)) %>%
      `[`(., ((metadata$pH_liquor %>% as.vector %>% `[`(. != "n.a.") %>% as.numeric) * 100) %>% 
            unique %>% sort %>% as.character) %>% as.vector %>% rev %>% c(., "grey92")

library(viridis)
source("http://pklab.med.harvard.edu/rasmus/scRNA_helper.R")
mitpercfactor <- mitoFraction(con, species="human")
mitpercfactor <- mitpercfactor*100
plot.batch.mitperc <- 
      con$plotGraph(groups = as.factor(round(mitpercfactor,0)),
                    alpha = 0.05, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "% of mitochondrial genes") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12))+ guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=2)) + scale_color_manual(values = viridis(32))

plot.batch.ph <- 
      con$plotGraph(groups = phfactor,
                    alpha = 0.05, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Liquor pH") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
 scale_color_manual(values = color.ph)


setNames(names(diseasef), diseasef)
celldiseasef <- setNames(diseasef[as.vector(con$getDatasetPerCell())], names(con$getDatasetPerCell()))

plot.batch.disease <- 
      con$plotGraph(groups = celldiseasef,
                    alpha = 0.015, size = 0.01,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Condition") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = palette_45_2[c(8, 15)])

sexfactor <- setNames(metadata$Gender,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

plot.batch.sex <- 
      con$plotGraph(groups = sexfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Gender") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = c('cyan', 'magenta'))


#Age
agefactor <- setNames(metadata$Age,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
color.age <-
    colorRamps::matlab.like2(n = length(min(metadata$Age) : max(metadata$Age))) %>%
      setNames(min(metadata$Age) : max(metadata$Age)) %>%
      `[`(., metadata$Age %>% unique %>% sort %>% as.character) %>% as.vector

plot.batch.age <- 
      con$plotGraph(groups = agefactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Age") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=2)) +
  scale_color_manual(values = color.age)


#PMI
pmifactor <- setNames(metadata$PMI_hrs,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#colors PMI
color.PMI <-
    colorRamps::matlab.like2(n = (length((min(metadata$PMI_hrs) * 100) : (max(metadata$PMI_hrs) * 100)))) %>%
      setNames((min(metadata$PMI_hrs) * 100) : (max(metadata$PMI_hrs) * 100)) %>%
      `[`(., (metadata$PMI_hrs * 100) %>% unique %>% sort %>% as.character) %>% as.vector

plot.batch.pmi <- 
      con$plotGraph(groups = pmifactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "PMI") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=2)) +
  scale_color_manual(values = color.PMI)


sourcefactor <- setNames(metadata$Source,
                         metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
library(pals)
plot.batch.source <- 
      con$plotGraph(groups = sourcefactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = col27, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Brain bank") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))
 

deathfactor <- 
  setNames(metadata$Cause_of_death %>% factor(., levels = c(`[`(unique(.), unique(.) != "n.a."), "n.a.")),
           metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Colors cause of death
color.death <- c(col27[1:length(metadata$Cause_of_death %>% unique) - 1], "grey92")

plot.batch.death <- 
      con$plotGraph(groups = deathfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Cause of death") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = color.death)

facsfactor <- setNames(metadata$FACS_10x_1st_day_date,
                       metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
plot.batch.facs <- 
      con$plotGraph(groups = facsfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = col27, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "FACS date") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))
 
neunfactor <- setNames(metadata$`%_NeuN_positive`,
                       metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#colors NeuN
color.neun <-
      colorRamps::matlab.like2(n = (length((min(metadata$`%_NeuN_positive`) * 10) : 
                                             (max(metadata$`%_NeuN_positive`) * 10)))) %>%
      setNames((min(metadata$`%_NeuN_positive`) * 10) : (max(metadata$`%_NeuN_positive`) * 10)) %>%
      `[`(., (metadata$`%_NeuN_positive` * 10) %>% unique %>% sort %>% as.character) %>% as.vector %>% rev

plot.batch.neun <- 
      con$plotGraph(groups = neunfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "% NeuN+ nuclei during FACS") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=2)) +
  scale_color_manual(values = color.neun)

gemfactor <- setNames(metadata$`10x_GEM_kit_lot`,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

rpal <- function(n) {sample(rainbow(n))}

plot.batch.gem <- 
      con$plotGraph(groups = gemfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = rpal, 
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "GEM kit lot") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))
 
preampfactor <- setNames(metadata$cDNA_cleanup_preAmp,
                         metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

plot.batch.preamp <- 
      con$plotGraph(groups = preampfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = rpal, 
        raster = T,  raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "cDNA pre-amplification date") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))

libdatefactor <- setNames(metadata$library_prep_date,
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

plot.batch.libprep <- 
      con$plotGraph(groups = libdatefactor[
        names(con$getDatasetPerCell()) %in% names(annotations_scz2$high)],
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = rpal, 
        raster = T,  raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Library prep date") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))
 

libfactor <- setNames(metadata$`10x_libr_prep_kit_batch`,
                      metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

plot.batch.liblot <- 
      con$plotGraph(groups = libfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F, palette = rpal, 
        raster = T,  raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Library prep kit lot") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1))
 

amyloidfactor <- setNames(factor(metadata$`Dementia_amyloid_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                 levels = c("No", "Low", "Med", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#Color Amyloid
color.amyloid <- c(colorRamps::matlab.like2(n = (metadata$`Dementia_amyloid_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")

plot.batch.amyloid <- 
      con$plotGraph(groups = amyloidfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Amyloid") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = color.amyloid)


atherofactor <- setNames(factor(metadata$`Brain_vasculature_atherosclerosis_N(No)_Low_Mid_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                levels = 
                                  c("No", "Low", "Med", "High", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))

color.athero <- c(colorRamps::matlab.like2(n = (metadata$`Brain_vasculature_atherosclerosis_N(No)_Low_Mid_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")

plot.batch.athero <- 
      con$plotGraph(groups = atherofactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T, raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Atherosclerosis") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = color.athero)

tanglesfactor <- setNames(factor(metadata$`Dementia_tangles_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .), 
                                 
                                 levels = c("No", "Low", "Med", "n.a.")), 
                          metadata$Identifier)[as.vector(con$getDatasetPerCell()[names(annotations_scz2$high)])] %>%
  setNames(names(con$getDatasetPerCell()[names(annotations_scz2$high)]))
#color tangles
color.tangles <- c(colorRamps::matlab.like2(n = (metadata$`Dementia_tangles_N(no)_Low_Mid_High_Very_High` %>%
                                   gsub("N", "No", .) %>% gsub("Mid", "Med", .) %>% 
                                                   unique %>% length) -1), "grey92")

plot.batch.tangles <- 
      con$plotGraph(groups = tanglesfactor,
                    alpha = 0.03, size = 0.1,
                    font.size = 4, show.legend = F, plot.na = F, mark.groups = F,
        raster = T,  raster.dpi = 600) +
  theme_void() +
  labs(subtitle = "Protein tangles") +
  theme(legend.position = "right", 
        plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1)) +
  scale_color_manual(values = color.tangles)
```

``` r
#Extended data figure 4
fig.E4_1 <-
  (plot.anno.high + theme(plot.margin = unit(c(0, 1.5, 0, 2), "cm"))) + 
  (plot.batch.samples + theme(plot.margin = unit(c(0, 1.5, 0, 0), "cm"))) + 
  (plot.batch.ph + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) +
  (plot.batch.neun + theme(plot.margin = unit(c(1.4, 1.5, 0, 2), "cm"))) +
  (plot.batch.disease + theme(plot.margin = unit(c(1.4, 1.5, 0, 0), "cm"))) +
  (plot.batch.sex + theme(plot.margin = unit(c(1.4, 0, 0, 0), "cm"))) +
  (plot.batch.age + theme(plot.margin = unit(c(1.4, 1.5, 0, 2), "cm"))) +
  (plot.batch.pmi + theme(plot.margin = unit(c(1.4, 1.5, 0, 0), "cm"))) +
  (plot.batch.source + theme(plot.margin = unit(c(1.4, 0, 0, 0), "cm"))) +
  (plot.batch.death + theme(plot.margin = unit(c(1.4, 1.5, 0, 2), "cm"))) + 
  (plot.batch.facs + theme(plot.margin = unit(c(1.4, 1.5, 0, 0), "cm"))) +
  (plot.batch.preamp + theme(plot.margin = unit(c(1.4, 0, 0, 0), "cm"))) +
  (plot.batch.mitperc + theme(plot.margin = unit(c(1.4, 1.5, 0, 2), "cm"))) +
  (plot.batch.libprep + theme(plot.margin = unit(c(1.4, 1.5, 0, 0), "cm"))) +
  (plot.batch.liblot + theme(plot.margin = unit(c(1.4, 0, 0, 0), "cm"))) +
  (plot.batch.amyloid + theme(plot.margin = unit(c(1.4, 1.5, 0, 2), "cm"))) +
  (plot.batch.athero + theme(plot.margin = unit(c(1.4, 1.5, 0, 0), "cm"))) +
  (plot.batch.tangles + theme(plot.margin = unit(c(1.4, 0, 0, 0), "cm"))) +
  plot_layout(ncol = 3) + theme(plot.subtitle = element_text(size = 20), legend.text = element_text(size = 17))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig4.jpg" width="80%" style="display: block; margin: auto;" />
