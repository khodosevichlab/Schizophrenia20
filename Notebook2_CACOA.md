Notebook 2: Cacoa plots - cell proportions, expression shifts, cell
density, cell loadings
================

NOTE:

most of this code was made with “Cacoa”: “Case-Control Analysis of
scRNA-seq experiments”, a package not published officially yet but
available on github:
<https://github.com/kharchenkolab/cacoa/tree/master>

Please refer to the github repository of Cacoa for most functions in
this notebook: estimate expression shift magnitudes, estimate cell
density, estimate cell loadings, estimate per cell type de

``` r
library(cacoa)
library(conos)
library(Matrix)
library(ggplot2)
library(ggrastr)
library(Cairo)
library(sccore)
library(tidyverse)
library(tidyr)
library(data.table)
library(ggsignif)
```

### 1\. Create Cacoa object with high level annotation

``` r
samplegroups <- readRDS("samplegroups.RDS") #for samplegroups please refer to Notebook 1
con #conos object from Notebook1
tans #annotation from Notebook1
cao_high <- Cacoa$new(con, 
                      ref.level="Ctr", 
                      target.level="Scz", 
                      sample.groups=samplegroups, 
                      cell.groups=tans$high.clean,
                      cell.groups.palette = high.pal,
                      n.cores=30)
```

### 2\. Estimate cell proportions

``` r
#Subtype orders
msub_order <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/msub_order.rds")
hsub_order <- readRDS("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/ctx_layers_SCZ_ALLEN/rds_objects/hsub_order.rds")


df.cellfrac.high <- data.frame(cell_name = names(celldiseasef), 
                               condition = celldiseasef,  
                               high_clust = factor(highanno[names(con$getDatasetPerCell())], 
                                                   levels = hsub_order),
                               sample = con$getDatasetPerCell())

df.cellfrac.med <- data.frame(cell_name = names(celldiseasef), 
                              condition = celldiseasef, 
                              med_clust = factor(mediumanno[names(con$getDatasetPerCell())], 
                                                 levels = msub_order),
                              sample = con$getDatasetPerCell())
```

``` r
#NA come from extra scrublet filter
#Remove NAs
df.cellfrac.high <- df.cellfrac.high[!is.na(df.cellfrac.high$high_clust), ]
df.cellfrac.med <- df.cellfrac.med[!is.na(df.cellfrac.med$med_clust), ]
#Rename MB8 as MB8-2 - same sample in the end:
df.cellfrac.med$sample %<>% gsub("^MB8$", "MB8-2", .)
df.cellfrac.high$sample %<>% gsub("^MB8$", "MB8-2", .)
#Remove glia - CONTAMINATION DURING SORT - DISTRIBUTION STOCHASTIC BASED ON SORT GATES
df.cellfrac.med <- df.cellfrac.med[!(df.cellfrac.med$med_clust == "Glia"), ]
df.cellfrac.high <- df.cellfrac.high[!(df.cellfrac.high$high_clust == "Glia"), ]
df.cellfrac.high$high_clust %<>% droplevels
df.cellfrac.med$med_clust %<>% droplevels
#Summarize data
df.cellrac.summarized.high <- df.cellfrac.high %>% group_by(high_clust, sample) %>% summarize(n =  n()) %>% group_by(sample) %>% mutate(freq = n / sum(n))
df.cellrac.summarized.med <- df.cellfrac.med %>% group_by(med_clust, sample) %>% summarize(n =  n()) %>% group_by(sample) %>% mutate(freq = n / sum(n))
df.cellrac.summarized.med$condition <- diseasef[df.cellrac.summarized.med$sample]
df.cellrac.summarized.high$condition <- diseasef[df.cellrac.summarized.high$sample]
#Log-transform - explore possibility of t-test
df.cellrac.summarized.med %<>% mutate(freq_log = log1p(freq))
df.cellrac.summarized.high %<>% mutate(freq_log = log1p(freq))
#log(max(x+1) - x) transformation
df.cellrac.summarized.med %<>% mutate(freq_logmax = log(max(freq + 1) - freq))
df.cellrac.summarized.high %<>% mutate(freq_logmax = log(max(freq + 1) - freq))
#sqrt transform
df.cellrac.summarized.med %<>% mutate(freq_sqrt = sqrt(freq))
df.cellrac.summarized.high %<>% mutate(freq_sqrt = sqrt(freq))
#sqrt(max(x+1) - x) transform
df.cellrac.summarized.med %<>% mutate(freq_sqrtmax = sqrt(max(freq+1) - freq))
df.cellrac.summarized.high %<>% mutate(freq_sqrtmax = sqrt(max(freq+1) - freq))
#reciprocal transform
df.cellrac.summarized.med %<>% mutate(freq_recip = 1/freq)
df.cellrac.summarized.high %<>% mutate(freq_recip = 1/freq)
#1/(max(x+1) - x) reciprocal transform
df.cellrac.summarized.med %<>% mutate(freq_recipmax = 1/(max(freq+1) - freq))
df.cellrac.summarized.high %<>% mutate(freq_recipmax = 1/(max(freq+1) - freq))
```

``` r
#Define Significance function for changing p values to *
signif <- function(x) {
  y <- x
  y[x <= 0.001] <- "***"
  y[x <= 0.01 & x > 0.001] <- "**"
  y[x <= 0.05 & x > 0.01] <- "*"
  y[x > 0.05] <- "NS"
  y[is.nan(x)] <- NaN
  y[is.na(x) & !is.nan(x)] <- NA
  return(y)}
```

``` r
#High resolution
df.cellrac.summarized.high.list2 <- 
    split(df.cellrac.summarized.high, f = df.cellrac.summarized.high$high_clust)
df.cellrac.summarized.med.list2 <- 
    split(df.cellrac.summarized.med, f = df.cellrac.summarized.med$med_clust)
wilcox.high <-
  lapply(df.cellrac.summarized.high.list2, function(x) {wilcox.test(freq ~ condition, data = x)})
wilcox.med <-
  lapply(df.cellrac.summarized.med.list2, function(x) {wilcox.test(freq ~ condition, data = x)})
```

``` r
df.cellrac.summarized.med$condition %<>% gsub("^schizo.*", "schizophrenia", .)
wilcox.med.adj <-
lapply(wilcox.med, function(x) x$p.value) %>% unlist %>% 
  p.adjust(method = "BH") 
#wilcox.med.adj %>% `[`(., . <= 0.05)

wilcox.high.adj <-
lapply(wilcox.high, function(x) x$p.value) %>% unlist %>% 
  p.adjust(method = "BH") 
#wilcox.high.adj %>% `[`(., . <= 0.05)
```

``` r
cell.fraction.plot.med.adj <-
  ggplot(data = df.cellrac.summarized.med, aes(x = med_clust, 
                                               y = freq, 
                                               dodge = condition, 
                                               fill = condition)) +
  geom_point(aes(y = freq, 
                 color = condition), 
             position = position_jitterdodge(dodge.width = 0.7, 
                                             jitter.width = 0.5),
               size = 1.1, 
             alpha = 0.2, 
             show.legend = F) +
  geom_boxplot(notch = F, 
               outlier.shape = NA, 
               alpha = 0.3, 
               lwd = 0.4) +
  geom_signif(annotations = signif(wilcox.med.adj)[!wilcox.med.adj > 0.05],
             xmin = c(0.75:13.75)[!wilcox.med.adj > 0.05], 
             xmax = c(1.25:14.25)[!wilcox.med.adj > 0.05], 
             y_position = c(rep(sqrt(max(
                df.cellrac.summarized.med$freq)*0.9), 14))[!wilcox.med.adj > 0.05],
             textsize = 5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(), , plot.margin = unit(c(0.2,0.2,1,0.2), "cm")) +
          labs(x = NULL, y = "Fraction of nuclei") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  scale_y_sqrt(breaks = c(10e-4, 0.05, 0.1, 0.2, 0.3,0.4), 
               limits = c(0.00, max(df.cellrac.summarized.med$freq)*1.15), 
          labels = c("0.1%", "5%", "10%", "20%", "30%", "40%"), 
          expand=c(0,0)) +
  labs(subtitle = "Medium resolution")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig5a.jpg" width="80%" style="display: block; margin: auto;" />

``` r
cell.fraction.plot.high.adj <- 
  ggplot(data = df.cellrac.summarized.high, aes(x = high_clust, 
                                                y = freq, 
                                                dodge = condition, 
                                                fill = condition)) +
  geom_point(aes(y = freq, color = condition), 
             position = position_jitterdodge(dodge.width = 0.7, 
                                             jitter.width = 0.5), 
             size = 1.1, 
             alpha = 0.2, 
             show.legend = F) +
  geom_boxplot(notch = F, 
               outlier.shape = NA, 
               alpha = 0.3, 
               lwd = 0.4) +
  geom_signif(annotations = signif(wilcox.high.adj)[!wilcox.high.adj > 0.05],
             xmin = c(0.75:34.75)[!wilcox.high.adj > 0.05], 
             xmax = c(1.25:35.25)[!wilcox.high.adj > 0.05], 
             y_position = c(rep(sqrt(max(
                df.cellrac.summarized.med$freq)*0.9), 35))[!wilcox.high.adj > 0.05],
             textsize = 5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(), plot.margin = unit(c(0.2,0.2,1,0.2), "cm")) +
          labs(x = NULL, y = "Fraction of nuclei") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  scale_y_sqrt(breaks = c(10e-4, 0.05, 0.1, 0.2, 0.3, 0.4), 
               limits = c(0.00, max(
          df.cellrac.summarized.med$freq)*1.15), 
          labels = c("0.1%", "5%", "10%", "20%", "30%", "40%"), 
          expand=c(0,0)) +
  labs(subtitle = "High resolution")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig5b.jpg" width="80%" style="display: block; margin: auto;" />

``` r
#Extended data figure 5.A and B
#NOTE: Figure 2.B was made from the same data, just with reduced number of cell types
cell.fraction.plot.full.med.high.adj <-
cell.fraction.plot.med.adj + cell.fraction.plot.high.adj +
  plot_layout(guides = 'collect', width = c(15, 25)) & theme(legend.position = "top",
                      legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) & plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig5ab.jpg" width="80%" style="display: block; margin: auto;" />

### 3\. Cell density (Figure 3.A)

``` r
#Creating some vectors for easier dealing with ggplot layers in cacoa objects
highanno <- tans$high.clean
highnames <- highanno %>% unique
cdmediumnames <- c("L2", "L2", "L2", "L5_6_THEMIS", "ID2_LAMP5", "VIP", "VIP","VIP", "ID2_PAX6", "ID2_NCKAP5","PVALB", "PVALB","PVALB","VIP", "VIP", "VIP", "SST", "SST", "SST", "SST", "L5_6_FEZF2_TLE4", "L5_6_FEZF2_TLE4","SST", "L2", "L4_5_FEZF2_LRRK1",  "L4_RORB_SCHLAP1", "L2", "ID2_LAMP5", "ID2_LAMP5", "SST", "L5_FEZF2_ADRA1A", "L5_6_THEMIS", "L5_6_FEZF2_TLE4", "L4_RORB_SCHLAP1", "L4_RORB_SCHLAP1")
fac <- setNames(cdmediumnames,highnames)
cdmediumanno <- setNames(recode_factor(highanno, !!!fac), names(highanno))

caocd <- Cacoa$new(con, ref.level="Ctr", 
                   target.level="Scz", 
                   sample.groups=samplegroups, 
                   cell.groups=cdmediumanno, n.cores=5)

caocd$estimateCellDensity(n.cores = 10,
                          name = "default", 
                          verbose = TRUE)
```

``` r
#get list of plots
pdens <- caocd$plotCellDensity(name = "default", add.points=TRUE, show.grid=TRUE,
                          show.cell.groups=FALSE, 
                          contours = c("PVALB", "SST", "L2", "L4_5_FEZF2_LRRK1", "VIP"), 
                          contour.conf = 98, 
                          show.legend = FALSE, cell.groups = cdmediumanno) 

#Make plot for Ctr
pdens1 <- pdens$Ctr + theme(panel.background = element_rect(fill = "black",
                                colour = "black",
                                size = 0, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.grid = element_blank(),
          panel.border = element_blank(), 
          text = element_text(size = 15)) + 
  ggtitle("control")

#Make plot for Scz
pdens2 <- pdens$Scz + theme(panel.background = element_rect(fill = "black",
                                colour = "black",
                                size = 0, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.grid = element_blank(),
          panel.border = element_blank(), 
          text = element_text(size = 15)) + 
  ggtitle("schizophrenia")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2dens.jpg" width="40%" style="display: block; margin: auto;" />

``` r
library(igraph)
caocd$estimateDiffCellDensity(type = "t.test", verbose=TRUE, n.cores = 5,name = "default", smooth = TRUE)

pcompdif <- caocd$plotDiffCellDensity(type='t.test', name = "default", adjust.pvalues = FALSE,contours =  c("PVALB", "SST", "L2", "L4_5_FEZF2_LRRK1", "VIP"), contour.conf = 98, contour.color = "white", show.legend = FALSE) + scale_colour_gradient2(
    low = "purple",
    mid = "black",
    high = "yellow", breaks = c(-4,0,4)) + theme(panel.background = element_rect(fill = "black",
                                colour = "black",
                                size = 0, linetype = "solid"), 
                                panel.grid = element_blank(),
        panel.border = element_blank(), 
        text = element_text(size = 15)) + 
  ggtitle("compositional difference")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2cd.jpg" width="40%" style="display: block; margin: auto;" />

``` r
#Create disease UMAP embedding from cacoa and add annotation to it

levels(cao$sample.groups) <-c("control", "schizophrenia")
cao$sample.groups.palette <-  c("#980e5c","#435790")
pcao1 <- cao$plotEmbedding(groups = cao$getConditionPerCell(), 
                           plot.na = FALSE, 
                           color.by = "condition", show.legend = TRUE, 
                           legend.position = c(0.1,0.1), 
                           plot.theme = theme_bw(), 
                           font.size = c(4,5))
pcao <- cao$plotEmbedding(groups = highanno, plot.na = FALSE, palette  = rep("white",35)) #for annotation purpose
pcao1$layers[[2]] <- pcao$layers[[2]] #add annotation layer
pcao1 <- pcao1 + guides(color = guide_legend(override.aes = list(size = 7, alpha = 1))) + theme(legend.text =  element_text(size = 15), legend.title = element_blank(), legend.background = element_rect(fill=alpha('white', 0.2)))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2umap.jpg" width="40%" style="display: block; margin: auto;" />

### 4\. Estimate cell loadings (Extended figure 5.C)

``` r
cao$estimateCellLoadings(n.seed = 1, n.cores = 10, coda.test = "significance")
celloadh <- cao$plotCellLoadings(show.pvals = FALSE)
col <- c("blue", "blue", "red", "red", "blue", "red", "red", "red", "red", "blue", "blue","red", "blue","blue","blue", "blue", "red", "red", "blue", "red", "blue", "blue", "red", "blue", "red", "blue", "blue", "red", "blue", "blue", "red", "blue", "blue", "blue", "red") %>% rev #hardcoded color vector for easier distinction between inh vs exc neurons

ggplot(celloadh$data, aes(y = ind, x = values, fill = ind)) + geom_boxplot(outlier.shape = NA, notch = TRUE) + geom_hline(yintercept = 9.5, color = "red") + scale_fill_manual(values = high.pal) + 
  theme_bw() +  
  theme(axis.text.y = element_text(size = 12, color = col),
        legend.position = "None", axis.title.x = element_text(size = 15), 
        plot.margin = unit(c(1,0.2,0.5,0.2), "cm")) + 
  xlab("separating coefficient") + ylab("") + xlim(c(-0.65,0.65)) + 
  geom_segment(aes(x = 0.1, y = 37, xend = 0.65, yend = 37), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  scale_y_discrete(expand = expansion(add = c(0.5, 4))) + 
  geom_segment(aes(x = -0.1, y = 37, xend = -0.65, yend = 37), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 0.35, y=38, 
           label = "higher in Scz", 
           size = 5) + 
  annotate("text", x = -0.35, y=38, 
           label = "higher in Ctr", 
           size = 5) + 
  geom_vline(xintercept = 0, color = "#484848")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig5c.jpg" width="40%" style="display: block; margin: auto;" />

Same for medium anno (Figure 2.B)

``` r
caom <- Cacoa$new(con, ref.level="Ctr", target.level="Scz", sample.groups=samplegroups, cell.groups=tans$med.clean, n.cores=5)
caom$estimateCellLoadings(n.seed = 1, n.cores = 10, coda.test = "significance")

celload <- caom$plotCellLoadings(show.pvals = FALSE)
col2 <- c("blue", "red", "red", "red", "blue", "red", "blue", "red", "blue", "blue", "blue", "red", "red", "red") %>% rev

ggplot(celload$data, aes(y = ind, x = values, fill = ind)) + 
  geom_boxplot(outlier.shape = NA, notch = TRUE) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  scale_fill_manual(values = med.pal) + 
  theme_bw() +  
  theme(axis.text.y = element_text(size = 12, color = col2),
        legend.position = "None", axis.title.x = element_text(size = 15), 
        plot.margin = unit(c(1,0.2,0.5,0.2), "cm")) + 
  xlab("separating coefficient") + 
  ylab("") + 
  xlim(c(-0.8,0.8)) + 
  geom_segment(aes(x = 0.1, y = 16, xend = 0.8, yend = 16), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  scale_y_discrete(expand = expansion(add = c(0.5, 4))) + 
  geom_segment(aes(x = -0.1, y = 16, xend = -0.8, yend = 16), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 0.45, y=17, 
           label = "higher in Scz", 
           size = 5) + 
  annotate("text", x = -0.45, y=17, 
           label = "higher in Ctr", 
           size = 5) + 
  geom_vline(xintercept = 0, color = "#484848")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2b.jpg" width="40%" style="display: block; margin: auto;" />

### 5\. Expression shift magnitued/distance

``` r
cao$estimateExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                      min.samp.per.type=5, min.gene.frac=0.05, dist.type="cross.ref")
caom$estimateExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                      min.samp.per.type=5, min.gene.frac=0.05, dist.type="cross.ref")
cao$plotExpressionShiftMagnitudes(show.jitter = FALSE) + scale_fill_manual(values =med.pal) + ylab("normalized distance")
caom$plotExpressionShiftMagnitudes(show.jitter = FALSE) + scale_fill_manual(values =med.pal) + ylab("normalized distance")
```

function for ploting estimation of layers

``` r
#caoinput is the cacoa object
#layerdtf is the dataframe of medium and high transfered annotations from Allen Brain data
#pal are the palettes

plotLayers <- function(caoinput, layerdtf, pal){
    x <- caoinput
    mv <- tapply(x$data$val,x$data$Type,median)
      df <- layerdtf; df$ndem <- mv[as.character(df$Subtypes)]; df <- na.omit(df);

        p <- ggplot(df, aes(x = reorder(Subtypes,ndem,mean), y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.15, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 15), axis.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
  
  return(p)
  
}
layhigh <- plotLayers(cao$plotExpressionShiftMagnitudes(jitter.alpha = 0.02), ldf_high, high.pal)
laymed <- plotLayers(caom$plotExpressionShiftMagnitudes(jitter.alpha = 0.02), ldf_med, high.pal)
```

``` r
#Figure 2.E
ggpubr::ggarrange(cao$plotExpressionShiftMagnitudes(show.jitter = FALSE) + scale_fill_manual(values =high.pal) + ylab("normalized distance") + theme(axis.text.x=element_blank()),
                  layhigh,ncol=1,nrow=2,heights=c(0.5,1),align = 'v')
```

``` r
#Figure 2.D
ggpubr::ggarrange(caom$plotExpressionShiftMagnitudes(show.jitter = FALSE) + scale_fill_manual(values =med.pal) + ylab("normalized distance") + theme(axis.text.x=element_blank()),
                  laymed,ncol=1,nrow=2,heights=c(0.5,1),align = 'v')
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2d.jpg" width="40%" style="display: block; margin: auto;" />

joint expression shifts

``` r
#Figure 2.F
cao$plotExpressionDistance(joint = TRUE, notch = TRUE, show.significance = TRUE)+ scale_fill_manual(labels = c("Control", "Schizophrenia"), values = c("#980e5c","#435790")) +theme(legend.text = element_text(size = 15), axis.title.y = element_text(size = 15))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2f.jpg" width="40%" style="display: block; margin: auto;" />

MDS plot

``` r
#Figure 2.G
p1 <- cao$plotExpressionDistanceEmbedding(method = "MDS",show.sample.size = FALSE, font.size = 0, size = 5, palette = c("#980e5c","#435790"))
p1$data$condition <- ifelse(p1$data$condition =="Scz", "schizophrenia", "control")
p1 + theme(legend.title = element_text(size = 0), legend.text = element_text(size = 15))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/fig2g.jpg" width="40%" style="display: block; margin: auto;" />

expression distance

``` r
#Extended data figure 5.E
cao$plotExpressionDistance(notch = TRUE, show.significance = FALSE, alpha = 0.1)+ scale_fill_manual(labels = c("Control", "Schizophrenia"), values = c("#980e5c","#435790"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/efig5e.jpg" width="60%" style="display: block; margin: auto;" />

common expression shifts

``` r
cao$estimateCommonExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                            min.samp.per.type=5, min.gene.frac=0.05)
caom$estimateCommonExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                            min.samp.per.type=5, min.gene.frac=0.05)
```

Extended data figure 5.H

``` r
plotLayersCommon <- function(caoinput, layerdtf, pal){
    x <- caoinput
    mv <- tapply(x$data$val,x$data$Type,median)
      df <- layerdtf; df$ndem <- mv[as.character(df$Subtypes)]; df <- na.omit(df);

        p <- ggplot(df, aes(x = reorder(Subtypes,ndem,mean), y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.15, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 15), axis.title = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
  
  return(p)
  
}
```

``` r
c1 <- plotLayersCommon(cao$plotCommonExpressionShiftMagnitudes(type = "box", show.subsampling.variability = FALSE), ldf_high,high.pal)
ggpubr::ggarrange(cao$plotCommonExpressionShiftMagnitudes(type = "box", show.subsampling.variability = FALSE) + ylab("normalized distance (common)") + theme(axis.text.x=element_blank()),
                  c1,ncol=1,nrow=2,heights=c(0.5,1),align = 'v')
```

``` r
c4 <- plotLayersCommon(caom$plotCommonExpressionShiftMagnitudes(type = "box", show.subsampling.variability = FALSE), ldf_med,med.pal)

ggpubr::ggarrange(caom$plotCommonExpressionShiftMagnitudes(type = "box", show.subsampling.variability = FALSE) + ylab("normalized distance (common)") + theme(axis.text.x=element_blank()),
                  c4,ncol=1,nrow=2,heights=c(0.5,1),align = 'v')
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/commonm.jpg" width="40%" style="display: block; margin: auto;" />

### 6\. Creating pseudobulk matrix for plotting proportions of chosen genes per each cell type

``` r
cms.psbulk <- lapply(setNames(names(cms0), names(cms0)), function(smpl) {
  lapply(setNames(levels(as.factor(annotations_scz2_no_MB8$high.clean)),
                  levels(as.factor(annotations_scz2_no_MB8$high.clean))), 
         function(subt) {
          #get nuclei names for sample/subtype combination    
          cells.subt <- split(as.factor(annotations_scz2_no_MB8$high.clean), 
          gsub("(^MB[[:digit:]][[:digit:]]?-?[[:digit:]]?)(.*)", "\\1", 
               names(as.factor(annotations_scz2_no_MB8$high.clean))))
          cells <- cells.subt[[smpl]][cells.subt[[smpl]] == subt] %>% names
          #filter cell/subtype combinations with < 10 cells + get pseudobulk
          if(length(cells) >= 10) {  
          pseudobulk <- rowSums(cms0[[smpl]][, cells])
          #normalize
          pseudobulk * 1000000 /sum(pseudobulk)}}) %>% 
          #remove NULL objects and merge to matrix
          plyr::compact() %>% 
          sapply(function(x) x)})

psbulk.df <- lapply(setNames(names(cms.psbulk), names(cms.psbulk)), 
                    function(smpl) {
                    lapply(setNames(cms.psbulk[[smpl]] %>% colnames, cms.psbulk[[smpl]] %>% colnames), 
                           function(subt) {
                            data.frame(sample = smpl, 
                                       subtype = subt, 
                                       condition = as.character(diseasef[smpl]), 
                                       t(cms.psbulk[[smpl]][, subt]))}) %>%
                        do.call(rbind, .)}) %>% 
  do.call(rbind, .)

psbulk.df$subtype <- factor(psbulk.df$subtype, levels = hsub_order)
psbulk.df$condition <- factor(psbulk.df$condition, levels = c("control", "schizo"))
psbulk.df$condition %<>% gsub("schizo", "SZ", .) %>% factor
```

``` r
#Extended data figure 5.I
ggplot(data = psbulk.df, aes_(y = as.name("CALB2"), x = ~subtype, dodge = ~condition, fill = ~condition)) +
    geom_point(aes(color = condition), 
               position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.5), 
               size = 1.1, 
               alpha = 0.4, 
               show.legend = F) +
  geom_boxplot(notch = F, 
               outlier.shape = NA, 
               alpha = 0.4, 
               lwd = 0.4, 
               width = 0.7) +
  stat_summary(fun = "median", 
               geom = "point",
               size = 2,
               show.legend = F,
               position = position_dodge(width = 0.7),
               mapping = aes(group = condition), 
               fill = "grey",
               shape = 23) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, 
                                   hjust=1, 
                                   vjust=0.5, 
                                   size = 12),
          axis.text.y=element_text(size = 15),
          axis.title.x = element_blank(),
          text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "top",
          plot.subtitle = element_text(size = 15, 
                                       hjust = 0.5, 
                                       face = "bold"),
          panel.grid.minor = element_blank()) +
          labs(x = NULL, y = "Normalized UMIs", subtitle = "CALB2") +
    scale_color_manual(values = palette_45_2[c(8, 15)], 
                       labels = c("control", "schizophrenia")) +
    scale_fill_manual(values = palette_45_2[c(8, 15)], 
                      labels = c("control", "schizophrenia"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/calb2.jpg" width="60%" style="display: block; margin: auto;" />

``` r
null.length.signif <- function(x) {
  if (signif(x)[
        x <= 0.05 & !is.na(x)] %>% length > 0) {
  signif(x)[x <= 0.05 & !is.na(x)]} else {NULL}}

null.length <- function(vect, x) {
  if (x[x <= 0.05 & !is.na(x)] %>% length > 0) {vect[x <= 0.05 & !is.na(x)]} else {NULL}}
```

``` r
#Extended data figure 6.E

inhib.subs <- 
  psbulk.df[grepl("^ID2|^SST|^PVALB|^VIP", psbulk.df$subtype), ]$subtype %>% unique %>% sort %>% droplevels

gn2 <-lapply(c("DRD2", "GABRA1", "GRIK3", "GRIN2A", "GRIN2C", "OXTR"), function(G) {
ps <-
  lapply(setNames(as.vector(inhib.subs), inhib.subs), function(x) {
    de.high.wald[[x]]$res[G, "pvalue"]}) %>% unlist(recursive = F)
ps %<>% p.adjust(method = "BH")
temp.plot <-
ggplot(data = psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ], 
       aes_(y = as.name(G), x = ~subtype, dodge = ~condition, fill = ~condition)) +
  geom_point(aes(color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.5), 
               size = 1.1, alpha = 0.2, show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.4, lwd = 0.4, width = 0.8) +
  stat_summary(fun = "median", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.8), 
                 mapping = aes(group = condition), fill = "grey",
                 shape = 23) + 
   geom_signif(annotations = null.length.signif(ps),
             xmin = null.length(c(0.75:19.75), ps), 
             xmax = null.length(c(1.25:20.25), ps), 
             y_position = null.length(rep(
               max(psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ][[G]]) * 1.07, 20), 
               ps), textsize = 5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_blank(),
          text = element_text(size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top", panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -3.7, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = "UMI counts", subtitle = G) +
    scale_color_manual(values = palette_45_2[c(8, 15)], labels = c("control", "schizophrenia")) +
    scale_fill_manual(values = palette_45_2[c(8, 15)], labels = c("control", "schizophrenia")) +
  
   ylim(-0.2, max(
     psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ][[G]]) * 1.25)})

gn2 <- setNames(gn2, c("DRD2", "GABRA1", "GRIK3", "GRIN2A", "GRIN2C", "OXTR"))
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/mgenes.jpg" width="60%" style="display: block; margin: auto;" />

``` r
#Figure 5.C

pvals <- lapply(setNames(as.vector(inhib.subs), inhib.subs), 
                function(x) {
                  de.high.wald[[x]]$res["CHRFAM7A", "pvalue"]}) %>% 
                unlist(recursive = F)

#correct p-vals for  multiple comparisons
pvals %<>% p.adjust(method = "BH")
```

``` r
#Figure 5.C
CHRFAM7A.plot <-
ggplot(data = psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ], 
       aes(x = subtype, 
           y = CHRFAM7A, 
           dodge = condition, 
           fill = condition)) +
  geom_point(aes(y = CHRFAM7A, 
                 color = condition),
             position = position_jitterdodge(dodge.width = 0.7, 
                                             jitter.width = 0.5), 
             size = 1.1,
             alpha = 0.2, 
             show.legend = F) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 0.4, lwd = 0.4, width = 0.8) +
  stat_summary(fun = "median", geom = "point",
                 size = 2,
                 show.legend = F, position = position_dodge(width = 0.8), 
                 mapping = aes(group = condition), fill = "grey",
                 shape = 23) + 
   geom_signif(annotations = null.length.signif(pvals),
             xmin = null.length(c(0.75:19.75), pvals), 
             xmax = null.length(c(1.25:20.25), pvals), 
             y_position = null.length(rep(max(psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ][["CHRFAM7A"]]) * 1.05, 20), pvals),
             textsize = 5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 11),
        axis.text.y=element_text(size = 13),
        axis.title.x = element_blank(),
        text = element_text(size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top", panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -3.7, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = "Normalized UMIs", subtitle = "CHRFAM7A") +
  scale_color_manual(values = palette_45_2[c(8, 15)]) +
  scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(-0.2, max(psbulk.df[psbulk.df$subtype %in% as.vector(inhib.subs), ][["CHRFAM7A"]]) * 1.25)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/chrfam7a.jpg" width="60%" style="display: block; margin: auto;" />

``` r
fam_expr <- read.csv("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/2ndary/schizo_notebook_2/fam_expr.csv")
fam_prop <- read.csv("/d0/home/mbatiuk/projects/human_schizo/MB6-MB23/2ndary/2ndary/schizo_notebook_2/fam_prop.csv")
fam_expr$condition %<>% as.factor

#average technical replicates
fam_prop %<>% group_by(sample, condition) %>% summarize_at(c(
  "total_CR", "FAM_pos_CR_pos", "CR_pos_only", "prop_FAM_CR"),
  mean)
fam_prop$prop_FAM_CR %<>% as.numeric
fam_prop$condition %<>% as.factor

fam_expr.median <- group_by(fam_expr, sample, condition) %>% summarize_at("expression", median)
fam_expr.median$expression %<>% as.numeric 
fam_expr.log1p <- fam_expr
fam_expr.log1p$expression %<>% log1p
```

Normality Shapiro-Wilk FAM proportions

``` r
fam.prop.list <- split(fam_prop, f = fam_prop$condition)
#shapiro test
shapiro.fam.prop <-
  lapply(fam.prop.list, function(x) shapiro.test(x$prop_FAM_CR))
#p<0.05 - non-normal distribution; p>0.05 - normal distribution
#normal
```

Normality Shapiro-Wilk FAM express PER CELL

``` r
fam.expr.list <- split(fam_expr, f = fam_expr$condition)
#shapiro test
shapiro.fam.expr <-
  lapply(fam.expr.list, function(x) shapiro.test(x$expression))
#p<0.05 - non-normal distribution; p>0.05 - normal distribution
#non normal
```

Normality Shapiro-Wilk FAM express PER SAMPLE MEDIAN

``` r
fam.expr.list2 <- split(fam_expr.median, f = fam_expr.median$condition)
#shapiro test
shapiro.fam.expr2 <-
  lapply(fam.expr.list2, function(x) shapiro.test(x$expression))
#p<0.05 - non-normal distribution; p>0.05 - normal distribution
#normal
```

Equality of variances (Levene’s test) FAM PROPS

``` r
library(car)
levene.fam.prop <- leveneTest(fam_prop$prop_FAM_CR ~ fam_prop$condition)
#p < 0.05 - unequal variances
#equal variances
```

Equality of variances (Levene’s test) FAM EXPR PER CELL

``` r
levene.fam.expr <- leveneTest(fam_expr$expression ~ fam_expr$condition)
#p < 0.05 - unequal variances
#non-equal variances
```

Equality of variances (Levene’s test) FAM EXPR MEDIAN

``` r
levene.fam.expr2 <- leveneTest(fam_expr.median$expression ~ fam_expr.median$condition)
#p < 0.05 - unequal variances
#equal variances
```

T-test - equal variances FAM PROPS

``` r
ttest.fam.prop <- t.test(fam_prop$prop_FAM_CR ~ fam_prop$condition, var.equal = T)
#non-significant
```

Independent 2-group Mann-Whitney U test FAM EXPR PER CELL

``` r
wilcox.fam.expr <- wilcox.test(fam_expr$expression ~ fam_expr$condition)
```

T-test - equal variances FAM EXPRESSION PER SAMPLE

``` r
ttest.fam.expr <- t.test(fam_expr.median$expression ~ fam_expr.median$condition, var.equal = T)
#non-significant difference
```

``` r
fam_expr.median$condition <- ifelse(fam_expr.median$condition == "SZ", "Schizophrenia", "Control")
fam.expr.plot <-
  ggplot(data = fam_expr.median, aes(x = condition, y = expression, fill = condition)) +
  
  
  rasterise(geom_sina(size = 1.4, alpha = 0.3, aes(color = condition), show.legend = T, scale = F), dpi = 600) +
 
  stat_boxplot(geom ='errorbar',  width = 0.4) +
  
  
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1, width = 0.13, color = "black", fill = "white", coef = 0) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 13),
         # axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
         legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
  labs(x = NULL, y = "CHRFAM7A mRNA levels\nin CR+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  
ylim(-0.0001, max(fam_expr.median$expression)*1.3)
fam.expr.plot


#CHRFAM7A mRNA levels in CR+ neurons
```

``` r
fam_prop$condition <- ifelse(fam_prop$condition == "SZ", "Schizophrenia", "Control")

fam.prop.plot <-
ggplot(data = fam_prop, aes(x = condition, y = prop_FAM_CR, fill = condition)) +
  geom_point(aes(y = prop_FAM_CR, color = condition), 
               size = 1.5, alpha = 0.5, show.legend = T, position = position_jitter(width = 0.2)) +
  
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = F, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey") +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(size = 13),
          
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = "CR+ neuron fraction\nexpressing CHRFAM7A") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  ylim(0.000, max(fam_prop$prop_FAM_CR)*1.1)
fam.prop.plot


library(cowplot)
  legend <- get_legend(
  # create some space to the left of the legend
  fam.prop.plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

pl <- plot_grid(plotlist = list(plot_grid(plotlist = list(NULL, legend, NULL), nrow = 1), 
                          plot_grid(plotlist = list(fam.prop.plot + theme(legend.position = "None"), 
                                                    fam.expr.plot + theme(legend.position = "None")), nrow = 1)), 
                          nrow = 2,
          rel_heights = c(0.2,2))

title <- ggdraw() + 
  draw_label("CHRFAM7A in Layer 2",
    fontface = 'bold',
    x = 0,
    hjust = -0.8)


pl <- plot_grid(
  title, pl,
  ncol = 1,
  rel_heights = c(0.1, 1)
); pl

#Figure 5.E
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/chrfam7a2.jpg" width="40%" style="display: block; margin: auto;" />
