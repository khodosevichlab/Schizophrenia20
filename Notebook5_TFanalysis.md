Notebook 6: Transcription factor analysis
================

``` r
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(data.table)
library(colorRamps)
library(circlize)
```

### 1\. Load DE genes calculated in Notebook3

``` r
de <- readRDS("~/RProjects/schizo/sp2/del.fixedW.RDS")
dehigh <- de$high
dehigh <- lapply(dehigh, "[[", "res")
```

### 2\. Create DE object for viper input

``` r
dehm <- lapply(dehigh, function(x){ 
    data.table(x) %>% 
    dplyr::select(Gene, Z) %>% 
    dplyr::filter(!is.na(Z)) %>% 
    column_to_rownames(var = "Gene") %>%
    as.matrix()})
```

### 3\. Extract dorothea regulons

``` r
## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))
```

### 4\. Run dorothea + viper

``` r
tf_activities_stat <- lapply(dehm, function(x) {
  dorothea::run_viper(x, regulons,
    options =  list(minsize = 3, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, nes = TRUE))})
```

### 5\. extract top 10 +NES and -NES TFs for each cell type and make union for HM plot

``` r
tf_activities_stat_top10_up <- lapply(tf_activities_stat, function(x){
  x%>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "Z") %>%
    dplyr::top_n(10, wt = NES) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))})


tf_activities_stat_top10_down <- lapply(tf_activities_stat, function(x){
  x%>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "Z") %>%
    dplyr::top_n(10, wt = -NES) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))})

top10up <- rbindlist(tf_activities_stat_top10_up)$GeneID %>% unique %>% as.character
top10down <- rbindlist(tf_activities_stat_top10_down)$GeneID %>% unique %>% as.character

m <- do.call(plyr::rbind.fill.matrix, lapply(tf_activities_stat, t))
rownames(m) <- names(tf_activities_stat)
m <- t(m)
m <- m[c(top10up, top10down),] %>% unique
```

### 6\. plot big heatmap

``` r
library(circlize)

col_fun = colorRamp2(c(-2,0,2), c("darkblue", "white" ,"red4"))
hm <- Heatmap(scales::rescale(scale(m[rowSums(ifelse(is.na(m) == TRUE, 1,0))<10,]), to = c(-2,2)), clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", na_col = "gray", column_split = 8, row_split = 6, row_title = NULL, column_title = NULL, heatmap_legend_param = list(title = "Score"), col = col_fun)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/tf.jpg" width="60%" style="display: block; margin: auto;" />

### 7\. Plot picked TFs

``` r
col_fun = colorRamp2(c(-8, -2, 1), c("steelblue", "white", "firebrick"))


downgenes <- c("FOXP1", "MYC", "HMBOX1", "ZEB2", "HIF1A", "SREBF1", "SREBF2")
set.seed(1)
hmdown <- Heatmap(m[rownames(m) %in% downgenes,], clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", na_col = "gray", column_km = 6 , row_title = NULL, column_title = NULL, heatmap_legend_param = list(title = "NES"), 
        column_names_max_height = unit(15, "cm"), column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10), col = col_fun)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/negtf.jpg" width="60%" style="display: block; margin: auto;" />

### 8\. Plot picked TFs

``` r
upgenes <- c("TEAD1", "POU5F1", "TCF12", "TCF4", "TEAD4", "PRDM14", "ASCL1")
col_fun2 = colorRamp2(c(0, 5, 10), c("steelblue", "white", "firebrick"))

set.seed(1)
hmup <- Heatmap((m[rownames(m) %in% upgenes,]), clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", na_col = "gray", row_title = NULL, column_title = NULL, heatmap_legend_param = list(title = "NES"), 
        column_names_max_height = unit(15, "cm"), column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10), col = col_fun2, column_km = 3)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/postf.jpg" width="60%" style="display: block; margin: auto;" />
