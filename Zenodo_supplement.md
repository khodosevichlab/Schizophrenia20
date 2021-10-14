Notebook+: some supplement from zenodo
================

Objects in this notebooks were made in previous ones, please refer to
those notebooks for how they were made.

1.  IQR plot

<!-- end list -->

``` r
con <- readRDS("con_integrated.RDS")
high.pal <- readRDS("high.pal.RDS")
sample_groups <- readRDS("samplegroups.RDS")

highanno <- anno$high.clean
```

``` r
testdata <- merge(data.frame("cell" = names(highanno), "cluster" = highanno, row.names = NULL),data.frame("cell" = names(con$getDatasetPerCell()), "sample" = con$getDatasetPerCell(), row.names = NULL))
testdata$group <- paste0(testdata$sample, "_", testdata$cluster)

library(data.table)
setDTthreads(20)
testdata <- data.table(testdata)
testdata$condition <- recode_factor(testdata$sample, !!!samplegroups) 
testdata[, sample(cell, replace = TRUE), by = .(condition, cluster)]
```

``` r
cao <- Cacoa$new(con, ref.level="control", target.level="schizo", 
                       sample.groups=sample.groupsa, 
                       cell.groups=highanno,
                       cell.groups.palette = high.pal, n.cores=1)
```

``` r
cols <- testdata$cluster %>% unique
library(dplyr)
devtools::load_all("~/cacoa/cacoa")
library(data.table)
highclean <- highanno


ll2<- plapply(seq(1:500), function(x){
  ttd <- testdata[sample %in% c(sample(cao$sample.groups[cao$sample.groups == "control"] %>% names, 5), sample(cao$sample.groups[cao$sample.groups == "schizo"] %>% names, 5)),]
   res <- Cacoa$new(con, ref.level="control", target.level="schizo", sample.groups=sample_groups,
                    cell.groups=cell_groups$high.clean[names(setNames(ttd$sample, ttd$cell))],
                    cell.groups.palette = high.pal, n.cores=1)
    
    res$sample.per.cell <- setNames(ttd$sample, ttd$cell)
  res <- res$estimateExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                      min.samp.per.type=5, dist.type="cross.ref", 
                                      min.gene.frac=0.05)
  res <- melt(res$dists.per.type) %>% with(., tapply(X = value, INDEX = L1,FUN = mean))
    return(res[cols])
}, n.cores = 25, mc.preschedule = TRUE, fail.on.error = TRUE) %>% do.call(rbind,.)
```

``` r
dtf_s <- cbind(ll_samples %>% 
                 colMedians(na.rm = TRUE) %>% 
                 rank %>% setNames(colnames((ll_samples))), 
               (ll_samples) %>% 
                 rowRanks() %>%  
                 colQuantiles(na.rm = TRUE) %>% .[,c(2,4)])
colnames(dtf_s) <- c("R" , "L", "U")
library(ggplot2)
dtf_s <- data.frame(ID = rownames(dtf_s), dtf_s)
dtf_s <- dtf_s[order(dtf_s$R),]
dtf_s$ID <- factor(dtf_s$ID, levels = dtf_s$ID)

ggplot(dtf_s[order(dtf_s$R),], aes(y = R, x = ID))+
  geom_errorbar(aes(ymax = U, ymin = L), color = high.pal[-which( names(high.pal) %in% c("Glia", "Other", NA))][levels(dtf_s$ID)])   +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = seq(1:35)) + 
  ylab("Rank") + 
  xlab("")
```

``` r
df <- ldf_high
df$ndem <- setNames(dtf_s$R, dtf_s$ID)[as.character(df$Subtypes)]
df <- na.omit(df)
layhigh <- df

layhigh$Subtypes <- layhigh$Subtypes %>% as.character
layhigh$Subtypes <- layhigh$Subtypes %>% factor(., levels = dtf_s[order(dtf_s$R),]$ID %>% as.character)
levels(layhigh$Subtypes)


 ggplot(layhigh, aes(x = Subtypes, y = Layers,  color = Subtypes)) +
  geom_jitter(width = 0.35, height = 0.45, alpha = 0.4, size = 0.15, show.legend = F) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2, size = 0.6, colour = "grey") +
  scale_color_manual(values = high.pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "Helvetica"), axis.title.x = element_blank())
```

2.  PMI CANCER plot

<!-- end list -->

``` r
library(dplyr)
library(ggplot2)
library(cowplot)

metadata_new <- metadata_new[,c("Identifier","PMI_hrs", "Cancer")]
metadata_new <- metadata_new[metadata_new$Identifier %in% c("MB20", "MB50", "MB52") != TRUE,]

    panel_pmi <- panel_pmi %>% setNames(con$samples%>% names)
    panel_pmi <- panel_pmi[setNames(c(metadata_new$PMI_hrs), c(metadata_new$Identifier)) %>% sort %>% names]
    metadata_new <- data.table::data.table(metadata_new)
    metadata_new <- metadata_new[order(PMI_hrs)]
    
    panel_pmi <- mapply(function(x,y,z){
      rasterize(panel_pmi[[x]], dpi = 300) + 
      geom_text(x=50, y=70, label= paste0("PMI: ", y)) +
      geom_text(x=45, y=50, label= paste0("Cancer: ", z)) + 
        ggtitle(names(panel_pmi[x]))}, 
                        x = c(1:length(panel_pmi)),
                        y = metadata_new$PMI_hrs, 
                        z = metadata_new$Cancer,
                        SIMPLIFY = FALSE)
    
    pmidf <- data.frame(x = rep(1,5), y = c(1,2,3,4,5), group = nm)
    legendpmi <- ggplot(pmidf, 
                        aes(x = x, y = y, group = group, color = group), 
                        show.legend = FALSE) +
      geom_point(size = 3) + 
      theme(legend.title = element_blank())
    leg <- ggdraw(get_legend(legendpmi))
    panel_pmi[[25]] <- leg

panel_pmi_grid <- (plot_grid(plotlist = panel_pmi, ncol = 5))
```

3.  NULL DISTRIBUTION plot

<!-- end list -->

``` r
cao$estimateExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, 
                                      min.samp.per.type=5, dist.type="cross.ref", 
                                      min.gene.frac=0.05)

nd <- cao$test.results$expression.shifts$randomized.dists %>% stack
colnames(nd) <- c("value", "key")
nd2 <- cao$test.results$expression.shifts$dists.per.type %>% stack
nd2 <- cao$test.results$expression.shifts$dists.per.type %>% sapply(function(x) mean(x, trim=0.2)) %>% data.frame(key = names(.), obsdif = .)

p <- ggplot(nd, 
       aes(x = value, 
           group = key)) + 
  geom_density(adjust = 4,
               aes(y = ..scaled..)) + 
  geom_vline(data = nd2, 
             aes(xintercept = obsdif), 
             linetype='dashed', 
             color = high.pal[-which( names(high.pal) %in% c("Glia", "Other", NA))], size = 1) +
  facet_wrap(~key, nrow = 5) + theme_bw() + theme(strip.text = element_text(size=7)) + ylab("Scaled density") + xlab("Shifts") + ggtitle("Calculated expression shift position in relation to normal distribution estimated by permutation test")
```
