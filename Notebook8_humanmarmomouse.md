Notebook 8: Human marmoset mouse
================

``` r
library(Seurat)
library(pagoda2)
library(conos)
library(dplyr)
library(ggplot2)
library(ggridges)
library(magrittr)
setwd("/home/qiwenhu/schizo/data/BICCN_data")
dir <- "/home/qiwenhu/schizo/data/BICCN_data"
source("~/processing.func.R")
source("~/expression_shifts.R")
```

### 1\. Load datasets, integrate with conos

``` r
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
```

### 2\. Load annotation file, estimate expression shifts using ‘expression\_shifts.R’

``` r
human.annot <- readRDS("/home/qiwenhu/schizo/data/BICCN_data/human.annot.integrate.rds")
mouse.annot <- readRDS("/home/qiwenhu/schizo/data/BICCN_data/mouse.annot.integrate.rds")
marmo.annot <- readRDS("/home/qiwenhu/schizo/data/BICCN_data/marmoset.annot.integrate.rds")
annot <- c(as.character(human.annot), as.character(mouse.annot), as.character(marmo.annot))
names(annot) <- c(names(human.annot), names(mouse.annot), names(marmo.annot))
integrated.con <- readRDS("/home/qiwenhu/schizo/data/BICCN_data/human.marmo.mouse.con.rds")
samplegroups <- list(
  nonhuman = c('mouse', "marmoset"),
  human = c('human')
)
names(integrated.con$samples) <- c("human", "mouse", "marmoset")
sgf <- setNames(rep(names(samplegroups),unlist(lapply(samplegroups,length))),unlist(samplegroups))

mag.shift <- estimateExpressionShiftMagnitudes(lapply(integrated.con$samples, conos:::getRawCountMatrix), 
                                               sample.groups=sgf, cell.groups=annot, dist='cor',n.cells=500, 
                                               n.subsamples=20, min.cells=10, n.cores=10, within.group.normalization=FALSE)
saveRDS(mag.shift, "exp.shift.cor.rds")
```

### 3\. Expression shifts: human - mouse

``` r
# human mouse
human.mouse.group <- list(
  mouse = c("mouse"),
  human = c('human')
)
names(integrated.con$samples) <- c("human", "mouse", "marmoset")
sgf <- setNames(rep(names(human.mouse.group),unlist(lapply(human.mouse.group, length))),unlist(human.mouse.group))
human.mouse.mag.shift <- estimateExpressionShiftMagnitudes(lapply(integrated.con$samples, conos:::getRawCountMatrix), 
                                               sample.groups=sgf, cell.groups=annot, dist='cor',n.cells=500, 
                                               n.subsamples=20, min.cells=10, n.cores=10, within.group.normalization=FALSE)
saveRDS(human.mouse.mag.shift, "human.mouse.mag.shift.rds")
```

### 4\. Expression shifts: human - marmo

``` r
#human marmo
human.marmo.group <- list(
  marmoset = c("marmoset"),
  human = c('human')
)
names(integrated.con$samples) <- c("human", "mouse", "marmoset")
sgf <- setNames(rep(names(human.marmo.group),unlist(lapply(human.marmo.group, length))),unlist(human.marmo.group))
human.marmo.mag.shift <- estimateExpressionShiftMagnitudes(lapply(integrated.con$samples, conos:::getRawCountMatrix), 
                                                           sample.groups=sgf, cell.groups=annot, dist='cor',n.cells=500, 
                                                           n.subsamples=20, min.cells=10, n.cores=10, within.group.normalization=FALSE)
saveRDS(human.marmo.mag.shift, "human.marmoset.mag.shift.rds")
```

### 5\. create expr shift dtfs, plot them

``` r
human.marmo.df <- human.marmo.mag.shift$df
human.mouse.df <- human.mouse.mag.shift$df
#df$Type <- as.factor(df$Type)
human.marmo.df.mean <- human.marmo.df %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(mean = mean(value, na.rm = TRUE))
human.marmo.df.sd <- human.mouse.df %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(sd = sd(value, na.rm = TRUE))

human.mouse.df.mean <- human.mouse.df %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(mean = mean(value, na.rm = TRUE))
human.mouse.df.sd <- human.mouse.df %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(sd = sd(value, na.rm = TRUE))

# extract only inter-neurons
interneuron.type <- c("Chandelier", "Lamp5", "Pvalb", "Sncg", "Sst", "Vip")
human.marmo.df.mean <- human.marmo.df.mean[gsub("_.*", "", human.marmo.df.mean$Type) %in% interneuron.type, ]
human.marmo.df.sd <- human.marmo.df.sd[gsub("_.*", "", human.marmo.df.sd$Type) %in% interneuron.type, ]
human.mouse.df.mean <- human.mouse.df.mean[gsub("_.*", "", human.mouse.df.mean$Type) %in% interneuron.type, ]
human.mouse.df.sd <- human.mouse.df.sd[gsub("_.*", "", human.mouse.df.sd$Type) %in% interneuron.type, ]


p1 <- ggplot(human.marmo.df.mean, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#E69F00") +
  geom_errorbar(aes(ymin=mean-human.marmo.df.sd$sd, ymax=mean+human.marmo.df.sd$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") + 
  ggtitle(" ") +  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) + ggtitle("human vs marmoset")

p2 <- ggplot(human.mouse.df.mean, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#56B4E9") +
  geom_errorbar(aes(ymin=mean-human.mouse.df.sd$sd, ymax=mean+human.mouse.df.sd$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") + 
  ggtitle(" ") +  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.7)) + ggtitle("human vs mouse")



## major subtypes
interneuron.type <- c("Chandelier", "Lamp5", "Pvalb", "Sncg", "Sst", "Vip")
human.marmo.df.mean <- human.marmo.df.mean[gsub("_.*", "", human.marmo.df.mean$Type) %in% interneuron.type, ]
human.marmo.df.mean$Type <- gsub("_.*", "", human.marmo.df.mean$Type)

human.marmo.df.mean.major <- human.marmo.df.mean %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(mean = mean(mean, na.rm = TRUE))
human.marmo.df.sd.major <- human.marmo.df.mean %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(sd = sd(mean, na.rm = TRUE))
human.marmo.df.sd.major[1, 2] <- 0

human.mouse.df.mean <- human.mouse.df.mean[gsub("_.*", "", human.mouse.df.mean$Type) %in% interneuron.type, ]
human.mouse.df.mean$Type <- gsub("_.*", "", human.mouse.df.mean$Type)

human.mouse.df.mean.major <- human.mouse.df.mean %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(mean = mean(mean, na.rm = TRUE))
human.mouse.df.sd.major <- human.mouse.df.mean %>% dplyr::group_by(Type) %>% 
    dplyr::summarise(sd = sd(mean, na.rm = TRUE))
human.mouse.df.sd.major[1, 2] <- 0

p1 <- ggplot(human.marmo.df.mean.major, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#E69F00") + ggtitle(" ") +  theme_bw() +
  geom_errorbar(aes(ymin=mean-human.marmo.df.sd.major$sd, ymax=mean+human.marmo.df.sd.major$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) + ggtitle("human vs marmoset")

p2 <- ggplot(human.mouse.df.mean.major, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#56B4E9") +
  geom_errorbar(aes(ymin=mean-human.mouse.df.sd.major$sd, ymax=mean+human.mouse.df.sd.major$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") + 
  ggtitle(" ") +  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.7)) + ggtitle("human vs mouse")
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/humanmarmo1.jpg" width="60%" style="display: block; margin: auto;" />

<img src="C:/Users/Katarina/Desktop/scznotebooks/humanmouse1.jpg" width="60%" style="display: block; margin: auto;" />

\#with scz labels

``` r
new.label.info <- readRDS("/home/qiwenhu/schizo/data/prop.label.rds")
scores <- new.label.info$label.distribution
schizo.scores <- scores[grep("MB", rownames(scores)), ]
label.cleaned <- CellAnnotatoR:::expandAnnotationToClusters(schizo.scores, tans$med.clean)
mapping <- read.table("/home/qiwenhu/schizo/shared_figs/annotation.all.txt", sep="\t", header=TRUE)
cleaned.label.df <- data.frame(cells=names(label.cleaned), annot=label.cleaned)
cleaned.label.df <- merge(cleaned.label.df, mapping[, c(1, 4)], by=c("cells"))


schizo.human.marmo.df.mean <- human.marmo.df.mean
schizo.human.marmo.df.mean$Type <- cleaned.label.df[match(schizo.human.marmo.df.mean$Type, cleaned.label.df$annot), ]$schizo_annot
schizo.human.marmo.df.mean <- schizo.human.marmo.df.mean[!is.na(schizo.human.marmo.df.mean$Type), ]
schizo.human.marmo.df.sd <- human.marmo.df.sd
schizo.human.marmo.df.sd$Type <- cleaned.label.df[match(schizo.human.marmo.df.sd$Type, cleaned.label.df$annot), ]$schizo_annot
schizo.human.marmo.df.sd <- schizo.human.marmo.df.sd[!is.na(schizo.human.marmo.df.sd$Type), ]

schizo.human.mouse.df.mean <- human.mouse.df.mean
schizo.human.mouse.df.mean$Type <- cleaned.label.df[match(schizo.human.mouse.df.mean$Type, cleaned.label.df$annot), ]$schizo_annot
schizo.human.mouse.df.mean <- schizo.human.mouse.df.mean[!is.na(schizo.human.mouse.df.mean$Type), ]
schizo.human.mouse.df.sd <- human.mouse.df.sd
schizo.human.mouse.df.sd$Type <- cleaned.label.df[match(schizo.human.mouse.df.sd$Type, cleaned.label.df$annot), ]$schizo_annot
schizo.human.mouse.df.sd <- schizo.human.mouse.df.sd[!is.na(schizo.human.mouse.df.sd$Type), ]

p1 <- ggplot(schizo.human.marmo.df.mean, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#E69F00") +
  geom_errorbar(aes(ymin=mean-schizo.human.marmo.df.sd$sd, ymax=mean+schizo.human.marmo.df.sd$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") + 
  ggtitle(" ") +  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) + ggtitle("human vs marmoset")

p2 <- ggplot(schizo.human.mouse.df.mean, aes(x=reorder(Type, mean), y=mean)) + 
   geom_bar(stat="identity", position=position_dodge(), alpha=0.6, fill="#56B4E9") +
  geom_errorbar(aes(ymin=mean-schizo.human.mouse.df.sd$sd, ymax=mean+schizo.human.mouse.df.sd$sd), width=.2,
                 position=position_dodge(.9)) + xlab(" ") + ylab("normalized distance") + 
  ggtitle(" ") +  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.7)) + ggtitle("human vs mouse")

cowplot::plot_grid(plotlist=list(p1, p2), ncol=1)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/inhhumanmarmomouse.jpg" width="40%" style="display: block; margin: auto;" />
