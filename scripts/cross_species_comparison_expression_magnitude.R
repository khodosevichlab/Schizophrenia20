library(Seurat)
library(pagoda2)
library(conos)
library(dplyr)
library(ggplot2)
library(ggridges)
library(magrittr)
source("/home/qiwenhu/pip/scripts/active_script/processing.func.R")
#source("/home/qiwenhu/software/speciesTree/R/evol.R")
#source("/home/qiwenhu/schizo/scripts/shift.func.R")
source("/home/qiwenhu/schizo/scripts/cocoa/R/expression_shifts.R")
setwd("/home/qiwenhu/schizo/data/BICCN_data/")

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
