#Please refer to cacoa repository for full code 


#' Get loadings by different methods
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param criteria Method to get loadings
#' @param ref.cell.type Reference cell type
#' @return Updated data frame with Z scores
getLoadings <- function(cnts, groups, criteria = 'lda', ref.cell.type = NULL) {
  discriminant.methods <- c('lda', 'svm', 'cda', 'cda.std')
  if(!(criteria %in% discriminant.methods)) stop(paste('The discriminant method', criteria, 'is not supported'))
  if(!is.null(ref.cell.type) && !(ref.cell.type %in% colnames(cnts)))
    stop(paste('Reference cell type', ref.cell.type, 'is not correct. Correct cell types are:',
               paste0(colnames(cnts), collapse = ', ') ))
  #Get freqs
  cnts[cnts == 0] <- 0.1
  # cnts[cnts == 0] <- 1
  freqs <- cnts/rowSums(cnts)
  
  # Get ilr
  psi <- coda.base::ilr_basis(ncol(freqs), type = "default")
  rownames(psi) <- colnames(cnts)
  b <- log(freqs) %*% psi
  
  # PCA
  b.norm <-  apply(b, 2, function(y) y - mean(y))
  pca.res <- prcomp(b.norm)
  pca.loadings <- psi %*% pca.res$rotation
  
  psi <- psi %*% pca.res$rotation
  b <- log(freqs) %*% psi
  b <- apply(b, 2, scale)
  
  
  idx.na = which(colSums(is.nan(b)) == 0)
  b <- b[,idx.na]
  psi <- psi[,idx.na]
  
  # Create dataframe
  b.df <- data.frame(b)
  b.df$groups <- 1*groups
  
  # Optimization
  if(criteria == 'lda') {
    if(is.null(ref.cell.type)) {  # pure linear regression which is proportional LDA when number of classes = 2
      res.lm <- lm(groups ~ ., data = b.df)
      w <- res.lm$coefficients[-1]
      w[is.na(w)] <- 0
      
      res.sum = summary(res.lm)
      # print(res.sum)
      # w[(res.sum$coefficients[-1,4] > 0.01) & (res.sum$coefficients[-1,4] != min(res.sum$coefficients[-1,4]))] = 0
      # w[(res.sum$coefficients[-1,4] != min(res.sum$coefficients[-1,4]))] = 0
      
      # # Compare with lda - proportional!
      # res.lda <- lda(groups ~ ., data = b.df)
      # w / res.lda$scaling
      
    } else {  # Quadratic optimization with the linear condition corresponding to the reference level
      X <- b
      Y <- groups * 1
      X2 <- cbind(X, 1)
      Rinv <- solve(chol(t(X2) %*% X2))
      dvec <- t(Y) %*% X2
      Amat <- t(psi[ref.cell.type,,drop=FALSE])
      Amat <- rbind(Amat, 0)
      res.qp <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)
      
      n.qp <- ncol(X2)
      w <- res.qp$solution[-n.qp]
      
      # # Compare with lm
      # res.lm <- lm(groups ~ ., data = b.df)
      # res.qp$unconstrained.solution[-n.qp] / res.lm$coefficients[-1]
    }
  } else if(criteria == 'svm') {
    # ---- SVM
    b.model <- e1071::svm(groups~., data = b.df, kernel = "linear", scale = FALSE)
    # Get hyperplane parameters
    w <- t(b.model$SV) %*% b.model$coefs  # slope
    # b.svm <- -b.model$rho # intercept
    # # create new score
    # v <- b %*% w + b.svm
  } else if(criteria == 'cda') {  # Canonical discriminant analysis
    model <- lm(b ~ groups)
    cda <- candisc::candisc(model, ndim=1)
    # w <- cda$structure
    w <- cda$coeffs.raw
  } else if(criteria == 'cda.std') {
    b.norm <-  apply(b, 2, function(y) y - mean(y))
    
    # PCA
    pca.res <- prcomp(b.norm)
    pca.loadings <- psi %*% pca.res$rotation
    
    # ANNA: please remain the following commented code,
    # in case someone wants to understand what is going on, and what was attempted
    # # CDA
    # df.pca <- as.data.frame(pca.res$x)
    #
    # model <- lm(pca.res$x ~ groups)
    # cda <- candisc::candisc(model, ndim=1)
    #
    # w <- pca.res$rotation  %*% as.matrix(cda$structure)  # standardized regression coefficients - just correlations
    # # w <- pca.res$rotation  %*% as.matrix(cda$coeffs.raw)  # pure regression coefficients
    
    # # Properties
    # # 1
    # b.pca = pca.res$x
    # scores <- cda$scores$Can1
    # cor(groups, b.pca) / cor(scores, b.pca)  # the same number
    # # 1
    # t(cda$structure) / cor(groups, b.pca)  # the same number
    # # What we can do:
    # w.tmp <- pca.res$rotation %*% t(cor(groups, b.pca))
    # w.tmp / w  # the same number
    
    # New version without CDA
    # To makes things comparable between bootstrap iterations,
    # we need standardized regression coefficients,
    # which are just correlations, when regressors are independent, as after PCA.
    # Standardization is needed,
    # because a regression coefficient reflects the variance of the corresponding regressor,
    # and it significantly depends on bootstrap subsample.
    w <- pca.res$rotation %*% t(cor(groups, pca.res$x))
    
  }
  
  w <- w / sqrt(sum(w ^ 2))
  
  scores <- b %*% w
  if(mean(scores[groups]) < mean(scores[!groups])) w = -w
  # scores <- b %*% w
  # print(c(mean(scores[groups]), mean(scores[!groups])))
  
  loadings <- psi %*% w
  names(loadings) <- colnames(cnts)
  
  return(loadings)
}

#' This function produces the resampling dataset
#' @param cnts Counts of cell types in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param n.iter Number of permutations
#' @param remain.groups TRUE - if the labels of groups remain, FALSE - if to produce the null distribution
#' @param replace.samples TRUE - is bootstrap on samples, FALAE - if to remain samples
#' @param seed random seed
#' @return Updated data frame with Z scores
produceResampling <- function(cnts, groups, n.perm = 1000, remain.groups = TRUE, replace.samples = TRUE, seed = 239) {
  
  if(!(remain.groups %in% c(NULL, TRUE, FALSE))) stop('Parameter remain.groups should be in {NULL, TRUE, FALSE}')
  
  cnts.perm <- list()
  groups.perm <- list()
  set.seed(239)
  for(i in 1:n.perm) {
    # Bootstrap samples
    
    samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = replace.samples)
    if (remain.groups) {  # Remain labels
      groups.tmp <- groups[samples.tmp]
    } else {  # Null distribution
      groups.tmp <- sample(groups, length(groups), replace = replace.samples)
      names(groups.tmp) <- samples.tmp
    }
    # Check that both groups are presented
    if((sum(groups.tmp) == 0) || (sum(!groups.tmp) == 0)) next
    
    # Bootstrap cell types
    cnts.resampling <- apply(cnts[samples.tmp,], 1, function(v) {
      table( c(names(v), sample(names(v), size = sum(v), replace = TRUE, prob = v / sum(v))) ) - 1   })
    
    cnts.tmp <- t(cnts.resampling)
    
    cnts.perm[[length(cnts.perm) + 1]] <- cnts.tmp
    groups.perm[[length(groups.perm) + 1]] <- groups.tmp
    
  }
  
  return(list(cnts = cnts.perm, groups = groups.perm))
}

runCoda <- function(cnts, groups, n.seed=239, n.boot=1000, ref.cell.type=NULL) {
  # Apply bootstrap
  loadings <- lapply(1:n.boot, function(ib) {
    # Create samples by bootstrap
    set.seed(n.seed+ib)
    samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = TRUE)
    groups.tmp <- groups[samples.tmp]
    
    
    # Check that both groups are presented
    while ((sum(groups.tmp) == 0) || (sum(!groups.tmp) == 0)) {
      # Create samples by bootstrap
      samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = TRUE)
      groups.tmp <- groups[samples.tmp]
    }
    cnts.tmp <- cnts[samples.tmp,]
    
    init.tmp <- produceResampling(cnts = cnts.tmp, groups = groups.tmp, n.perm = 1,
                                  replace.samples = FALSE,
                                  remain.groups = TRUE, seed = n.seed+ib)
    
    res <- getLoadings(init.tmp$cnts[[1]], init.tmp$groups[[1]])
    res[colnames(cnts), , drop=F]
  })
  
  # Calculate p-values
  loadings.init <- do.call(cbind, loadings)
  
  # ld <- loadings.init
  # threshold <- rowMeans(ld)[abs(rowMeans(ld)) == min(abs(rowMeans(ld)))]
  # threshold <- 0
  
  # # Calculate p-values of confidence interval by bootstrap
  # if(is.null(ref.cell.type)){
  #   tmp <- sapply(1:nrow(loadings.init), function(i) sum(threshold > loadings.init[i,])) / ncol(loadings.init)
  #   pval <- apply((cbind(tmp, 1-tmp)), 1, min) * 2
  # } else {
  #   ref.level = abs(mean(loadings.init[ref.cell.type,]))
  #   tmp1 <- sapply(1:nrow(loadings.init), function(i) sum(loadings.init[i,] > ref.level)) / ncol(loadings.init)
  #   tmp2 <- sapply(1:nrow(loadings.init), function(i) sum(loadings.init[i,] < -ref.level)) / ncol(loadings.init)
  #   pval = (1 -apply((cbind(tmp1, tmp2)), 1, max)) * 2
  # }
  
  # Calculate p-values of confidence interval by bootstrap
  
  
  tmp <- referenceSet(cnts, groups, p.thresh = 0.1)
  cell.list <- tmp$cell.list
  
  sdt.list <- c()  # standard deviation in a list
  mean.list <- c()  # mean value of loadings in a list
  min.list <- c()  # min values of mean loadings
  n.list <- c()  # number of cell types in a list
  for (i.list in 1:length(cell.list)) {
    loadings.list <- loadings.init[cell.list[[i.list]],]
    sdt.list <- c(sdt.list, sd(c(loadings.list)))
    mean.list <- c(mean.list, mean(c(loadings.list)))
    
    
    n.list <- c(n.list, length(cell.list[[i.list]]))
    if (length(cell.list[[i.list]]) == 1) {  # if a list contains only one sample - it cannot be a reference group
      mean.list[i.list] <- 10
      min.list <- c(min.list, 10)
    } else {
      min.list <- c(min.list, min(abs(apply(loadings.list, 1, mean)))  )
      # min.list <- c(min.list, min(abs(apply(loadings.list, 1, mean))) / n.list[i.list] )
    }
  }
  # print(sdt.list)
  # print(mean.list)
  # print(cell.list)
  
  # Define a cluster with reference cell type
  if (is.null(ref.cell.type)) {
    # id.ref.cluster <- which(sdt.list == max(sdt.list))
    
    # id.ref.cluster <- which(n.list == max(n.list))
    # id.ref.cluster <- id.ref.cluster[abs(mean.list[id.ref.cluster]) == min(abs(mean.list[id.ref.cluster]))]
    id.ref.cluster <- which(abs(mean.list) == min(abs(mean.list)))  # <- working version
    # id.ref.cluster <- which(min.list == min(min.list))  # 
  } else {
    id.ref.cluster <- -1
    for (i.list in 1:length(cell.list)) {
      if (ref.cell.type[1] %in% cell.list[[i.list]]) {
        id.ref.cluster <- i.list
      }
    }
    if (id.ref.cluster == -1) stop('Wrong name of reference cell type')
  }
  ref.cell.type <- cell.list[[id.ref.cluster]]
  
  # sorting of cell types
  cell.types.order <- c()
  for(i.list in order(-abs(mean.list - mean.list[id.ref.cluster]))){
    cell.types.tmp <- cell.list[[i.list]]
    cell.types.tmp <- cell.types.tmp[order(-abs(rowMeans(loadings.init[cell.types.tmp,, drop=FALSE]) - mean.list[id.ref.cluster])  ) ]
    cell.types.order <- c(cell.types.order, cell.types.tmp)
  }
  loadings.init <- loadings.init[cell.types.order,]
  
  # pval <- pvalInLoadingsOrder(cell.types.order, cnts, groups)
  
  ref.load.level <- mean(loadings.init[ref.cell.type,])
  tmp <- rowSums(ref.load.level > loadings.init)
  pval <- (pmin(tmp, ncol(loadings.init) - tmp) + 1) / (ncol(loadings.init) + 1) * 2
  
  names(pval) <- rownames(loadings.init)
  #
  # # ----------------------------
  # # Additional correction of p-values
  # ld <- loadings.init
  # # ld.means <- rowMeans(ld)
  #
  # idx <- order(abs(rowMeans(ld)))
  # ld <- ld[idx,]
  # pval <- pval[rownames(ld)]
  #
  # idx <- names(pval)[order(-pval)]
  # ld <- ld[idx,]
  # ld.means <- rowMeans(ld)
  # pvals_tmp <- pval[idx]
  # pvals_tmp[1] <- 1
  # current.mean <- ld.means[1]
  #
  # for(i in 2:length(ld.means)){
  #   # Fid the index of previous closest
  #   d.prev <- abs(ld.means[1:(i-1)] - ld.means[i])
  #   j <- which(d.prev == min(d.prev))
  #
  #   tmp <- sum(ld[i,] > ld.means[j]) / ncol(ld)
  #   tmp <- min(tmp, 1-tmp)
  #   pvals_tmp[i] <- max(min(tmp, pvals_tmp[j]), pvals_tmp[i])
  # }
  # names(pvals_tmp) <- rownames(ld)
  # # print(pvals_tmp)
  # # Combining p-values
  # p.names <- names(pval)
  # pval <- pmax(pval, pvals_tmp[p.names])
  # names(pval) <- p.names
  # #s----------------------------
  
  padj <- p.adjust(pval, method = "fdr")
  p.names <- names(pval)
  loadings.init <- loadings.init[p.names,]
  return(list(padj=padj,
              loadings.init=loadings.init,
              pval=pval,
              ref.load.level=ref.load.level,
              cell.list=cell.list))
}

referenceSet <- function(freqs, groups, p.thresh=0.05) {
  checkPackageInstalled("psych", cran=TRUE)
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)
  
  for (it in 1:(length(cell.list) - 2)) {
    
    mx <- matrix(0, nrow=length(cell.list), ncol=length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next
        
        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))
        # print('---')
        # print(freqs)
        # print(ratio)
        # print(groups)
        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
        # mx[j, i] <- pval
      }
    }
    if(max(max(mx)) < p.thresh) break
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    # print(c(i, j))
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]
    
    if(it == 1){
      mx.first <- mx
      rownames(mx.first) <- cell.types
      colnames(mx.first) <- cell.types
    }
  }
  len.cell.list <- sapply(cell.list, length)
  id.ref <- which(len.cell.list == max(len.cell.list))
  
  # Make mx.first symmetric
  for(i in 1:nrow(mx.first)){
    for(j in 1:ncol(mx.first)){
      if (j <= i) next
      mx.first[j, i] <- mx.first[i, j]
    }
  }
  
  return(list(ref.set = cell.list[[id.ref[1]]],
              cell.list = cell.list,
              mx.first = mx.first,
              mx.final = mx))
}



sbpDiff <- function(freqs, groups){
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)
  
  sbp <- matrix(nrow = 0, ncol = length(cell.types), dimnames = list(c(), cell.types))
  for(it in 1:(length(cell.list) - 1)){
    
    mx <- matrix(0, nrow = length(cell.list), ncol = length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next
        
        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))
        
        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
        # mx[j, i] <- pval
      }
    }
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    
    sbp.tmp <- rep(0, length(cell.types))
    sbp.tmp[which(cell.types %in% cell.list[[i]])] <- 1
    sbp.tmp[which(cell.types %in% cell.list[[j]])] <- -1
    sbp <- rbind(sbp.tmp, sbp)
    # print(c(i, j))
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]
    
  }
  
  return(t(sbp))
}

pvalInLoadingsOrder <- function(cell.types.order, cnts, groups){
  p <- c()
  cnts[cnts <= 0] <- 0.5
  freqs <- (cnts)/rowSums(cnts)
  for(i in 0:(ncol(freqs) - 2)){
    
    if(i == 0){
      freqs1 <- freqs
    } else {
      freqs1 <- freqs[, -which(colnames(freqs) %in% cell.types.order[1:i] )]
    }
    
    psi <- coda.base::ilr_basis(ncol(freqs1), type = "default")
    rownames(psi) <- colnames(freqs1)
    b <- log(freqs1) %*% psi
    
    x <- summary(lm(groups ~ b))
    pval.lm <- pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
    p <- c(p, pval.lm)
    
  }
  p <- c(p, 1)
  names(p) <- cell.types.order
  return(p)
}


estimateCdaSpaceNew <- function(cnts, groups){
  # ----- Get data -----
  
  cnts[cnts == 0] = 0.5
  freqs <- cnts / rowSums(cnts)
  psi <- coda.base::ilr_basis(ncol(freqs), type = "default")
  rownames(psi) <- colnames(cnts)
  b <- log(freqs) %*% psi
  
  f1 = apply(exp(b %*% t(psi)),1, function(x) {x/sum(x)})
  
  # PCA - to get new ilr basis of principal components
  b.norm <-  apply(b, 2, function(y) y - mean(y))
  pca.res <- prcomp(b.norm)
  pca.loadings <- psi %*% pca.res$rotation
  
  psi <- psi %*% pca.res$rotation  # new Psi matrix
  
  b <- log(freqs) %*% psi  # New ilr coordinates
  b.init <- log(freqs.init) %*% psi
  
  
  
  # Standardization
  b.mean <- colMeans(b)
  b.sd <- apply(b, 2, sd)
  b.init <- (b.init - b.mean) / b.sd
  b <- apply(b, 2, scale)
  
  b.plot <- rbind(b, b.init)
  
  # ------------------------------
  # Find the first projection
  
  # Create dataframe
  b.df <- data.frame(b)
  b.df$groups <- 1*groups
  
  # Linear model
  res.lm <- lm(groups ~ ., data = b.df)
  w <- res.lm$coefficients[-1]  # coefficients without intercept
  
  w <- as.matrix(w / sqrt(sum(w ^ 2)))  # not required
  
  
  # Projections of b.plot on w - is the first coordinate
  b.proj = as.matrix(apply(b.plot, 1, function(x) sum(x*w)))
  b.proj.vec = b.proj %*% t(w)
  b.ort <- b.plot - b.proj.vec
  
  b.ort.proj = as.matrix(apply(b.ort, 1, function(x) sum(x*w))) # <- check that all are small
  if(sum(b.ort.proj > 10^(-10)) != 0) stop('Something is going wrong')
  
  # ------------------------------
  # Find the second projection - from PCA,
  # because the all differences between groups are in the first coordinate
  
  pca.ort <- prcomp(b.ort)
  
  
  b.plot <- data.frame(cbind(b.proj, pca.ort$x[,1]))
  colnames(b.plot) <- c('S1', 'S2')
  # b.plot$group <- c(as.character(uniq.combo$site), 'init')
  # colnames(b.plot)
  
  return(b.plot)
}


estimateCdaSpace <- function(d.counts, d.groups, thresh.pc.var = 0.95, n.dim = 2){
  cell.loadings <- c()
  sample.pos <- c()
  
  bal <- getRndBalances(d.counts)
  d.used <- bal$norm
  
  if(ncol(d.used) != 2){
    for(i in 1:n.dim){
      res.remove <- removeGroupEffect(d.used, d.groups, thresh.pc.var = 0.9)
      cell.loadings <- cbind(cell.loadings, bal$psi %*% res.remove$rotation)
      sample.pos <- cbind(sample.pos, res.remove$scores)
      d.used <- d.used - res.remove$used.part
    }
  } else {
    cell.loadings <- bal$psi
    sample.pos <- d.used
  }
  
  cn <- paste('S', 1:n.dim, sep = '')
  df.cda <- as.data.frame(sample.pos) %>% set_colnames(cn)
  df.loadings <- as.data.frame(cell.loadings * 8) %>% set_colnames(cn)
  return(list(red=df.cda, loadings=df.loadings))
}

plotCodaSpaceInner <- function(df.space, df.loadings, d.groups, ref.level, target.level, font.size=3, palette=NULL) {
  group.names <- c(ref.level, target.level)
  rda.plot <- ggplot(df.space, aes(x=S1, y=S2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Condition")
  
  if(!is.null(palette)) rda.plot <- rda.plot + scale_color_manual(values=palette)
  
  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=S1, y=0, yend=S2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    geom_text(data=df.loadings,
              aes(x=S1, y=S2, label=rownames(df.loadings)),
              color="black", size=font.size)
  
  dx <- max(diff(range(df.space$S1)), diff(range(df.loadings$S1)))
  dy <- max(diff(range(df.space$S2)), diff(range(df.loadings$S2)))
  
  rda.biplot <- rda.biplot + coord_fixed(ratio=dy/dx)
  
  return(rda.biplot)
}

# helper function for creating dendograms
ggdend <- function(dend.data, angle=90, plot.theme=theme_get(), font.size=3) {
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + plot.theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.border=element_blank(),
          axis.line=element_blank()) +
    geom_text(data = dend.data$labels, aes(x, y, label = label),
              hjust=1, angle=angle, size=font.size) + ylim(-0.5, NA)
}

plotContrastTree <- function(d.counts, d.groups, ref.level, target.level, plot.theme,
                             p.threshold=0.05, adjust.pvalues=TRUE, h.methods='both', font.size=3) {
  log.f <- getLogFreq(d.counts)
  
  if(h.methods == 'up'){
    print('up')
    t.cur <- constructTreeUp(d.counts, d.groups)
  } else if (h.methods == 'down'){
    print('down')
    t.cur <- constructTree(d.counts, d.groups)
  } else {
    print('up and down')
    t.cur <- constructTreeUpDown(d.counts, d.groups)
  }
  
  
  # t.cur <- constructBestPartitionTree(d.counts, d.groups)
  
  tree <- t.cur$tree
  sbp <- sbpInNodes(tree)
  # sbp = t.cur$sbp
  
  # ---------------------------------
  # Positions of the dendrogram
  dend.data <- ggdendro::dendro_data(t.cur$dendro, type = "rectangle")
  types.order <- dend.data$labels$label
  # change the sign of sbp corresponding to the dendrodgram
  for (k in 1:ncol(sbp)) {
    p <- sbp[, k]
    type.plus <- rownames(sbp)[p > 0]
    type.minus <- rownames(sbp)[p < 0]
    if(which(types.order == type.plus[1]) < which(types.order == type.minus[1])){
      sbp[, k] <- -sbp[, k]
    }
  }
  
  # Get positions of segments (nodes) of the dendrogram
  # only inner nodes
  node.pos <- dend.data$segments %$% .[(y == yend) & (yend != 0),]
  node.pos$id <- tree$edge[,1]  # id of the inner node
  node.pos$to <- tree$edge[,2]
  
  # Positions of inner nodes
  innode.pos <- unique(node.pos[,c('x','y','id')])
  rownames(innode.pos) <- innode.pos$id
  # Ranges of inner nodes (length of branches)
  innode.pos$range <- -1
  for (i in 1:nrow(innode.pos)){
    tmp <- node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] <- max(tmp) - min(tmp)
  }
  innode.pos = innode.pos[order(innode.pos$id),]
  
  # ----------------------------------------
  
  # Balances
  balances <- getNodeBalances(log.f, sbp)
  colnames(balances) <- rownames(innode.pos)
  
  p.val <- c()
  for (i in 1:ncol(balances)) {
    aov.data <- data.frame(balance = balances[,i], group = d.groups)
    # anova
    # res <- aov(balance ~ group, data=aov.data)
    # res <- aov(group ~ balance, data=aov.data)
    # p.val <- c(p.val,summary(res)[[1]][1,5])
    
    mod <- lm(group ~ balance, data=aov.data)
    res <- summary(mod)
    p.val <- c(p.val, res$coefficients[2,4])
    
  }
  p.val[is.na(p.val)] <- 1
  if(adjust.pvalues){
    p.adj <- p.adjust(p.val, method = 'fdr')
  } else {
    p.adj = p.val
  }
  
  px.init <- ggdend(dend.data, plot.theme=plot.theme, font.size=font.size)
  
  if(sum(p.adj < p.threshold) == 0)
    return(px.init)
  
  df.pval <- data.frame()
  df.bals <- data.frame()
  df.bal.range <- data.frame()
  df.bal.quartile <- data.frame()
  df.bal.median <- data.frame()
  df.bal.range <- data.frame()
  group.levels <- c(ref.level, target.level)
  
  for(id.node in 1:ncol(balances)){
    
    if(p.adj[id.node] < p.threshold){
      df.pval <- rbind(df.pval, innode.pos[id.node, c('x', 'y')])
    }else
      next
    
    # Normalization of values of balances according to the  range between nodes
    x.tmp <- balances[,id.node]
    x.tmp <- x.tmp - mean(x.tmp)
    
    # DO NOT MOVE THIS LINE DOWN:
    df.bal.range <- rbind(df.bal.range, data_frame(x = innode.pos$x[id.node] + innode.pos$range[id.node]/2,
                                                   y = innode.pos$y[id.node],
                                                   val = max(abs(x.tmp))))
    
    x.tmp <- x.tmp / max(abs(x.tmp)) / 2 * innode.pos$range[id.node] * 0.9
    x.tmp <- x.tmp + innode.pos$x[id.node]
    y.tmp <- d.groups * 0.03 + innode.pos$y[id.node] - 0.05
    
    q.case <- quantile(x.tmp[d.groups], c(0.25, 0.75))
    q.control <- quantile(x.tmp[!d.groups], c(0.25, 0.75))
    
    df.bal.quartile <- rbind(df.bal.quartile,
                             data_frame(x = c(q.case, q.control),
                                        y = c(rep(y.tmp[d.groups][1],2), rep(y.tmp[!d.groups][1],2)),
                                        group = group.levels[1 + c(T, T, F, F)], node = id.node))
    
    df.bal.median <- rbind(df.bal.median,
                           data_frame(x = c(median(x.tmp[d.groups]), median(x.tmp[!d.groups])),
                                      y = c(y.tmp[d.groups][1], y.tmp[!d.groups][1]),
                                      group = group.levels[1 + c(T, F)], node = id.node))
    
    
    
    df.bals <- rbind(df.bals, data.frame(x=x.tmp, y=y.tmp, group=group.levels[1 + d.groups], node = id.node))
    
  }
  
  px <- px.init + geom_point(data = df.bals,
                             aes(x=x, y=y, col = as.factor(group), group=as.factor(node)), alpha = 0.1, size = 1) +
    geom_point(data = df.bal.median,
               aes(x=x, y=y, col = as.factor(group)),
               size = 2.5, shape = 18) +
    geom_point(data = df.pval, aes(x=x, y=y)) +
    geom_line(data=df.bal.quartile, aes(x = x, y = y,
                                        col = as.factor(group),
                                        group=interaction(group, node)), size = 0.75) +
    geom_text(data=df.bal.range, mapping=aes(x=x, y=y, label=sprintf('%2.1f', val)), vjust=0, hjust=0, size=font.size) +
    labs(col=" ")
  
  return(px)
}


plotCellLoadings <- function(loadings, pval, signif.threshold, jitter.alpha, palette,
                             show.pvals, ref.level, target.level, plot.theme,
                             ordering=c("pvalue", "loading"), ref.load.level=0) {
  ordering <- match.arg(ordering)
  yintercept <- ref.load.level
  
  loading.order <- order(abs(rowMeans(loadings)))
  loadings <- loadings[loading.order, ]
  if (!is.null(pval)) {
    # if some p-values are the same - then order by mean, therefore a prior sorting is required
    pval <- pval[loading.order]
    
    if (ordering == "pvalue") {
      # additional ordering by p-value
      loadings <- loadings[order(-pval), ]
      pval <- pval[order(-pval)]
    }
    
    # Get significant cells
    n.significant.cells <- sum(pval <= signif.threshold)
  } else {
    n.significant.cells <- 0
  }
  
  res.ordered <- t(loadings) %>% as.data.frame()
  ymax <- max(loadings)
  
  p <- ggplot(stack(res.ordered), aes(x = ind, y = values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA) + geom_jitter(aes(x = ind, y = values), alpha = jitter.alpha, size=1) +
    geom_hline(yintercept = yintercept, color = "grey37") +
    coord_flip() + xlab('') + ylab('loadings') + plot.theme + theme(legend.position = "none") +
    scale_x_discrete(position = "top") + ylim(-1, 1)
  
  # Add text
  p <- p +
    annotate('text', x = 1, y = -ymax, label = paste('\u2190', ref.level), hjust = 'left') +
    annotate('text', x = 1, y = ymax, label = paste(target.level, '\u2192'), hjust = 'right')
  
  if (!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  if ((n.significant.cells > 0) && (ordering == "pvalue")) {
    p <- p + geom_vline(xintercept=nrow(loadings) - n.significant.cells + 0.5, color='red')
  }
  
  
  if (show.pvals) {
    d <- data.frame(y=-log(pval, base=10), x=names(pval), row=1:length(pval))
    p.pval <- ggplot(d, aes(x=reorder(x, row), y=y, fill=factor(x))) +
      geom_bar(stat="identity") +
      geom_hline(yintercept=-log(signif.threshold, base=10)) +
      scale_y_continuous(expand=c(0, 0)) +
      coord_flip() + labs(x='', y='-log10(adj. p-value)') +
      plot.theme + theme(legend.position = "none") + theme(axis.text.y = element_blank())
    
    if (!is.null(palette)) p.pval <- p.pval + scale_fill_manual(values=palette)
    
    p.combo <- cowplot::plot_grid(p, p.pval, nrow=1, rel_widths=c(2, 1))
    return(p.combo)
  }
  
  return(p)
  
}