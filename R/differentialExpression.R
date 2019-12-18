# Script: differentialExpression.R
# Author: Jasen Jackson
# Description: Script contains functions to perform differential expression
#              between cell lines that differentially depend on a provided gene
#              or set of genes
# require(dplyr)
# require(matrixStats)
# require(genefilter)

# plot heat map
plotHeatmapDE <- function(heatmap.dat, meta.cols=NULL) {

  # rotate column labels
  assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
  )

  # plot heatmap
  colors = list(Dependence = c(high="green", low="yellow"))
  pheatmap(heatmap.dat,
           col= colorRampPalette(c("blue", "white", "red"))(n = 299),
           cluster_rows=FALSE,
           fontsize=11,
           cluster_cols=TRUE,
           scale='row',
           annotation_col=meta.cols,
           annotation_names_col=FALSE,
           annotation_colors= colors,
           #border_color=NA,
           breaks=seq(from=-1.5,to=1.5,length.out=300),
           show_rownames=isTRUE(nrow(heatmap.dat) < 35),
           show_colnames=isTRUE(ncol(heatmap.dat) < 100))
} # plotHeatmap()


# requires getLineageName() function ~ in network.R
differentialExpression <- function(achilles, ccleExpr, genes, context, pCutoff=0.01, Q=25){

  n_tiles = floor(100/Q)
  if (context == "Pan-Cancer"){
    # Pan-Cancer DE (TODO: put repeated code in function)# get data frame for "average" signature
    signature <- achilles[,c(genes)] %>%
      mutate_all(rank) %>%
      as.matrix %>% rowMeans
    signature.df <- as_tibble(cbind(achilles[,"DepMap_ID"], as.data.frame(signature))) %>%
      arrange(signature)

    # Rank cell line IDs (HIGH = most essential, LOW = least essential)
    signature.df <- signature.df %>% mutate(quantile=ntile(signature,n_tiles))

    # get IDs for cell lines highly dependent on signature
    high_ceres <- signature.df %>%
      filter(quantile==1) %>%
      select(DepMap_ID) %>%
      unlist

    # get IDs for cell lines lowly dependent on signature
    low_ceres <- signature.df %>%
      filter(quantile==n_tiles) %>%
      select(DepMap_ID) %>%
      unlist %>% rev

    # Filter to expression of test/control cell lines
    achilles <- achilles %>% tibble::column_to_rownames("DepMap_ID") # for naming
    RNA.high <- ccleExpr %>% filter(DepMap_ID %in% high_ceres)
    RNA.high.mtx <- RNA.high[,2:19142] %>% as.matrix %>% t
    colnames(RNA.high.mtx) <- achilles[RNA.high[,"DepMap_ID"] %>% unlist,"stripped_cell_line_name"]

    RNA.low <- ccleExpr %>% filter(DepMap_ID %in% low_ceres)
    RNA.low.mtx <- RNA.low[,2:19142] %>% as.matrix %>% t
    colnames(RNA.low.mtx) <-  achilles[RNA.low[,"DepMap_ID"] %>% unlist,"stripped_cell_line_name"]

    # Exclude genes with mean(TPM) < 0.1
    RNA.high.means <- rowMeans(RNA.high.mtx)
    RNA.low.means <- rowMeans(RNA.low.mtx)

    high.filtered.genes <- RNA.high.means[RNA.high.means > 0.1] %>% names
    low.filtered.genes <- RNA.low.means[RNA.low.means > 0.1] %>% names
    filtered.genes <- intersect(high.filtered.genes, low.filtered.genes)

    RNA.high.mtx <- RNA.high.mtx[c(filtered.genes),]
    RNA.low.mtx <- RNA.low.mtx[c(filtered.genes),]

    # Run differential expression
    test.matrix <- cbind(RNA.high.mtx, RNA.low.mtx)
    g <- factor( c(rep(0, ncol(RNA.high.mtx)), rep(1, ncol(RNA.low.mtx))) )
    test <- rowttests(test.matrix, g)
    test[,"gene"] <- rownames(test)
    test <- test %>% filter(p.value<pCutoff) %>% arrange(desc(statistic))
    test <- test[,c("gene","dm","statistic","p.value")]
    colnames(test) <- c("gene","mean.high_dep.minus.mean.low_dep","statistic","p.value")

    # Make df for table/heatmap (faster way to compute this?)
    diffGenes <- test[,"gene"] %>% unlist
    results.high <- RNA.high.mtx[c(diffGenes),] # c(high_ceres[1:4])
    results.low  <-  RNA.low.mtx[c(diffGenes),] # rev(c(low_ceres[1:4]))

    # add column annotations for heatmap
    hm.meta.cols <- data.frame(Dependence=c(rep("high", ncol(results.high)),
                                            rep("low", ncol(results.low))))
    rownames(hm.meta.cols) <- c(colnames(results.high), colnames(results.low))

    heatmap.data <- cbind(results.high, results.low)
  } else {
    # context-specific DE
    achilles.filtered <- achilles %>% filter(lineage==getLineageName(context))

    # get data frame for "average" signature
    signature <- achilles.filtered[,c(genes)] %>%
      mutate_all(rank) %>%
      as.matrix %>% rowMeans
    signature.df <- as_tibble(cbind(achilles.filtered[,"DepMap_ID"], as.data.frame(signature))) %>%
      arrange(signature)

    # Rank cell line IDs (HIGH = most essential, LOW = least essential)
    signature.df <- signature.df %>% mutate(quantile=ntile(signature,n_tiles))

    # get IDs for cell lines highly dependent on signature
    high_ceres <- signature.df %>%
      filter(quantile==1) %>%
      select(DepMap_ID) %>%
      unlist

    # get IDs for cell lines lowly dependent on signature
    low_ceres <- signature.df %>%
      filter(quantile==n_tiles) %>%
      select(DepMap_ID) %>%
      unlist %>% rev

    # Filter to expression of test/control cell lines
    achilles.filtered <- achilles.filtered %>% tibble::column_to_rownames("DepMap_ID") # for naming
    RNA.high <- ccleExpr %>% filter(DepMap_ID %in% high_ceres)
    RNA.high.mtx <- RNA.high[,2:19142] %>% as.matrix %>% t
    colnames(RNA.high.mtx) <- achilles.filtered[RNA.high[,"DepMap_ID"] %>% unlist,"stripped_cell_line_name"]

    RNA.low <- ccleExpr %>% filter(DepMap_ID %in% low_ceres)
    RNA.low.mtx <- RNA.low[,2:19142] %>% as.matrix %>% t
    colnames(RNA.low.mtx) <-  achilles.filtered[RNA.low[,"DepMap_ID"] %>% unlist,"stripped_cell_line_name"]

    # Exclude genes with mean(TPM) < 0.1
    RNA.high.means <- rowMeans(RNA.high.mtx)
    RNA.low.means <- rowMeans(RNA.low.mtx)

    high.filtered.genes <- RNA.high.means[RNA.high.means > 0.1] %>% names
    low.filtered.genes <- RNA.low.means[RNA.low.means > 0.1] %>% names
    filtered.genes <- intersect(high.filtered.genes, low.filtered.genes)

    RNA.high.mtx <- RNA.high.mtx[c(filtered.genes),]
    RNA.low.mtx <- RNA.low.mtx[c(filtered.genes),]

    # Run differential expression
    test.matrix <- cbind(RNA.high.mtx, RNA.low.mtx)
    g <- factor( c(rep(0, ncol(RNA.high.mtx)), rep(1, ncol(RNA.low.mtx))) )
    test <- rowttests(test.matrix, g)
    test[,"gene"] <- rownames(test)
    test <- test %>% filter(p.value<pCutoff) %>% arrange(desc(statistic))
    #test <- test %>% select(-dm)
    test <- test[,c("gene","dm","statistic","p.value")]
    colnames(test) <- c("gene","mean.high_dep.minus.mean.low_dep","statistic","p.value")

    # Make df for table/heatmap (faster way to compute this?)
    diffGenes <- test[,"gene"] %>% unlist
    results.high <- RNA.high.mtx[c(diffGenes),] # c(high_ceres[1:4])
    results.low  <-  RNA.low.mtx[c(diffGenes),] # rev(c(low_ceres[1:4]))

    # add column annotations for heatmap
    hm.meta.cols <- data.frame(Dependence=c(rep("high", ncol(results.high)),
                                        rep("low", ncol(results.low))))
    rownames(hm.meta.cols) <- c(colnames(results.high), colnames(results.low))

    heatmap.data <- cbind(results.high, results.low)
  }


  #  Return test output and heatmap data
  return(list(heatmap.data, test, hm.meta.cols))
} # differentialExpression()
