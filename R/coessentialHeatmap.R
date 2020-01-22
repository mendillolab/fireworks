# Script: coessentialHeatmap.R
# Author: Jasen Jackson

draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates(length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- grid::textGrob(
      coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
      vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
    )
    return(res)
}

plotHeatmapCEHM <- function(heatmap, meta.cols=NULL, meta.rows=NULL) {
  # rotate column labels
  assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
  )

  # plot heatmap
  main_title = paste("Relative gene essentiality across", ncol(heatmap), "cancer cell lines", sep=" ")
  colors = list(Direction = c(positive="green", negative="yellow"))
  pheatmap(heatmap,
           main=main_title,
           fontsize=11,
           col= colorRampPalette(c("blue", "white", "red"))(n = 299),
           cluster_rows=FALSE,
           cluster_cols=TRUE,
           annotation_col=meta.cols,
           annotation_row=meta.rows,
           annotation_legend=TRUE,
           annotation_colors = colors,
           annotation_names_row=FALSE,
           scale='row',
           #border_color=NA,
           breaks=seq(from=-1.5,to=1.5,length.out=300),
           show_rownames=isTRUE(nrow(heatmap) < 25),
           show_colnames=isTRUE(ncol(heatmap) < 55))
} # plotHeatmap()

coessentialHeatmapSQL <- function(table,pool,achilles,gene,pos,neg,rank){
  SQLtable <- tbl(pool, table)
  gene.corrMat <- SQLtable %>% filter(gene1 == gene | gene2 == gene) %>% collect()
  codeps <- getCodeps(corrmat=gene.corrMat, gene=gene, pos=pos, neg=neg, k=rank)

  # get list of co-essential genes
  posGenes <- NULL
  negGenes <- NULL
  if (!is.null(codeps[[1]])){posGenes <- codeps[[1]] %>% select(gene) %>% unlist}
  if (!is.null(codeps[[2]])){negGenes <- codeps[[2]] %>% select(gene) %>% unlist}
  coessentialGenes <- c(posGenes, negGenes)

  # format heatmap output
  achilles.filtered <- achilles %>% filter(lineage==table)
  hm.mtx <- achilles.filtered[,c(gene, coessentialGenes)] %>% as.matrix
  rownames(hm.mtx) <- achilles.filtered[,"stripped_cell_line_name"] %>% unlist
  hm.mtx <- t(hm.mtx) * -1

  # add column annotations (lineage)
  hm.meta.cols <- achilles.filtered[,c("stripped_cell_line_name","lineage")]  %>%
    tibble::column_to_rownames("stripped_cell_line_name")

  # add row annotations (pos vs. neg correlation)
  hm.meta.rows <- data.frame(Direction=c(rep("positive", length(posGenes)+1), rep("negative", length(negGenes))))
  rownames(hm.meta.rows) <- c(gene, coessentialGenes)

  # format table for output
  tbl.output <- bind_rows(codeps[[1]], codeps[[2]])
  tbl.output <- tbl.output %>% select(-gene1, -gene2)
  tbl.output <- tbl.output[,c("gene","cor","p")]
  return(list(hm.mtx, hm.meta.rows, tbl.output))
}

coessentialHeatmap_local <- function(corrMat, achilles, sourceGenes, gene, pos, neg, rank){
  gene.corrMat <- corrMat %>% filter(gene1 == gene | gene2 == gene)
  codeps <- getCodeps(corrmat=gene.corrMat, gene=gene, pos=pos, neg=neg, k=rank)

  # get list of co-essential genes
  posGenes <- NULL
  negGenes <- NULL
  if (!is.null(codeps[[1]])){posGenes <- codeps[[1]] %>% select(gene) %>% unlist}
  if (!is.null(codeps[[2]])){negGenes <- codeps[[2]] %>% select(gene) %>% unlist}
  coessentialGenes <- c(posGenes, negGenes)

  # format heatmap output
  hm.mtx <- achilles[,c(gene, coessentialGenes)] %>% as.matrix
  rownames(hm.mtx) <- achilles[,"stripped_cell_line_name"] %>% unlist
  hm.mtx <- t(hm.mtx) * -1

  # add column annotations (lineage)
  hm.meta.cols <- achilles[,c("stripped_cell_line_name","lineage")]  %>%
    tibble::column_to_rownames("stripped_cell_line_name")

  # add row annotations (pos vs. neg correlation)
  hm.meta.rows <- data.frame(Direction=c(rep("positive", length(posGenes)+1), rep("negative", length(negGenes))))
  rownames(hm.meta.rows) <- c(gene, coessentialGenes)

  # format table for output
  tbl.output <- bind_rows(codeps[[1]], codeps[[2]])
  tbl.output <- tbl.output %>% select(-gene1, -gene2)
  tbl.output <- tbl.output[,c("gene","cor","p")]
  return(list(hm.mtx, hm.meta.rows, tbl.output))
}
