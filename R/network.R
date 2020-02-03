# Script: networkSQL.R
# Author: Jasen Jackson
# Description: Script contains functions used to build fitness interaction
#              networks, using connection to external SQL database. SQL
#              access info removed.
#
#
# require(RPostgreSQL)
# require(dbConnection)
# require(DPLYR)
# require(JSONlite)

getLineageName <- function(choice){

  # TODO: convert to dictionary-like object
  if (choice=="Colon"){return("colorectal")}
  if (choice=="Urinary"){return("urinary_tract")}
  if (choice=="Lung"){return("lung")}
  if (choice=="Ovary"){return("ovary")}
  if (choice=="Skin"){return("skin")}
  if (choice=="Breast"){return("breast")}
  if (choice=="Pancreas"){return("pancreas")}
  if (choice=="CNS"){return("central_nervous_system")}
  if (choice=="Bone"){return("bone")}
  if (choice=="Stomach"){return("gastric")}
  if (choice=="Autonomic Ganglia"){return("peripheral_nervous_system")}
  if (choice=="Kidney"){return("kidney")}
  if (choice=="Aerodigestive"){return("upper_aerodigestive")}
  if (choice=="Liver"){return("liver")}
  if (choice=="Endometrium"){return("uterus")}
  if (choice=="Esophagus"){return("esophagus")}
  if (choice=="Soft Tissue"){return("soft_tissue")}

} # getLineageNameFunction

getCodeps <- function(corrmat, gene, pos, neg, k){

  # get positive co-dependent genes
  posCodeps <- c()
  if (pos==TRUE){
    posCodeps <- corrmat %>%
     filter(cor > 0) %>% arrange(desc(cor)) %>%
     head(k) %>%
     mutate(gene=ifelse(gene1==gene, gene2, gene1))
  }

  # grab negative co-dependent genes
  negCodeps <- c()
  if (neg==TRUE){
    negCodeps <- corrmat %>%
    filter(cor < 0) %>% arrange(cor) %>%
    head(k) %>% mutate(gene=ifelse(gene1==gene, gene2, gene1))
  }

  return(list(posCodeps, negCodeps))
} # getCodeps function

updateNetwork <- function(nodes, edges, codeps, gene){

  # prepare data for addition to network
  posGenes <- NULL
  posCors <- NULL
  posPvalues <- NULL
  negGenes <- NULL
  negCors <- NULL
  negPvalues <- NULL

  if (!is.null(codeps[[1]])){
    posGenes <- codeps[[1]] %>% select(gene) %>% unlist
    posCors <- codeps[[1]] %>% select(cor) %>% unlist
    posPvalues <- codeps[[1]] %>% select(p) %>% unlist
  }
  if (!is.null(codeps[[2]])){
    negGenes <- codeps[[2]] %>% select(gene) %>% unlist
    negCors <- codeps[[2]] %>% select(cor) %>% unlist
    negPvalues <- codeps[[2]] %>% select(p) %>% unlist
  }

  coessentialGenes <- c(posGenes, negGenes)
  cors <- c(posCors, negCors)
  pValues <- c(posPvalues, negPvalues)
  edgeColors <- c( rep("#ff5f62",length(posGenes)),
                   rep("#5fafff", length(negGenes)) )
  newNodes <- setdiff(coessentialGenes, nodes[,"gene"])

  # update nodes
  if (length(newNodes)>0){
    nodes <- bind_rows(nodes, data.frame(
      gene=newNodes,
      type=rep("secondary",length(newNodes)),
      color=rep("#dedede", length(newNodes)) ) %>%
      mutate_all(as.character)
    )
  }

  # update edges
  edges <- bind_rows(edges, data.frame(
    source=rep(gene, length(coessentialGenes)),
    target=coessentialGenes,
    correlation=cors,
    p=pValues,
    color=edgeColors)
    %>% mutate_all(as.character)
  )

  # return updated network
  return(list(nodes, edges))

} # updateNetwork

# remove isolated secondary nodes
removeISN <- function(nodes, edges, source_genes, color) {
    isolatedNodes <- c()
    secondaryNodes = nodes %>% filter(color=="#dedede")
    secondaryNodes = secondaryNodes[,"gene"]

    for (gene in secondaryNodes){
      if(edges %>% filter(source == gene | target == gene) %>% nrow < 2){
        isolatedNodes <- c(gene, isolatedNodes)}}

    isolatedNodes <- unique(isolatedNodes)
    nodes.filtered <- nodes %>% filter(!(gene %in% isolatedNodes))
    edges.filtered <- edges %>% filter(!(target %in% isolatedNodes))
    network <- list(nodes.filtered, edges.filtered)
}

# identify isolated primary nodes
findIPN <- function(nodes, edges, primaryNodes) {

  # find isolated primary nodes hierarchically
  isolatedPrimaryNodes <- c()
  for (gene in primaryNodes) {

    # gather edges and initialize boolean flags for each primary node
    geneEdges = edges %>% filter(source==gene | target==gene)
    oneParentalEdge = FALSE
    isIPN = FALSE

    # if the primary node has only one parental edge
    oneParentalEdge <- (geneEdges %>% filter(target==gene) %>% nrow == 1)
    if (oneParentalEdge) {

      # then.. test each descendant for connectivity.
      isIPN <- TRUE
      childNodes <- geneEdges %>% filter(target != gene) %>% select(target) %>% unlist
      for (childNode in childNodes) {

        # if it's connected to anything else, the primary node isn't isolated (TODO: not *anything* else.. what if it connects back to sourceGene?)
        if (edges %>% filter(target==childNode|source==childNode, source!=gene) %>% nrow > 0) {
          isIPN <- FALSE
          break}}

      # if it passes all of these tests (or fails them?) it's an IPN.
      if (isIPN) {isolatedPrimaryNodes <- c(gene, isolatedPrimaryNodes)}
    }
  } # end for loop

  return(isolatedPrimaryNodes)
}

# remove isolated primary nodes
removeIPN <- function(nodes, edges, source_genes, color) {

  # find isolated primary nodes
  primaryNodes <- nodes %>% filter(type=="primary") %>% select(gene) %>% unlist
  isolatedPrimaryNodes <- findIPN(nodes, edges, primaryNodes)

  # filter edges
  edges.filtered <- edges %>% filter(!(source %in% isolatedPrimaryNodes),
                                     !(target %in% isolatedPrimaryNodes))
  # filter nodes not in edges
  keptNodes <- c(edges.filtered[,"source"] %>% unlist, edges.filtered[,"target"] %>% unlist) %>% unique
  nodes.filtered <- nodes %>% filter(gene %in% keptNodes)

  network <- list(nodes.filtered, edges.filtered)
  return(network)
}

getGeneInfo <- function(nodes, edges, geneinfo, achilles, context="pan_cancer") {

  #### merge nodes df with gene annotation dataframe
  nodes <- merge(nodes, geneinfo, by="gene")

  #### merge nodes df with dependency data
  genes <- nodes[,"gene"] %>% unlist
  if (context != "pan_cancer") {
    # TODO: test this!
    achilles.nodes <- achilles %>% filter(lineage==context)
    achilles.nodes <- achilles.nodes[,genes]
  } else { achilles.nodes <- achilles[,genes] }

  ## flip signs on achilles dataframe (more positive = more essential)
  achilles.nodes <- achilles.nodes*(-1)

  # grab median, min, max, and 25th/75th quartiles
  medianApply <- function(x){return(median(x, na.rm=TRUE))}
  minApply <- function(x){return(min(x, na.rm=TRUE))}
  maxApply <- function(x){return(max(x, na.rm=TRUE))}
  quantile25 <- function(x){return(quantile(x,probs=0.25, na.rm=TRUE))}
  quantile75 <- function(x){return(quantile(x,probs=0.75, na.rm=TRUE))}
  nodes.withAchilles <- data.frame(gene=genes,
                              median=apply(achilles.nodes, 2, medianApply),
                              min=apply(achilles.nodes,2,minApply),
                              max=apply(achilles.nodes,2,maxApply),
                              Q25=apply(achilles.nodes,2,quantile25),
                              Q75=apply(achilles.nodes,2,quantile75))

  # merge with other nodes w/ nodes.withAchilles
  nodes <- merge(nodes, nodes.withAchilles, by="gene")

  return(list(nodes, edges))
}


# (IN PROGRESS) build a network using pre-loaded correlation matrix (for pan-cancer analyses)
buildNetwork_local <- function(corrMat, sourceGenes, k1=10,
                         k2=0, pos1=TRUE, neg1=TRUE, secondOrder=FALSE,
                         pos2=TRUE, neg2=FALSE,
                         showIPN=TRUE, showISN=TRUE, exampleNetwork=FALSE,
                         geneinfo, achilles) {
  #' @corrMat(df): melted correlation matrix (row, column, cor, p)
  #' @gene(chr): list of codependent genes, chr
  #' @k1(int): rank cutoff for primary connection
  #' @k2(int): rank cutoff for secondary connections
  #' @pos1(boolean): return positive correlations (first order connections)
  #' @pos2(boolean): return positive correlations (second order connections)
  #' @neg1(boolean): return negative correlations (first order connections
  #' @neg2(boolean): return negative correlations (second order connections)
  #' @secondOrder(boolean): include second order connections
  #' @showIPN(boolean): show isolated primary nodes
  #' @showISN(boolean): show isolated secondary nodes
  #' @geneinfo(df): gene annotation data
  #' @achilles(df): gene dependency data


  if ((sourceGenes == c("C16orf72")) &
      (pos1==TRUE) & (neg1==TRUE) &
      (pos2==TRUE) & (neg2==FALSE) &
      (secondOrder==TRUE) &
      (k1==30) & (k2==5) &
      (showIPN==FALSE) & (showISN==FALSE)){
    # load example network
    nodes <- read_csv("data/C16orf72.30.nodes.csv") #nodes2
    edges <- read_csv("data/C16orf72.30.edges.csv")
  } else {
      # initialize network
      nodes <- data.frame(gene=sourceGenes) %>% mutate_all(as.character)
      edges <- data.frame(source=sourceGenes,
                          target=sourceGenes,
                          color=rep("#", length(sourceGenes))) %>% mutate_all(as.character)

      # fill out network for each source gene
      withProgress(message = "Adding primary connections for: ", value=0,{
      n <- length(sourceGenes)
      for (gene in sourceGenes){

        # update progress bar
        incProgress(1/n, detail = paste(gene))

        # find top co-essential genes (codeps)
        gene.corrMat <- corrMat %>% filter(gene1 == gene | gene2 == gene)
        codeps <- getCodeps(corrmat=gene.corrMat, gene=gene, pos=pos1, neg=neg1, k=k1)
        network.new <- updateNetwork(nodes, edges, codeps, gene)
        nodes <- network.new[[1]]
        edges <- network.new[[2]]

      } # for (gene in sourceGenes)
    }) # withProgress
      edges <- edges %>% filter(source != target)
      nodes <- nodes %>% mutate(
        type = ifelse(gene %in% sourceGenes, "source", "primary"),
        color = ifelse(gene %in% sourceGenes, "#ffbf8f", "#e6dcfc"))

      # find secondary pools
      if (secondOrder){
        primaryGenes = edges[,'target'] %>% unique
        primaryGenes = primaryGenes[!(primaryGenes %in% sourceGenes)]
        edgesPerGene = nrow(edges) / length(sourceGenes)

        withProgress(message = "Adding second order connections for: ", value=0,{

        n <- length(primaryGenes)
        progress = 0
        sourceIndex = 1
        for (codep in primaryGenes){

        # update progress
        progress = progress + 1
        incProgress(1/n, detail = paste(sourceGenes[sourceIndex], ">", codep))
        if (progress %% edgesPerGene == 0){sourceIndex = sourceIndex + 1}


         # find top co-essential genes (codeps)
         gene.corrMat <- corrMat %>% filter(gene1 == codep | gene2 == codep)
         codeps <- getCodeps(corrmat=gene.corrMat, gene=codep, pos=pos2, neg=neg2, k=k2)
         network.new <- updateNetwork(nodes, edges, codeps, codep)
         nodes <- network.new[[1]]
         edges <- network.new[[2]]


       }}) # for (gene in primaryGenes) / }) end w/ progress
      } # if (secondOrder)

      # trim isolated nodes as indicated by showISN, showIPN
      network <- list(nodes, edges)
      if (!(showISN)){network <- removeISN(network[[1]], network[[2]], sourceGenes)}
      if (!(showIPN)){network <- removeIPN(network[[1]], network[[2]], sourceGenes)}

      # add gene annotation data
      network <- getGeneInfo(network[[1]], network[[2]], geneinfo, achilles)

      # return network
      return(network)


    } # else: end network build

    # RETURN STATEMENT FOR DEFAULT LOADED NETWORK
    network <- list(nodes, edges)
    return(network)

} # buildNetwork_local

# (COMPLETE) build a network using SQL database connection (for lineage-specific analyses)
buildNetworkSQL2 <- function(table="pan_cancer_full", pool, sourceGenes, k1=10,
                         k2=0, pos1=TRUE, neg1=TRUE, secondOrder=FALSE,
                         pos2=TRUE, neg2=FALSE,
                         showIPN=TRUE, showISN=TRUE,
                         geneinfo, achilles) {

  # initialize network with temporary node/edge (removed after the loop)
  nodes <- data.frame(gene=sourceGenes) %>% mutate_all(as.character)
  edges <- data.frame(source=sourceGenes,
                      target=sourceGenes,
                      color=rep("#ccc", length(sourceGenes))) %>% mutate_all(as.character)

  # fill out network for each source gene
  SQLtable <- tbl(pool, table)
  withProgress(message = "Adding primary connections for: ", value=0,{
  n <- length(sourceGenes)
  for (gene in sourceGenes) {

     # update progress bar
     incProgress(1/n, detail = paste(gene))

     # find top co-essential genes (codeps)
     corrmat <- SQLtable %>% filter(gene1 == gene | gene2 == gene) %>% collect()
     codeps <- getCodeps(corrmat=corrmat, gene=gene, pos=pos1, neg=neg1, k=k1)
     network.new <- updateNetwork(nodes, edges, codeps, gene)
     nodes <- network.new[[1]]
     edges <- network.new[[2]]
  } # for (gene in sourceGenes)
  }) # withProgress

  edges <- edges %>% filter(source != target)
  nodes <- nodes %>% mutate(
    type = ifelse(gene %in% sourceGenes, "source", "primary"),
    color = ifelse(gene %in% sourceGenes, "#ffbf8f", "#e6dcfc"))

  # find secondary connections
  if (secondOrder){
    primaryGenes = edges[,'target'] %>% unique
    primaryGenes = primaryGenes[!(primaryGenes %in% sourceGenes)]
    edgesPerGene = nrow(edges) / length(sourceGenes)

    withProgress(message = "Adding second order connections for: ", value=0,{
    n <- length(primaryGenes)
    progress = 0
    sourceIndex = 1
    for (codep in primaryGenes){

    # update progress
    progress = progress + 1
    incProgress(1/n, detail = paste(sourceGenes[sourceIndex], ">", codep))
    if (progress %% edgesPerGene == 0){sourceIndex = sourceIndex + 1}

     # find top co-essential genes (codeps)
     corrmat <- SQLtable %>% filter(gene1 == codep | gene2 == codep) %>% collect()
     codeps <- getCodeps(corrmat=corrmat, gene=codep, pos=pos2, neg=neg2, k=k2)
     network.new <- updateNetwork(nodes, edges, codeps, codep)
     nodes <- network.new[[1]]
     edges <- network.new[[2]]
   }

 }) # for (gene in primaryGenes) / }) end w/ progress
  } # if (secondOrder)

  # trim isolated nodes as indicated by showISN, showIPN
  network <- list(nodes, edges)
  if (!(showISN)){network <- removeISN(network[[1]], network[[2]], sourceGenes)}
  if (!(showIPN)){network <- removeIPN(network[[1]], network[[2]], sourceGenes)}

  # add gene annotation data
  network <- getGeneInfo(network[[1]], network[[2]], geneinfo, achilles)

  # return network
  return(network)

} # buildNetworkSQL2
