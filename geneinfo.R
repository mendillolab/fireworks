library(feather)
library(JSONlite)

# load gene names
genes <- strsplit(read_lines("data/gene_names.txt"),",") %>% unlist

# initialize dataframe
geneinfo.df <- data.frame(id=character(),
                          gene=character(),
                          name=character(),
                          aliases=character())

# add each gene to dataframe
genes <- genes.all[3500:length(genes.all)]
count=0
for (gene in genes) {

  # progress update
  count = count + 1
  if (count %% 100 == 0){print(paste0(count, " genes annotated."))}

  # get gene annotation data
  query = paste0("http://mygene.info/v3/query?q=",gene,"&species=human&fields=name,id,alias")
  geneJSON <- fromJSON(query)

  # create new row based on annotation data
  if (is.null(geneJSON$max_score)){
    geneinfo <- data.frame("gene"=gene, "id"="","name"="","aliases"="")
  } else{
    geneinfo <- data.frame("gene"=gene,
                            "id"=geneJSON$hits$'_id'[1],
                            "name"=geneJSON$hits$name[1],
                            "aliases"=paste0(geneJSON$hits$alias[1] %>% unlist, collapse=","))
  }

  # add row to dataframe
  geneinfo.df <- rbind(geneinfo.df, geneinfo)
}
write_feather(geneinfo.df, path="data/geneinfo.feather")


getGeneInfo <- function(nodes, edges){
f <- function(gene) {
  query = paste0("http://mygene.info/v3/query?q=",gene,"&species=human&fields=name,id,alias")
  geneJSON <- fromJSON(query)
  # geneinfo <- fromJSON("http://mygene.info/v3/query?q=HUWE1&species=human&fields=name,id,genomic_pos_hg19")
  geneinfo <- data.frame("gene"=gene,
                          "id"=geneJSON$hits$'_id'[1],
                          "name"=geneJSON$hits$name[1],
                          "aliases"=paste0(geneJSON$hits$alias[1] %>% unlist, collapse=","))
  return(geneinfo)
}
