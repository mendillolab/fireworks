############## pre-load data & variables #########################
geneNames <<- strsplit(read_lines("data/gene_names.txt"),",") %>% unlist %>% as.data.table
colnames(geneNames) <<- "Gene:"

# boolean flags to check if dataset(s) have been loaded
panCorrMat <<- NULL
achilles <<- NULL
ccleExpr <<- NULL
plotRecord <<- NULL
geneinfo <<- NULL

# boolean flags to check if tab 'pages' are loaded
exampleNetwork <<- TRUE
CEHM.page.loaded <<- FALSE
RNA.page.loaded <<- FALSE
networkTabsLoaded <<- FALSE
CEHMtabsLoaded <<- FALSE
RNAtabsLoaded <<- FALSE

# Convenience functions for loading 'external' datasets into memory
readCorrMat <- function(session, panCorrMat) {
  progress <- Progress$new(session)
  progress$set(value = 0.5, message = 'Loading correlation matrix...')
  panCorrMat <<- read_feather("data/panCancer.corr.feather")
  progress$set(value = 1, message = 'Loading...')
  progress$close()
}

readAchilles <- function(session, achilles) {
  progress <- Progress$new(session)
  progress$set(value = 0.5, message = 'Loading genetic dependency data...')
  achilles <<- read_feather("data/achilles.feather")
  progress$set(value = 1, message = 'Loading...')
  progress$close()
}

readGeneInfo <- function(session, geneinfo) {
  progress <- Progress$new(session)
  progress$set(value=0.5, message = 'Loading gene annotation data...')
  geneinfo <<- read_feather("data/geneinfo.all.feather")
  progress$set(value = 1, message = 'Loading..')
  progress$close()
}

readExpression <- function(session, ccleExpr) {
  progress <- Progress$new(session)
  progress$set(value = 0.5, message = "Loading expression data...")
  ccleExpr <<- read_feather("data/CCLE_expression.feather")
  progress$set(value = 1, message = 'Loading...')
  progress$close()
}

# Server logic ----
server <- function(input, output, session) {

  ##################### server functions for network ##########################

  # build network with parameters
  codep_network <- eventReactive(input$buildNetwork,
                                 ignoreNULL=FALSE,{

         # load genetic dependency data
         if (is.null(achilles)){readAchilles(session, achilles)}

         # load gene annotation data
         if (is.null(geneinfo)){readGeneInfo(session, geneinfo)}

          # show tabset panel for output
          if (networkTabsLoaded==TRUE){
            removeTab(inputId="networkTabs", target="Nodes")
            removeTab(inputId="networkTabs", target="Edges")
            #removeTab(inputId="networkTabs", target="Essentiality")
            removeTab(inputId="networkTabs", target="Download")
            removeTab(inputId="networkTabs", target="Legend")
            networkTabsLoaded <<- FALSE
          }

          if (networkTabsLoaded==FALSE){
            appendTab(inputId="networkTabs", tab=tabPanel("Nodes", DT::dataTableOutput("nodesTable")))
            appendTab(inputId="networkTabs", tab=tabPanel("Edges", DT::dataTableOutput("edgesTable")), select = TRUE)
            #appendTab(inputId="networkTabs", tab=tabPanel("Essentiality",
            #                                               HTML('<br />'),
            #                                               plotlyOutput('depBoxPlot')), select = TRUE)
            appendTab(inputId="networkTabs", tab=tabPanel("Download",
                                                          HTML('<br />'),
                                                          selectInput("chooseDatasetNetwork", "Choose a dataset:",
                                                          choices = c("nodes(.csv)", "edges(.csv)"), selected="image(.svg)"),
                                                          downloadButton("downloadNetwork", "Download"),
                                                          actionButton("downloadVectorImage", "Download Vector Image"),

                                                          HTML('<br /><br /><br />')))
            appendTab(inputId="networkTabs", tab=tabPanel("Legend",
                                                          HTML('<br />'),
                                                          tags$div(id="legend2",
                                                            img(src="legend.png", width="300px", align="middle")
                                                          ),
                                                          HTML('<br /><br />')))
            networkTabsLoaded <<- TRUE
          }

          # read user parameters from sidebar panel
          sourceGenes <- input$genes_selected
          context <- input$context
          k1 <- input$k_primary
          direction1 <- input$corrType
          secondOrder <- input$second_order
          if (secondOrder) {
            k2 <- input$k_secondary
            showIPN <- input$showIsolatedPrimaryNodes
            showISN <- input$showIsolatedSecondaryNodes
            direction2 <- input$corrType2
            if (direction2=="pos") {pos2=TRUE; neg2=FALSE}
            if (direction2=="neg") {pos2=FALSE; neg2=TRUE}
            if (direction2=="both") {pos2=TRUE; neg2=TRUE}
          } else{showIPN <- input$showIN}

          if (direction1=="pos") {pos1=TRUE; neg1=FALSE}
          if (direction1=="neg") {pos1=FALSE; neg1=TRUE}
          if (direction1=="both") {pos1=TRUE; neg1=TRUE}

          # use pre-loaded matrix for pan-cancer analyses
          if (context=="Pan-Cancer") {

            # load pan-cancer correlation matrix
            if(is.null(panCorrMat)){readCorrMat(session, panCorrMat)}

            # build network from corr matrix
            if (secondOrder) {
              network <- buildNetwork_local(corrMat=panCorrMat, sourceGenes=sourceGenes, k1=k1, k2=k2,
                                      pos1=pos1, neg1=neg1, pos2=pos2, neg2=neg2,
                                      showIPN=showIPN, showISN=showISN, secondOrder=TRUE,
                                      exampleNetwork=exampleNetwork, geneinfo=geneinfo, achilles=achilles)
              exampleNetwork <<- FALSE
            } else {
              network <- buildNetwork_local(corrMat=panCorrMat, sourceGenes=sourceGenes, k1=k1,
                                      pos1=pos1, neg1=neg1, showIPN=showIPN, secondOrder=FALSE,
                                      exampleNetwork=exampleNetwork, geneinfo=geneinfo, achilles=achilles)
              exampleNetwork <<- FALSE
            }
          }

          # use SQL database for lineage-specific analyses
          else {
            print("Use SQL database")
            if (secondOrder) {
              network <- buildNetworkSQL2(table=getLineageName(context), pool=pool, sourceGenes=sourceGenes, k1=k1, k2=k2,
                                      pos1=pos1, neg1=neg1, pos2=pos2, neg2=neg2,
                                      showIPN=showIPN, showISN=showISN, secondOrder=TRUE,
                                      geneinfo=geneinfo, achilles=achilles)
            } else {
              network <- buildNetworkSQL2(table=getLineageName(context), pool=pool, sourceGenes=sourceGenes, k1=k1,
                                      pos1=pos1, neg1=neg1, showIPN=showIPN, secondOrder=FALSE,
                                      geneinfo=geneinfo, achilles=achilles)
            }
          }

         # return computed network
         network
       })

  # visualize network using vizNetwork
  output$network <- renderVisNetwork({

      # get network
      nodes <- codep_network()[[1]]
      edges <- codep_network()[[2]]


      # create 'width' column (for edges)
      sourceNodes <- nodes %>% filter(type=="source") %>% select(gene) %>% unlist
      primaryEdges <- edges %>% filter(   source %in% sourceNodes  )
      otherEdges  <- edges %>% filter( !(source %in% sourceNodes) )

      ## edge width by correlation quantile
      corrValues <- primaryEdges[,'correlation'] %>% unlist %>% as.double %>% abs
      primaryEdges[,'width'] <- rescale(corrValues, to = c(0.5, 5)) %>% as.double
      otherEdges[,'width'] <- rep(2, nrow(otherEdges))

      # add dashes for secondary edges
      primaryEdges[,'dashes'] <- rep(FALSE, nrow(primaryEdges))
        otherEdges[,'dashes'] <- rep(TRUE,  nrow(otherEdges))
      edges <- rbind(primaryEdges, otherEdges)

      # create 'label' column (for nodes)
      nodes[,'label'] <- nodes[,'gene']

      # create group column (for nodes)
      nodes[,'group'] <- nodes[,'type']

      # create title column (for nodes)
      nodes <- nodes %>% mutate(title=paste0("<p><b><a href='https://www.ncbi.nlm.nih.gov/gene/",id,"\' target='_blank'>",gene,"</a></b><br />",
                                              #aliases,"<br />",
                                              ifelse( ((aliases %>% is.na) | (nchar(aliases)==0)),"",paste0(aliases,"<br />") ),
                                              "<i>",name,"</i><br />",
                                              "<b>Dependency:</b><br />",
                                              "<b>Min: </b>",round(min,2)," <b>Median: </b>",round(median,2)," <b>Max: </b>",round(max,2),"<br />",
                                              "<b>Q25: </b>",round(Q25,2)," <b>Q75: </b>",round(Q75,2),"<br />"
                                              ))
      # create id > gene map.
      nodes[,"id"] <- c(1:nrow(nodes))
      id2geneMap <- nodes[,"id"] %>% unlist
      names(id2geneMap) <- nodes[,"gene"] %>% unlist

      # use map to make 'from'/'to' columns
      from <- id2geneMap[edges[,'source'] %>% unlist]
      to <- id2geneMap[edges[,'target'] %>% unlist]
      edges <- cbind(edges, from, to)

      print(dim(nodes))

    visNetwork(nodes, edges) %>%
         visNodes(shape = "dot",
                  labelHighlightBold=FALSE,
                  font = list(vadjust=-50,
                              size=30,
                              strokeWidth=0.5,
                              strokeColor='#000000')) %>%
         visGroups(groupname = "source", size=45, font=list(size=45,vadjust=-85)) %>%
         visPhysics(solver='barnesHut'

                    #repulsion=list(nodeDistance=1000,
                                  #springConstant=0.01)
                                ) %>%
         visEdges(smooth=list(enabled=TRUE,roundness=0)) %>%
         visInteraction(keyboard = TRUE,
                        tooltipDelay=800) %>%
         visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
         #visLegend(addEdges = ledges, addNodes = lnodes, useGroups = FALSE, width=0.3, ncol=1, zoom=TRUE) %>%
         visIgraphLayout(randomSeed = 7, smooth=TRUE, type="full", physics=FALSE, layout="layout_with_fr")

  })

# generate essentiality data
#   output$depBoxPlot <- renderPlotly({
#
#     # get network data
#     nodes <- codep_network()[[1]]
#     edges <- codep_network()[[2]]
#
#     # reshape essentiality data for plotting
#     achillesPlot <- achilles[,c("stripped_cell_line_name",
#                                 "disease",
#                                 nodes %>% select(gene) %>% unlist)] %>%
#       gather(nodes %>% select(gene) %>% unlist,
#              key="gene", value="dependency")
#
#     ## add color/name data to plot dataframe
#     # create map gene --> color
#     colMap <- nodes[,"color"]
#     names(colMap) <- nodes[,"gene"]
#
#     # create color column & add to plot df
#     color = colMap[achillesPlot[,'gene'] %>% unlist]
#     achillesPlot <- cbind(achillesPlot, color)
#
#     # plot boxplot
#     depBoxPlot <- plot_ly(achillesPlot,
#                           x = ~gene,
#                           y = ~dependency,
#                           color = ~color,
#                           type = 'box',
#                           text = ~paste(stripped_cell_line_name,
#                                         '<br>Disease:', disease),
#                           hoverinfo="text") %>%
#                   plotly::layout(title = "Essentiality of genes in network")
#   })

  output$nodesTable <- renderDataTable({
    nodes <- codep_network()[[1]][,c("gene","type","name","median","min","max")]
    datatable(nodes,
              options = list(pageLength=10, lengthChange=FALSE),
              rownames=FALSE)
  })

  output$edgesTable <- renderDataTable({
    edges <- codep_network()[[2]] %>% select(-color)
    datatable(edges,
              options = list(pageLength=15, lengthChange=FALSE),
              rownames=FALSE)
  })

  # Reactive value holding the selected dataset for download
  datasetNetwork <- reactive({
    switch(input$chooseDatasetNetwork,
           "nodes" = codep_network()[[1]],
           "edges" = codep_network()[[2]]
           #"vector" = TRUE
          )
  })

  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste(input$chooseDatasetNetwork, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetNetwork(), file, row.names = FALSE)
    }
  )

  ### load gene names for drop down
  updateSelectizeInput(session, "genes_selected", choices = geneNames, selected=c("C16orf72"))
  updateSelectizeInput(session, "geneCEHM", choices = geneNames, selected=c("C16orf72"))
  updateSelectizeInput(session, "genesRNA", choices = geneNames, selected=c("C16orf72"))

  # svg function 1: collect network proxy data
  observeEvent(input$downloadVectorImage,{
      visNetworkProxy("network") %>% visGetNodes()
      visNetworkProxy("network") %>% visGetEdges()
  })

  # svg function 2: send network data to javascript listener
  observe({
    if(!is.null(input$network_nodes)){
    # send network data
    session$sendCustomMessage("svg_handler_nodes", input$network_nodes)
    session$sendCustomMessage("svg_handler_edges", input$network_edges)
    message = "Network data sent to console"
    session$sendCustomMessage("svg_handler_update", message)
    }

  })

  ##################### end of server functions for co-essentiality network ###################
  ##################### beginning of server functions for co-essential heatmap ##########################

  # observe: clicking on co-essential heatmap
  check.cehm <- reactiveVal()
  observeEvent(input$navbar,{
    if ((req(input$navbar)=="Coessentiality Heatmap") & CEHM.page.loaded == FALSE){
      check.cehm(check.cehm(0)+1)
      CEHM.page.loaded <<- TRUE
    }
  })

  # update co-essential gene list
  CEHM.output <- eventReactive(c(input$goCEHM,check.cehm()),
                              ignoreNULL=TRUE,{

    if (CEHMtabsLoaded==TRUE){
      removeTab(inputId="CEHMtabs", target="Table")
      removeTab(inputId="CEHMtabs", target="Download")
      CEHMtabsLoaded <<- FALSE
    }

    if (CEHMtabsLoaded==FALSE){
      appendTab(inputId="CEHMtabs", tab=tabPanel("Table", DT::dataTableOutput("cehmTable")), select=TRUE)
      appendTab(inputId="CEHMtabs", tab=tabPanel("Download",
                                                  HTML('<br />'),
                                                  selectInput("chooseDatasetCEHM", "Choose a file:",
                                                  choices = c("coessential_genes.csv", "heatmap_data.csv")),
                                                  downloadButton("downloadCEHM", "Download"),
                                                  HTML('<br /><br /><br />')))

      CEHMtabsLoaded <<- TRUE
    }

    # Read user parameters from sidebar panel
    sourceGene <- input$geneCEHM
    context <- input$contextCEHM
    rank <- input$kCEHM
    direction <- input$directionCEHM
    if (direction=="pos"){pos=TRUE; neg=FALSE}
    if (direction=="neg"){pos=FALSE; neg=TRUE}
    if (direction=="both"){pos=TRUE; neg=TRUE}

    # Load required data
    if (is.null(achilles)){readAchilles(session, achilles)}

    # if pan-cancer: use pre-loaded matrix for pan-cancer analyses
    if (context=="Pan-Cancer"){
      if (is.null(panCorrMat)){readCorrMat(session, panCorrMat)}
      CEHM.output <- coessentialHeatmap_local(corrMat=panCorrMat,achilles=achilles,
                                                gene=sourceGene, pos=pos, neg=neg,rank=rank)
    } else {
      CEHM.output <- coessentialHeatmapSQL(table=getLineageName(context), pool=pool,
                                           achilles=achilles, gene=sourceGene,
                                           pos=pos, neg=neg, rank=rank)}

    # return output for co-essentiality heatmap
    CEHM.output

  }) # CEHM.output

  output$cehmHeatmap <- renderPlot({
    plotHeatmapCEHM(heatmap=CEHM.output()[[1]], meta.rows=CEHM.output()[[2]])
  })

  output$cehmTable <- renderDataTable({
    datatable(CEHM.output()[[3]],
              options = list(pageLength=10, lengthChange=FALSE),
              rownames=FALSE)
  })

  # Reactive value for selected dataset
  datasetCEHM <- reactive({
    switch(input$chooseDatasetCEHM,
           "coessential_genes.csv" = CEHM.output()[[3]],
           "heatmap_data.csv" = CEHM.output()[[1]])
  })

  output$downloadCEHM <- downloadHandler(
    filename = function() {
      input$chooseDatasetCEHM
    },
    content = function(file) {
      write.csv(datasetCEHM(), file, row.names = TRUE)
    }
  )


  ##################### server functions for differential expression (RNA) ##########################

  # observe: clicking on co-essential heatmap
  check.RNA <- reactiveVal()
  observeEvent(input$navbar,{
    if ((req(input$navbar)=="RNA")& RNA.page.loaded == FALSE){
      check.RNA(check.RNA(0)+1)
      RNA.page.loaded <<- TRUE
    }
  })

  diffExp <- eventReactive(c(input$goRNA, check.RNA()),
                           ignoreNULL=TRUE,{

     if (RNAtabsLoaded==TRUE){
       removeTab(inputId="RNAtabs", target="Table")
       removeTab(inputId="RNAtabs", target="Download")
       RNAtabsLoaded <<- FALSE
     }

     if (RNAtabsLoaded==FALSE){
       appendTab(inputId="RNAtabs", tab=tabPanel("Table", DT::dataTableOutput("RNAtable")), select=TRUE)
       appendTab(inputId="RNAtabs", tab=tabPanel("Download",
                                                   HTML('<br />'),
                                                   selectInput("chooseDatasetRNA", "Choose a file:",
                                                   choices = c("differential_expression_output.csv",
                                                               "expression_heatmap_data.csv")),
                                                   downloadButton("downloadRNA", "Download"),
                                                   HTML('<br /><br /><br />')))

       RNAtabsLoaded <<- TRUE
     }

    # Load required data
    if (is.null(ccleExpr)){readExpression(session, ccleExpr)}
    if (is.null(achilles)){readAchilles(session, achilles)}

    # Get DE parameters from input variables
    genesRNA = input$genesRNA
    contextRNA = input$contextRNA
    pCutoffRNA = input$pCutoffRNA
    Q = input$Q

    # filter df by disease
    DE.out <- differentialExpression(achilles, ccleExpr, genesRNA, contextRNA, pCutoffRNA, Q)
    DE.out
  })

  output$RNAheatmap <- renderPlot({
    plot <- plotHeatmapDE(heatmap.dat=diffExp()[[1]], meta.cols=diffExp()[[3]])
  })

  output$RNAtable <- DT::renderDataTable({
    DT::datatable(diffExp()[[2]],
                  options = list(pageLength=10,
                                 lengthChange=FALSE,
                                 columns = list(
                                   list(title = 'gene'),
                                   list(title = 'mean(high-dep) - mean(low-dep)'),
                                   list(title = 't-statistic'),
                                   list(title = 'p-value'))),
                  rownames=FALSE)
  })

  datasetRNA <- reactive({
    switch(input$chooseDatasetRNA,
           "differential_expression_output.csv" = diffExp()[[2]],
           "expression_heatmap_data.csv" = diffExp()[[1]])
  })

  output$downloadRNA <- downloadHandler(
    filename = function() {
      input$chooseDatasetRNA
    },
    content = function(file) {
      write.csv(datasetRNA(), file, row.names = TRUE)
    }
  )
}
