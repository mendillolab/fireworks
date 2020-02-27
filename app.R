library(DT)
library(plotly)
library(readr)
library(data.table)
library(pheatmap)
library(cytoscape)
library(shinythemes)
library(htmlwidgets)
library(rmarkdown)
library(DBI)
library(RPostgreSQL)
library(shinyjs)
library(pool)
library(feather)
library(genefilter)
library(matrixStats)
library(gtable)
library(grid)
library(tibble)
library(shinyWidgets)
library(dplyr)
library(tidyr)
library(visNetwork)
library(stringi)
library(scales)
source("R/network.R")
source("R/coessentialHeatmap.R")
source("R/differentialExpression.R")

############# connect to SQL database ##################
SQLdata <- readRDS("data/sqlData.Rds")
pool <- dbPool(
  drv = dbDriver("PostgreSQL"),
  host=SQLdata[[1]],
  port=as.numeric(SQLdata[[2]]),
  dbname=SQLdata[[3]],
  user=SQLdata[[4]],
  password=SQLdata[[5]]
)

onStop(function() {
  poolClose(pool)
})

############## pre-load data & variables #########################
geneNames <- strsplit(read_lines("data/gene_names.txt"),",") %>% unlist %>% as.data.table
colnames(geneNames) <- "Gene:"

# boolean flags to check if dataset(s) have been loaded
panCorrMat <- NULL
achilles <- NULL
ccleExpr <- NULL
plotRecord <- NULL
geneinfo <- NULL

# boolean flags to check if tab 'pages' are loaded
exampleNetwork <<- TRUE
CEHM.page.loaded <- FALSE
RNA.page.loaded <- FALSE
networkTabsLoaded <- FALSE
CEHMtabsLoaded <- FALSE
RNAtabsLoaded <- FALSE

# Convenience functions for loading 'external' datasets into memory
# External datasets (not included in repo):
#     - panCancer.corr.feather (pan-cancer correlation matrix)
#     - achilles.feather (gene x cell line essentiality data)
#     - CCLE_expression.feather (CCLE_expression.feather)

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

############# style specifications (TODO: move to separate CSS script) #############
appCSS <- "
#loading-content {
  position: absolute;
  background: #FFFFFF;
  opacity: 0.9;
  left: 48%;
  top: 25%;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}

.shiny-notification{
  position: fixed;
  top: 33%;
  left: 15%;
  width: 30%;
  font-size: 15px;
}

label{font-size:12px;}
.label2 {font-size:12px;}

.selectize-input {margin-bottom: 0px;}

.aboutText {
  width: 70%;
  font-size: 20px;
  margin: 0 auto;
  line-height: 1.4;
  text-align: justify;
}

.multiOmicsWarning {
  border-radius: 10px;
  border-style: solid;
  border-width: 1px;
  border-color: #fdc1ca;
  color: #7b001d;
  text-align: justify;
  padding: 0px 20px;
  margin: 20px;
  background-color: #fed4d9;

}

#corrType{width=100%}

#network_svg_canvas {
  width: 800px;
  height: 600px;
  position: absolute !important;
  top: -9999px !important;
  left: -9999px !important;

}

.progress-bar {
  background-color: red;
}'
"
############## convenience functions #########################
splitLines <- function(line){
    line <- paste0(stri_extract_all_regex(line, '.{1,20}')[[1]], collapse="<br />")
    return(line)
}

flatten <- function(title){
  title <- stringr::str_split(string=title, pattern="<br />" %>% unlist)
  title <- lapply(title, splitLines)
  title <- unlist(new_title) %>% paste0(collapse="<br />")
  return(title)
}

############# User interface #########################
ui <- fluidPage(theme=shinytheme("paper"),

  # top navbar
  tags$head(includeHTML("www/ga.html")),
  includeScript("www/canvas2svg.js"),
  includeScript("www/script.js"),
  inlineCSS(appCSS),
  HTML('<div id="network_svg_canvas"></div>'),
      navbarPage(id="navbar",

        # add fireworks logo to navbar (TODO: change logo?)
        title=div(img(src="network-icon.png",
                      style="margin-top: -14px; padding-right:5px; padding-bottom:15px",
                      height = 60)),

          ##################### UI functions for network panel ##########################
          tabPanel("Network",
          sidebarLayout(
            position = "right",
            sidebarPanel(width="4",

              # Gene input
              selectizeInput("genes_selected",
                          label = "Source node(s):",
                          choices = c("C16orf72"),
                          selected = c("C16orf72"),
                          # choices = c("EIF2AK1","EIF2AK2","EIF2AK3","EIF2AK4","TP53",
                          #              "ATM","HSF1","HSF2","ATF5","ERN1","XBP1","EIF2AK3",
                          #              "ATF6","HIF1A","ARNT","EPAS1","HIF3A","ATF3","EIF2AK2",
                          #              "ATF4","ATG7","ATG5","KEAP1","NFE2L2"),
                          #  selected = c("EIF2AK1","EIF2AK2","EIF2AK3","EIF2AK4","TP53",
                          #              "ATM","HSF1","HSF2","ATF5","ERN1","XBP1","EIF2AK3",
                          #              "ATF6","HIF1A","ARNT","EPAS1","HIF3A","ATF3","EIF2AK2",
                          #              "ATF4","ATG7","ATG5","KEAP1","NFE2L2"),
                          multiple = TRUE),

              # Context input
              selectizeInput("context",
                          label = "Context:",
                          choices = c("Pan-Cancer", "Colon", "Urinary", "Lung",
                                      "Ovary", "Skin", "Breast", "Pancreas",
                                      "CNS","Bone","Stomach","Soft Tissue",
                                      "Autonomic Ganglia", "Kidney", "Liver",
                                      "Endometrium", "Esophagus"),
                          selected = "Pan-Cancer",
                          multiple = FALSE),

              fluidRow(

              # Correlation direction (primary connections)
              column(4,
                numericInput("k_primary",
                             label = "Rank (primary):",
                             value=30)),

              # Rank for primary connections
              column(8,
              radioButtons("corrType", "Correlation direction:",
                          c("Positive" = "pos",
                          "Negative" = "neg",
                          "Both" = "both"), selected="both", inline=TRUE))),


             # Add second order connections
             tags$div(class="label2",
             materialSwitch("second_order",
                           label="Secondary nodes",
                           status="primary",
                           value=TRUE)),

              # Show isolated nodes (no second order connections)
              conditionalPanel(
                            condition = "input.second_order == false",
                            checkboxInput("showIN",
                              label="Show isolated nodes",
                              value=FALSE)),

              # Additional options for second order connections
              conditionalPanel(
                            condition = "input.second_order == true",

                            wellPanel(style="background:#e6e6e6; border-radius: 10px; padding-top: 5px; padding-right: 5px;",

                            fluidRow(
                            # Correlation direction for second order connections
                            column(4,
                              numericInput("k_secondary",
                                label = "Rank (secondary):",
                                value=5)),

                            # Rank for second order connections
                            column(8,
                              radioButtons("corrType2", "Correlation direction",
                                          c("Positive" = "pos",
                                          "Negative" = "neg",
                                          "Both" = "both"),
                                          selected="pos",
                                          inline=TRUE))
                            ),

                            # Show isolated primary nodes
                            checkboxInput("showIsolatedPrimaryNodes",
                              label="Show isolated primary nodes",
                              value=FALSE),

                            # Show isolated secondary nodes
                            checkboxInput("showIsolatedSecondaryNodes",
                              label="Show isolated secondary nodes",
                              value=FALSE))),

              # Build network
              fluidRow(
                column(6, align="center", offset=3,
                  actionBttn(inputId="buildNetwork",
                              label="Build Network",
                              style="simple",
                              color="danger",
                              size="sm"),
                  tags$style(type='text/css', "#buildNetwork {height: 50px; width=100%}") ))

            ),

            # View the network
            mainPanel(id="mainpanel",width="8",style="z-index:500",
              visNetworkOutput("network", height="600px"),
              #HTML('<hr>'),
              #plotlyOutput('depBoxPlot'),
              #HTML('<hr>'),
              tabsetPanel(id="networkTabs")
           ) # main panel - pan-cancer
        ) # sidebar layout - pan-cancer
      ), # tab panel - network

       ##################### UI functions for co-essential heatmap ##########################
        tabPanel("Coessentiality Heatmap",
                 sidebarLayout(
                   position="right",
                   sidebarPanel(

                     # Gene input
                     selectizeInput("geneCEHM",
                                 label = "Gene:",
                                 choices = c("C16orf72"),
                                 multiple=FALSE,
                                 selected="C16orf72"),

                     # Context input
                     selectizeInput("contextCEHM",
                                 label = "Context:",
                                 choices = c("Pan-Cancer", "Colon", "Urinary", "Lung",
                                             "Ovary", "Skin", "Breast", "Pancreas",
                                             "CNS","Bone","Stomach","Soft Tissue",
                                             "Autonomic Ganglia", "Kidney", "Liver",
                                             "Aerodigestive","Endometrium", "Esophagus"),
                                 selected = "Pan-Cancer",
                                 multiple = FALSE),

                     fluidRow(

                     column(4,
                     # Correlation direction
                     numericInput("kCEHM",
                                  label = "Rank (primary):",
                                  value=10)
                                 ),

                     column(8,
                       radioButtons("directionCEHM", "Correlation direction:",
                                   c("Positive" = "pos",
                                     "Negative" = "neg",
                                     "Both" = "both"), inline=TRUE, selected="both")
                     )
                      ), # end fluidRow


                     # Build network
                     fluidRow(
                       column(6, align="center", offset=3,
                         actionBttn(inputId="goCEHM",
                                     label="Submit",
                                     style="simple",
                                     color="danger",
                                     size="sm"),
                         tags$style(type='text/css', "#goCEHM {height: 50px; width=100%}") ))

                   ),

                   mainPanel(
                     plotOutput("cehmHeatmap"),
                     tabsetPanel(id="CEHMtabs"))
                 ) # end sidebarLayout
        ), # end co-essential heatmap layout


      ##################### UI functions for multi-omics exploration ##########################
      navbarMenu("Multiomics",

        ####### UI functions
        tabPanel("RNA",
          sidebarLayout(
            position = "right",
            sidebarPanel(

              # Gene input
              selectizeInput("genesRNA",
                        label = "Gene(s):",
                        choices = c("C16orf72"),
                        selected = "C16orf72",
                        multiple = TRUE),

              # Context input
              selectizeInput("contextRNA",
                          label = "Context:",
                          choices = c("Pan-Cancer", "Colon", "Urinary", "Lung",
                                      "Ovary", "Skin", "Breast", "Pancreas",
                                      "CNS","Bone","Stomach","Soft Tissue",
                                      "Autonomic Ganglia","Kidney", "Liver",
                                      "Endometrium", "Esophagus", "Aerodigestive"),
                          selected = "Liver",
                          multiple = FALSE),

              # P-value cutoff
              numericInput("pCutoffRNA",
                    label = "P-value cutoff:",
                    value=0.005,
                    min=0,max=1),

            sliderInput("Q", "Compare top and bottom N% (most vs. least dependent). N= ",
                min = 1, max = 50,
                value = 25, step = 1),

              # Generate Signature
              fluidRow(
                column(6, align="center", offset=3,
                  actionBttn(inputId="goRNA",
                              label="Submit",
                              style="simple",
                              color="danger",
                              size="sm"),
                  tags$style(type='text/css', "#buildNetwork {height: 50px; width=100%}") )),

              # Warning
              tags$div(class="multiOmicsWarning",
              HTML('<br/><b>IMPORTANT NOTE:</b> Signatures containing
              anti-correlated genes will still group cell lines by dependence and perform
              differential expression between these groups. This may not provide biologically relevant results. </p>')
            ) # end multi-Omics multiOmicsWarning


            ), # sidebarPanel

          # main panel (Multi-omics > RNA)
          mainPanel(id="mainpanel",
            plotOutput("RNAheatmap"),
            tabsetPanel(id="RNAtabs")
          ) # main panel (Multi-omics > RNA)
        ) # sidebarLayout
      ), # end tab panel (Multi-omics > RNA)


       tabPanel("Metabolites (Coming soon)",
        h4("Coming soon")
     ) # end tab panel (Metabolites)
   ), # end navbarMenu
      tabPanel("About",
      tags$div(class="aboutText",

      # FIREWORKS
      #h3("FIREWORKS (Fitness Interaction Ranked nEtWORKS)"),
      HTML('<h3>FIREWORKS (Fitness Interaction Ranked nEtWORKS) <a href="https://github.com/mendillolab/fireworks"><img src="github_icon.png"></a></h3>
      <p>FIREWORKS is an interactive web tool developed to facilitate
      interrogation of biological networks using unbiased
      genetic screens. Gene-gene relationships can be modeled by correlating essentiality
      scores across hundreds of cancer cell lines. This application allows you to visualize these
      relationships, build networks from them and leverage multi-omics data to investigate their underlying biology.</p>'),


      # Getting started
      h3("Getting started"),
      tags$ul(
        tags$li(HTML('<p><b>Network:<br /></b>
                      Any list of genes can be used as the
                      foundation of a network. These genes are referred to as
                      “source nodes”. The most coessential genes for each
                      source node form the “primary nodes".
                      The most coessential genes for each primary
                      node can be added by selecting “secondary nodes”.
                      <br /><Br />Adding secondary nodes increases the compute
                      time, but often results in genes organizing into functional
                      modules. Frequently, these networks return nodes that are not
                      connected to any other sub-networks. This can result in a
                      cluttered looking visualization. To simplify these networks,
                      we have provided an option to remove isolated primary and
                      secondary nodes.

                      </p>')),
        tags$li(HTML('<p><b>Coessentiality Heatmap: <br /></b>
                      In the “coessentiality heatmap” tab, you can investigate
                      the relative gene essentiality of a gene of interest and
                      its top coessential partners across cell lines. You must
                      specify the context, correlation direction and
                       rank cut-off for primary connections. Second order
                       connections are not included.
                      </p>')),
        tags$li(HTML('<p><b>Multi-omics: <Br /></b>
                      What determines which cell lines are dependent on a given
                      gene signature? The “multi-omics” tab aids in answering
                      this question, by integrating fitness data with other
                      functional characterizations of cancer cell lines.
                      Currently, you can use the “RNA” tool to perform
                      differential expression between the cell lines that are
                      the most and least dependent on a given gene, or
                      gene signature. The analysis depends on the provided context, p-value
                      cutoff and quantile cut-off. <br /><br />

                      <font color="#c72100"><b>IMPORTANT NOTE: </b>Signatures
                      containing anti-correlated genes will still group cell
                      lines by dependence and perform differential expression between these
                      groups. This may not provide biologically
                      relevant results. </font>
                      </p>'))
      ),

      # Contact us
      h3("Contact Us"),
      HTML('<p>This tool is a work in progress. Please feel free to contact us with any questions,
      comments or other related inquiries. We will get back to you as soon as possible (<a href="jasen.jackson@northwestern.edu">jasen.jackson@northwestern.edu</a>).
      </p>'),

      # Read the paper
      h3("Read the paper"),
      HTML('<p><b>Coessential Gene Networks Reveal the Organization and Constituents
      of a Dynamic Cellular Stress Response</b><br />
      David R Amici, Jasen M Jackson, Kyle A Metz, Daniel J Ansel, Byoung-Kyu Cho,
      Roger S Smith, Sonia Brockway, Seesha Takagishi, Shashank Srivastava,
      Brendan P. O’Hara, Young Ah Goo, Neil L Kelleher, Issam Ben-Sahra, Daniel
      R Foltz, Marc L Mendillo<br />
      <a href="https://www.biorxiv.org/content/10.1101/847996v1">Biorxiv 2019</a>
      </p><br /><br />')
       # END un-ordered list
    ) # end aboutText div

      )
    ) # end navbarPage
  #)) # end div, hidden
) # end fluid page

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

    visNetwork(nodes, edges, height="600px") %>%
         visNodes(shape = "dot",
                  labelHighlightBold=FALSE,
                  font = list(vadjust=-50,
                              size=30,
                              strokeWidth=0.5,
                              strokeColor='#000000')) %>%
         visGroups(groupname = "source", size=45, font=list(size=45,vadjust=-85)) %>%
         visEdges(smooth=TRUE, length=10) %>%
         visInteraction(keyboard = TRUE,
                        tooltipDelay=800) %>%
         visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
         visIgraphLayout(randomSeed = 6, smooth=TRUE, type="full")

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

# Run app ----
shinyApp(ui, server)
