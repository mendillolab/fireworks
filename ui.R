library(DT) #
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
pool <<- dbPool(
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

############# User interface #########################
ui <- fluidPage(theme=shinytheme("paper"),

  # import javascript & CSS
  tags$head(includeHTML("www/ga.html")),
  includeScript("www/canvas2svg.js"),
  includeScript("www/script.js"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  HTML('<div id="network_svg_canvas"></div>'),

  # use navigation bar layout
  navbarPage(id="navbar",

  # add fireworks logo to navbar
  title=div(img(src="network-icon.png",
                style="margin-top: -14px; padding-right:5px; padding-bottom:15px",
                height = 60)),






############################### UI functions for network panel ##########################
          tabPanel("Network",

          #imitate sidebar layout using fluidRow
          fluidRow(

          # sidepanel
           column(4, # used to be: sidebarPanel(width="4",
              fluidRow(
              wellPanel(

              # Gene input
              selectizeInput("genes_selected",
                          label = "Source node(s):",
                          choices = c("C16orf72"),
                          selected = c("C16orf72"),
                          # choices = c("EIF2AK1","EIF2AK2","EIF2AK3","EIF2AK4","TP53",
                          #               "ATM","HSF1","HSF2","ATF5","ERN1","XBP1","EIF2AK3",
                          #               "ATF6","HIF1A","ARNT","EPAS1","HIF3A","ATF3","EIF2AK2",
                          #               "ATF4","ATG7","ATG5","KEAP1","NFE2L2"),
                          #   selected = c("EIF2AK1","EIF2AK2","EIF2AK3","EIF2AK4","TP53",
                          #               "ATM","HSF1","HSF2","ATF5","ERN1","XBP1","EIF2AK3",
                          #               "ATF6","HIF1A","ARNT","EPAS1","HIF3A","ATF3","EIF2AK2",
                          #               "ATF4","ATG7","ATG5","KEAP1","NFE2L2"),
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
                          "Both" = "both"), selected="both", inline=TRUE))
             ),


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
                              value=FALSE))

               ), # end of conditional panel

              # Show legend
              #checkboxInput("showLegend", "Show legend", FALSE),

              # Build network
              fluidRow(
                column(6, align="center", offset=3,
                  actionBttn(inputId="buildNetwork",
                              label="Build Network",
                              style="simple",
                              color="danger",
                              size="sm"),
                  tags$style(type='text/css', "#buildNetwork {height: 50px; width=100%}") )
                )

            ) # end of well panel (for control panel)
          ), # end of control panel div

          # Legend (optional)
          # fluidRow(
          #   tags$div(id="legend",
          #     img(src="legend.png", height="200px"),
          #   )
          #  )
          ), # end of sidepanel

        # Main panel
        column(8,id="mainpanel",style="z-index:500",
          visNetworkOutput("network", width="100%", height="600px"),
          tabsetPanel(id="networkTabs")
       ) # end of main panel
     ) # end of sidebar layout
   ), # end of 'Network' tab

       ##################### UI functions for co-essential heatmap ##########################
        tabPanel("Coessentiality Heatmap",
                 sidebarLayout(
                   position="left",
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
            position = "left",
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
