#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(guessMyLOI)
# library(pheatmap)
# require(biovizBase)
# require(GenomicRanges)
# require(ggrepel)
# library(ggbio)


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Guess My LOI"),
  selectizeInput("species", label=NULL, choices = c("Hsapiens","Mmusculus"),
                           options = list(create = TRUE, placeholder = 'Hsapiens')),
  fluidRow(column(6, fileInput("file1", "Choose File",
                               accept = c("text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv/.txt/.ae"))),
           column(6, numericInput("depth", "Min Depth", 20, min = 2, max = 100))
           ),
  fluidRow(column(6, numericInput("countsToCall", "Min Counts in Minor allele To Call", 10, min = 2, max = 20)),
           column(6, numericInput("callheteroThr", "Minor / (Major + Minor) to call hetero", 0.2, min = 0.05, max = 0.45, step=0.05))
           ),

  tabsetPanel(type = "tabs",
              tabPanel("LOI",
                       selectizeInput("gene", label = NULL, choices = NULL,
                                      options = list(create = TRUE, placeholder = 'select a gene name')),
                       splitLayout(cellWidths = c("25%", "75%"),
                                   plotOutput("LOIheatmap", height="1200px"),
                                   plotOutput("LOIgene"))
              ),
              tabPanel("XCI",
                       selectizeInput("geneX", label = NULL, choices = NULL,
                                      options = list(create = TRUE, placeholder = 'select a gene name')),
                       splitLayout(cellWidths = c("25%", "75%"),
                                   plotOutput("XCIheatmap", height="1200px"),
                                   plotOutput("LOIgeneX"))
              ),
              tabPanel("LOI Gene Tab",
                       # downloadButton("downloadDataG", "Download"),
                       dataTableOutput("tableGene")),
              tabPanel("LOI SNP Tab",
                       # downloadButton("downloadDataS", "Download"),
                       dataTableOutput("tableSnp")),
              tabPanel("XCI Gene Tab",
                       # downloadButton("downloadDataXG", "Download"),
                       dataTableOutput("tableGeneX")),
              tabPanel("XCI Gene Tab",
                       # downloadButton("downloadDataXS", "Download"),
                       dataTableOutput("tableSnpX"))
  )
)

server <- function(input, output, session) {
  readData <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)

    annotCol = 7
    readGATK.AllelicImbalance(file=inFile$datapath, annotCol = annotCol)
  })

  processTheData <- reactive({
    data <- readData()
    if (is.null(data))
      return(NULL)

    annotCol = 7
    depth = input$depth
    if (is.null(depth))
      return(NULL)

    countsToCall = input$countsToCall
    callheteroThr = input$callheteroThr
    total=TRUE
    species = tolower(input$species)

    gatk_ai <- filter_by_overall_depth(data, depth = depth)
    S <- processNonSexualData(gatk_ai, countsToCall, total, callheteroThr, species, removeNonDBsnp = TRUE)
    X <- processSexualData(gatk_ai, countsToCall, total, callheteroThr, species, removeNonDBsnp = TRUE)

    updateSelectizeInput(session, 'gene', choices = S$genes$geneAnnot$gene, server = TRUE)
    updateSelectizeInput(session, 'geneX', choices = X$genes$geneAnnot$gene, server = TRUE)
    return(list(S=S, X=X))
  })

  output$LOIheatmap <- renderPlot({
    data <- processTheData()
    gr = NULL
    if (is.null(data))
      return(NULL)
    if (input$species=="mmusculus") {
      gr <- 93
    }
    plot_guessedLOI(data$S$genes, gaps_row = gr, plot=T)
  })

  output$LOIgene <- renderPlot({
    input$gene
    data <- processTheData()
    if (is.null(data) | is.null(input$gene))
      return(NULL)
    snpsGeneAi <- data$S$snps
    keep = snpsGeneAi$annot$gene == input$gene
    g <- filter_rows(snpsGeneAi, keep)
    plotAllelicRatios(g, na.rm = F)
  })

  output$XCIheatmap <- renderPlot({
    data <- processTheData()
    if (is.null(data))
      return(NULL)
    plot_guessedLOI(data$X$genes, plot=T)
  })

  output$LOIgeneX <- renderPlot({
    input$geneX
    data <- processTheData()
    if (is.null(data) | is.null(input$gene))
      return(NULL)
    snpsGeneAi <- data$X$snps
    keep = snpsGeneAi$annot$gene == input$geneX
    g <- filter_rows(snpsGeneAi, keep)
    plotAllelicRatios(g, na.rm = F)
  })

  output$tableGene <- renderDataTable({
    data <- processTheData()
    if (is.null(data))
      return(NULL)
    data$S$geneSummary
    }, options=list(pageLength=10))

  output$tableSnp <- renderDataTable({
    data <- processTheData()
    if (is.null(data))
      return(NULL)
    data$S$snpSummary
  }, options=list(pageLength=10))

  output$tableGeneX <- renderDataTable({
    data <- processTheData()
    if (is.null(data))
      return(NULL)
    data$X$geneSummary
  }, options=list(pageLength=10))

  output$tableSnpX <- renderDataTable({
    data <- processTheData()
    if (is.null(data))
      return(NULL)
    data$X$snpSummary
  }, options=list(pageLength=10))

  # output$downloadDataG <- downloadHandler(
  #   filename = function() {
  #     paste("geneSummary", ".txt", sep = "")
  #   },
  #   content = function(file) {
  #     data <- processTheData()
  #     write.csv(data$S$geneSummary, file, row.names = FALSE)
  #   }
  # )
  #
  # output$downloadDataS <- downloadHandler(
  #   filename = function() {
  #     paste("spnSummary", ".txt", sep = "")
  #   },
  #   content = function(file) {
  #     data <- processTheData()
  #     write.csv(data$S$snpSummary, file, row.names = FALSE)
  #   }
  # )
  #
  # output$downloadDataXG <- downloadHandler(
  #   filename = function() {
  #     paste(X-geneSummary, ".txt", sep = "")
  #   },
  #   content = function(file) {
  #     data <- processTheData()
  #     write.csv(data$X$geneSummary, file, row.names = FALSE)
  #   }
  # )
  #
  # output$downloadDataXG <- downloadHandler(
  #   filename = paste("X-spnSummary", ".txt", sep = ""),
  #   content = function(file) {
  #     data <- processTheData()
  #     write.csv(data$X$spnSummary, file, row.names = FALSE)
  #   }
  # )
}

# Run the application
shinyApp(ui = ui, server = server)

