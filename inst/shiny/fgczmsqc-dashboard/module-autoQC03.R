#R

source("helpers-ggplot2.R")

autoQC03UI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      htmlOutput(ns("ui")),
      width = 12)
  )
}

#' autoQC01Server shiny module
#' 
#' @description
#' handles DDA and DIA 
#' @importFrom rawDiag rawDiagServer
#'
#' @param id 
#' @param filterValues time and instrument information
autoQC03Server <- function(id, filterValues, BFabric, inputfile, readFUN, ggplot2FUN, title, footer){
  moduleServer(id,
               function(input, output, session) {
                 ns <- NS(id)
                 
                 rootdirraw <- reactive({
                   config <- configServer("config2")
                   config$rootdirraw
                 })
                 
                 vals <- reactiveValues(fn = NA,
                                        rawfile = NULL,
                                        mZ = .iRTmz(),
                                        hover = NA,
                                        userawrr = FALSE,
                                        userawDiag = FALSE)
                 
                 rawrrServer("rawrr-autoQC03", vals) 
                 rawDiag::rawDiagServer("rawDiag-autoQC03", vals)
                 
                 output$rawrrEnableUI <- renderUI({
                   shiny::req(vals$fn)
                   if (file.exists(vals$fn)){
                     rawrrUI(NS(id, "rawrr-autoQC03"))
                   }else{
                     warning(paste0("raw file ", vals$fn, " does not exist."))
                   }
                 })
                 
                 output$rawDiagEnableUI <- renderUI({
                   shiny::req(vals$fn)
                   if (file.exists(vals$fn)){
                     column(11, offset = 1,
                            fluidRow(
                              rawDiag::rawDiagUI(NS(id, "rawDiag-autoQC03"))
                            ))
                   }else{
                     warning(paste0("raw file ", vals$fn, " does not exist."))
                   }
                 })
                 
                 
                 output$hoverInfo <- renderUI({
                   shiny::req(vals$hover)
                   L <- tagList()
                   if (is.na(vals$hover$File.Name)){
                     L <- HTML("Hover over a point to see file information. <br>
                               Click on a point to trace peptides and gather raw file header information.")
                   }else{
                     L <- renderTable(t(vals$hover[, c('File.Name', 'time')]),
                                      rownames = T, colnames = F)
                   }
                   
                   shinydashboard::box(title = 'File Information',
                                       L,
                                       status = .fileStatus(rootdirraw(), vals$hover$File.Name),
                                       solidHeader = TRUE,
                                       collapsible = FALSE,
                                       width = 12,
                                       height = 150)
                 })
                 
                 observeEvent(input$rawrr, {
                   vals$userawrr <- input$rawrr
                   message("vals$rawrr: ", vals$userawrr)
                 })
                 
                 observeEvent(input$rawDiag,{
                   vals$userawDiag <- input$rawDiag
                   message("vals$rawDiag: ", vals$userawDiag)
                 })
                 observeEvent(input$hoverInfo, {
                   if(!is.null(input$hoverInfo)){
                     np <- nearPoints(data(), input$hoverInfo, xvar = "time",
                                      yvar = "value")
                     
                     cat("Hover (throttled):\n")
                     str(np$File.Name[1])
                     vals$hover <- np[1, ]
                   }
                 })
                 ## observeEvent =================
                 observeEvent(input$plotClick, {
                   np <- nearPoints(data(), input$plotClick, xvar = "time",
                                    yvar = "value")
                   
                   message(paste0("plotClick: ", np$filename, collapse = " | "))
                   
                   
                   vals$fn <- file.path(rootdirraw(), np$File.Name[1])
                   ## the rawDiag module will read the rawfile
                   vals$rawfile <- vals$fn
                 })
                 
                 
                 data <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = paste0("Reading ", title, "data ..."))
                   on.exit(progress$close())
                   
                   readFUN(inputfile)
                 })
                 
                 dataVariables <- reactive({
                   data()$variable |> unique()
                 })
                 output$variable <- renderUI({
                   dataVariables() %in% c('nConfidentProteins',
                                          'nConfidentPeptides',
                                          'nMS2', 'Precursors.Identified',
                                          'Proteins.Identified', 'AUC.lg2') |> 
                     which() -> defaulVariablesIdx 
                   
                   shinydashboard::box(
                     selectInput(ns('variables'), 'Variables',
                                 dataVariables(),
                                 multiple = TRUE,
                                 selected = dataVariables()[defaulVariablesIdx]),
                     footer = paste0("Select variables to ", title, "plots"),
                     status = "primary",
                     solidHeader = TRUE,
                     collapsible = FALSE,
                     width = 12,
                     height = 150,
                   )
                 })
                 dataFiltered <- reactive({
                   shiny::req(data())
                   
                   filter <- data()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin < data()$time & data()$time < filterValues$timeMax 
                   
                   if ("peptide" %in% colnames(data())){
                     #  message("filtering peptides ...")
                     filter <- data()$Instrument %in% filterValues$instrument &
                       filterValues$timeMin < data()$time & data()$time < filterValues$timeMax &
                       data()$peptide %in% filterValues$peptide
                   }
                   #S<-data()
                   #base::save(S, file = '/tmp/RRR.RData')
                   data()[filter, ] 
                 })
                 
                 nFacets <- reactive({
                   n <- 1
                   if ('peptide' %in% colnames(dataFiltered())){
                     n <- n * length(filterValues$peptide)
                   }
                   message("nFacets: ", n)
                   n
                 })
                 
                 output$ui <- renderUI({
                   
                   tl <- tagList(
                     fluidRow(
                       column(6,
                              tagList(
                                htmlOutput(ns("variable")),
                                fluidRow(
                                  shinydashboard::box(footer="Select modules to run when clicking",
                                                      tagList(
                                                        column(4, checkboxInput(ns('rawrr'), 'rawrr', value = vals$userawrr)),
                                                        column(4, checkboxInput(ns('rawDiag'), 'rawDiag', value = vals$userawDiag)),
                                                        column(4, checkboxInput(ns('timsR'), 'timsR', value = FALSE)),
                                                        
                                                      ),
                                                      status = "primary",
                                                      solidHeader = TRUE,
                                                      width = 12,
                                  )
                                )
                              )),
                       column(6,
                              htmlOutput(ns("hoverInfo"))
                       ),
                       shinydashboard::box(plotOutput(ns("plot"), 
                                                      height = 300 * nFacets(),
                                                      hover = hoverOpts(NS(id, "hoverInfo"),
                                                                        delay = 200, nullOutside = TRUE),
                                                      click = NS(id, "plotClick")),
                                           
                                           status = "primary",
                                           solidHeader = TRUE,
                                           collapsible = FALSE,
                                           width = 12)
                     )
                   )# tagList
                   if (filterValues$useBFabric){
                     tl <- append(tagList(htmlOutput(ns("instrumentEvents"))), tl)
                   }
                   
                   
                   
                   if (filterValues$instrument %in% c("TIMSTOF_1")){
                     ## TODO(cp): make a module module-tims.R
                     tl 
                   }else{
                     ## input$rawrr
                     if (vals$userawrr){
                       tl <- tagList(tl, htmlOutput(NS(id, "rawrrEnableUI")))
                     }
                     if (vals$userawDiag){
                       tl <- tagList(tl, 
                                     shinydashboard::box(title = "rawDiag module",
                                                         status = "primary",
                                                         solidHeader = TRUE,
                                                         collapsible = TRUE,
                                                         width = 12,
                                                         offset = 0,
                                                         htmlOutput(NS(id, "rawDiagEnableUI"))
                                     ))
                     }
                   }
                   shinydashboard::box(
                     title = paste0(title, " plots - mtime: ", file.mtime(inputfile) |> strftime("%a %F %T")),
                     footer = footer,
                     status = .status(inputfile),
                     solidHeader = TRUE,
                     collapsible = TRUE,
                     width = 12,
                     tl
                   )
                 })
                 ## render instrumentEvents
                 ## TODO(cp): add module parameter when the table should be collapsed
                 output$instrumentEvents <- renderUI({
                   shinydashboard::box(title = "Instrument events",
                                       status = "primary",
                                       solidHeader = TRUE,
                                       collapsible = TRUE,
                                       width = 12,
                                       tagList(
                                         DT::renderDataTable({ BFabric$bfabricInstrumentEventsFiltered() })
                                       ))
                 })
                 
                 ## -----------click table ----------------
                 output$plot <- renderPlot({
                   shiny::req(dataFiltered(), input$variables)
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = paste0("ggplotting ", title, "data ..."))
                   on.exit(progress$close())
                   
                   message(paste0("nrow dataFiltered(): ", nrow(dataFiltered())))
                   
                   if (nrow(dataFiltered()) == 0){
                     .missing()
                     return()
                   }
                   
                   ggplot2FUN(dataFiltered(), input$variables) -> gp
                   
                   if (filterValues$useBFabric){
                     gp + ggplot2::geom_vline(xintercept =  BFabric$bfabricInstrumentEventsFiltered()$time,
                                              linetype = "dashed", 
                                              color = "red",
                                              size = 1) -> gp
                   }
                   gp
                 }, res = 96)
                 return(data)
               })
}
