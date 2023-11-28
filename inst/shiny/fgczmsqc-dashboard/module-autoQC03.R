#R


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
#' 
#'
#' @param id 
#' @param filterValues time and instrument information
autoQC03Server <- function(id, filterValues, BFabric, inputfile, readFUN, title, footer){
  moduleServer(id,
               function(input, output, session) {
                 ns <- NS(id)
                 
                 rootdirraw <- reactive({
                   config <- configServer("config2")
                   config$rootdirraw
                 })
                 vals <- reactiveValues(fn = NA, mZ = .iRTmz(), hover = NA)
                 
                 rawrrServer("rawrr-autoQC03", vals) 
                 
                 output$rawrrEnableUI <- renderUI({
                   shiny::req(vals$fn)
                   if (file.exists(vals$fn)){
                     rawrrUI(NS(id, "rawrr-autoQC03"))
                   }else{
                     warning(paste0("raw file ", vals$fn, " does not exist."))
                   }
                 })
                 .fileStatus <- function(p, f){
                   message(f)
                   if (is.na(f)){
                     "primary"
                   }
                   else if (file.exists(file.path(p, f))){
                     "success"
                   }else{
                     "danger"
                   }
                 }
                 output$hoverInfo <- renderUI({
                   shiny::req(vals$hover)
                   L <- tagList()
                   if (is.na(vals$hover$File.Name)){
                     L <- HTML("Hover over a point to see file information. <br>
                               Click on a point to trace peptides and gather raw file header information.")
                   }else{
                     L <- renderTable(t(vals$hover[, c('File.Name', 'time')]), rownames = T, colnames = F)
                   }
                   
                   shinydashboard::box(title = 'file information',
                                       L,
                                       status = .fileStatus(rootdirraw(), vals$hover$File.Name),
                                       solidHeader = TRUE,
                                       collapsible = FALSE,
                                       width = 12)
                 })
                 observeEvent(input$hoverInfo, {
                   if(!is.null(input$hoverInfo)){
                     np <- nearPoints(data(), input$hoverInfo, xvar = "time", yvar = "value")
                     
                     cat("Hover (throttled):\n")
                     str(np$File.Name[1])
                     vals$hover <- np[1, ]
                   }
                 })
                 ## observeEvent =================
                 observeEvent(input$plotClick, {
                   np <- nearPoints(data(), input$plotClick, xvar = "time", yvar = "value")
                   
                   message(paste0("plotClick: ", np$filename, collapse = " | "))
                   
                   
                   vals$fn <- file.path(rootdirraw(), np$File.Name[1])
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
                                          'Proteins.Identified', 'FWHM.RT') |> 
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
                   )
                 })
                 dataFiltered <- reactive({
                   shiny::req(data())
                   #message(paste0("nrow data(): ", nrow(data())))
                   filter <- data()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin < data()$time & data()$time < filterValues$timeMax 
                   #data()$variable %in% input$Variables 
                   
                   data()[filter,]
                 })
                 output$ui <- renderUI({
                   
                   tl <- tagList(
                     fluidRow(
                       column(6,
                              htmlOutput(ns("variable")),
                       ),
                       column(6,
                              htmlOutput(ns("hoverInfo"))
                       ),
                       
                       shinydashboard::box(plotOutput(ns("plot"),
                                                      hover = hoverOpts(NS(id, "hoverInfo"), delay = 200, nullOutside = TRUE),
                                                      click = NS(id, "plotClick")),
                                           
                                           status = "primary",
                                           solidHeader = TRUE,
                                           collapsible = FALSE,
                                           width = 12)
                     )) # tagList
                   if (filterValues$useBFabric){
                     tl <- append(tagList(htmlOutput(ns("instrumentEvents"))), tl)
                   }
                   tagList(
                     shinydashboard::box(
                       title = paste0(title, " plots - mtime: ", file.mtime(inputfile) |> strftime("%a %F %T")),
                       footer = footer,
                       status = .status(inputfile),
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       width = 12,
                       tagList(tl,
                               htmlOutput(NS(id, "rawrrEnableUI"))
                       ),
                     ),
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
                   
                   dataFiltered() |> 
                     subset(variable %in% input$variables) |> 
                     ggplot2::ggplot(ggplot2::aes(time, value)) +
                     ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1)  -> gp
                   
                   if ("scanType" %in% names(dataFiltered())){
                     gp +
                       ggplot2::geom_point(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) +
                       ggplot2::geom_line(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) -> gp
                     
                   }else{
                     gp +
                       ggplot2::geom_point(ggplot2::aes(time, value), alpha = 0.4) +
                       ggplot2::geom_line(ggplot2::aes(time, value), alpha = 0.4) -> gp
                   }
                   if (filterValues$useBFabric){
                     gp + ggplot2::geom_vline(xintercept = BFabric$bfabricInstrumentEventsFiltered()$time,
                                              linetype="dashed", 
                                              color = "red", size = 1) -> gp
                   }
                   gp
                 }, res = 96)
                 return(data)
               })
}
