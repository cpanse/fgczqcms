#R


autoQC03UI <- function(id){
  ns <- NS(id)
  fluidRow(
    htmlOutput(ns("ui")),
    width = 12)
}


#' autoQC01Server shiny module
#' 
#' @description
#' handles DDA and DIA 
#' 
#'
#' @param id 
#' @param filterValues time and instrument information
autoQC03Server <- function(id, filterValues, BFabric, inputfile, readFUN, title){
  moduleServer(id,
               function(input, output, session) {
                 ns <- NS(id)
                 
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
                   #message(paste0("nrow data(): ", nrow(data())))
                   filter <- data()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin < data()$time & data()$time < filterValues$timeMax 
                   #data()$variable %in% input$Variables 
                   
                   data()[filter,]
                 })
                 
                 output$ui <- renderUI({
                   tl <- tagList(
                     htmlOutput(ns("variable")),
                     shinydashboard::box(plotOutput(ns("plot")),status = "primary",
                                         solidHeader = TRUE,
                                         collapsible = FALSE,
                                         width = 12)
                   )
                   
                   if (filterValues$useBFabric){
                     tl <- append(tl, tagList(htmlOutput(ns("instrumentEvents"))))
                   }
                   
                   tagList(
                     #h2(title),
                     shinydashboard::box(
                       title = paste0(title, " plots"),
                       status = "primary",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       width = 12,
                       tl
                     ),
                   )

                 })
                 
                 output$instrumentEvents <- renderUI({
                   #shiny::req(BFabric$bfabricInstrumentEventsFiltered())
                   shinydashboard::box(
                     title = "Instrument events",
                     status = "primary",
                     solidHeader = TRUE,
                     collapsible = TRUE,
                     width = 12,
                     tagList(
                       DT::renderDataTable({ BFabric$bfabricInstrumentEventsFiltered() })
                     ))

                 })
                 
                 
                 output$plot <- renderPlot({
                   shiny::req(dataFiltered())
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = paste0("GG-plotting ", title, "data ..."))
                   on.exit(progress$close())
                   
                   message(paste0("nrow dataFiltered(): ", nrow(dataFiltered())))
                   # dataFiltered() |> head() |> print()
                   
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
                     gp + ggplot2::geom_vline(xintercept = BFabric$bfabricInstrumentEventsFiltered()$time, linetype="dashed", 
                                              color = "red", size = 1) -> gp
                   }
                   gp
                 })
                 
                 return(data)
               })
}