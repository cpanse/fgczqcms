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
                   
                   defaulVariablesIdx <-  dataVariables() %in% c('nConfidentProteins',
                                                                 'nConfidentPeptides', 'nMS2', 'Precursors.Identified',
                                                                 'Proteins.Identified', 'FWHM.RT')  |> which()
                   
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
                   message(paste0("nrow data(): ", nrow(data())))
                   filter <- data()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin < data()$time & data()$time < filterValues$timeMax 
                   #data()$variable %in% input$Variables 
                   
                   data()[filter,]
                 })
                 
                 output$ui <- renderUI({
                   tagList(
                     h2(title),
                     fluidRow(htmlOutput(ns("variable"))),
                     shinydashboard::box(
                       title = paste0(title, " plots"),
                       status = "primary",
                       solidHeader = TRUE,
                       collapsible = TRUE,
                       width = 12,
                       plotOutput(ns("plot"))
                     ),
                   )
                 })
                 
                 output$plot <- renderPlot({
                   shiny::req(dataFiltered())
                   
                   message(paste0("nrow dataFiltered(): ", nrow(dataFiltered())))
                   dataFiltered() |> head() |> print()
                   
                   #lattice::xyplot(value ~ time | Instrument * variable,
                   #                 subset = dataFiltered()$variable %in% input$variables,
                   #                data = dataFiltered())
                   dataFiltered() |> 
                     subset(variable %in% input$variables) |> 
                     ggplot2::ggplot(ggplot2::aes(time, value)) +
                     ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1)  -> gp
                   
                   if ("scanType" %in% names(dataFiltered())){
                     gp +
                       ggplot2::geom_point(ggplot2::aes(time, value, colour = factor(scanType)), alpha = 0.4) +
                       ggplot2::geom_line(ggplot2::aes(time, value, colour = factor(scanType)), alpha = 0.4) -> gp
                     
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
                 
               })
}