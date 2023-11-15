#R


autoQC01UI <- function(id){
  ns <- NS(id)
  
  p <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
         "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
         "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
         "LFLQFGAQGSPFLK")
  
  v <- c("APEX", "AUC", "FWHM")
  
  tagList(
    selectInput(ns('peptides'), "peptides", multiple = TRUE, choices = p, selected = p[c(1, 6, 11)]),
    selectInput(ns('variables'), "variables", multiple = TRUE, choices = v, selected = "AUC"),
    tableOutput(NS(id, "nearAuc")),
    plotOutput(NS(id, "auc"), click = NS(id, "plot_click"), height = 600),
  )
}

#' autoQC01Server shiny module
#'
#' @param id 
#' @param ii instrument TODO(cp): should be added to filterValues by using observeEvent
#' @param filterValues time and instrument information
#'
#' @return data.frame 
#' @export
autoQC01Server <- function(id, ii,  filterValues){
  moduleServer(id,
               function(input, output, session) {
                 
                 instrument <- reactive({ii()})
                 
                 autoQC01APEXwide <- reactive({
                   fn <- "autoQC01-fit-apex-auc-fwhm.txt"
                   progress <- shiny::Progress$new(session = session)
                   msg <- paste0("Reading ", fn)
                   progress$set(message = msg, detail = "APEX | AUC | FWHM")
                   on.exit(progress$close())
                   
                   config <- configServer("config1")
                   
                   config$rootdir |>
                     file.path(fn) |>
                     readr::read_delim(
                       delim = ";",
                       escape_double = FALSE,           
                       col_names =c('filename', 'time', 'peptide', 'AUC', 'APEX', 'FWHM'),
                       col_types = readr::cols(time = col_datetime(format = "%s")),
                       trim_ws = TRUE) -> rv
                   
                   ## TODO(cp): perform the ordering before!
                   rv[order(rv$time), ] -> rv
                   message(paste0("autoQC01 module APEX wide nrow: ", nrow(rv)))
                   
                   rv$Instrument <- NA
                   rv |> .assignInstrument(coln = 'filename')
                 })
                 
                 data <- reactive({
                   shiny::req(autoQC01APEXwide())
                   shiny::req(input$peptides)
                   #shiny::req(input$variables)
                   
                   autoQC01APEXFilter <- autoQC01APEXwide()$Instrument %in% instrument() &
                     filterValues$timeMin < autoQC01APEXwide()$time & autoQC01APEXwide()$time < filterValues$timeMax &
                     autoQC01APEXwide()$peptide %in% input$peptides
                 
                   autoQC01APEXwide()[autoQC01APEXFilter, ] |>
                     reshape2::melt(id.vars = c("filename", "time", "Instrument", "peptide")) -> rv
                   
                   
                   rv$value[rv$variable == "AUC"] <- log(rv$value[rv$variable == "AUC"], 10) 
                   
                   message(paste0("autoQC01 module APEX long nrow: ", nrow(rv)))
                   message(paste0(input$peptides, collapse = ', '))
                   message(input$variables)
                   #rv[, c("time", "variable", "value", "Instrument", "peptide", input$variable)]
                   rv
                 })
                 
                 
                 output$nearAuc <- renderTable({
                   req(input$plot_click)
                   nearPoints(data(), input$plot_click, xvar = "time", yvar = "value")
                   # [data()$variable == 'AUC', ]
                 })
                 
                 output$auc <- renderPlot({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 data ...")
                   on.exit(progress$close())
                   
                   subset(data(), variable %in% input$variables) |> 
                   ggplot2::ggplot(ggplot2::aes(time, value)) +
                     ggplot2::geom_point(ggplot2::aes(color = peptide), alpha = 0.4) +
                     ggplot2::geom_line(ggplot2::aes(group = peptide, color = peptide), alpha = 0.4) +
                     #ggplot2::scale_y_log10() +
                     ggplot2::facet_wrap(. ~ variable, scales="free_y", ncol = 1) 
                 }, res = 96)

                #  return(data)
               }
  )
}