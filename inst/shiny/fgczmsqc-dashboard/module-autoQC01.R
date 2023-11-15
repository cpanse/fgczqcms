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
    htmlOutput(ns("plotSelection"), width = "100%")
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
                   message(paste0("autoQC01 module APEX long nrow: ", nrow(rv)))
                   message(paste0(input$peptides, collapse = ', '))
                   message(input$variables)
                   #rv[, c("time", "variable", "value", "Instrument", "peptide", input$variable)]
                   rv
                 })
                 
                 
                 output$plotSelection <- renderUI({
                   shiny::req(input$variables)
                   
                   L <- tagList()
                   if ("AUC" %in% input$variables){
                     L <- append(L, tagList(plotOutput(NS(id, "auc"))))
                   }
                   
                   if ("APEX" %in% input$variables){
                     L <- append(L, tagList(plotOutput(NS(id, "apex"))))
                   }
                   
                   if ("FWHM" %in% input$variables){
                     L <- append(L, tagList(plotOutput(NS(id, "fwhm"))))
                   }
                   
                   #tagList(plotOutput(NS(id, "auc")), plotOutput(NS(id, "fwhm")))
                  L
                 })

                 output$apex <- renderPlot({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 APEX module...")
                   on.exit(progress$close())
                   
                   lattice::xyplot(value ~ time | variable * Instrument,
                                   group = peptide,
                                   data = data(),
                                   type = 'b',
                                   auto.key = list(space = "right"),
                                   main = "MODULE",
                                   subset = variable == 'APEX'
                                   )
                 })
                 
                 output$auc <- renderPlot({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 AUC data ...")
                   on.exit(progress$close())
                   
                   lattice::xyplot(log(value, 10) ~ time | variable * Instrument,
                                   group = peptide,
                                   data = data(),
                                   type = 'b',
                                   #panel = .iqrPanel,
                                   auto.key = list(space = "right"),
                                   main = "MODULE",
                                   scales = list(y = list(relation = "free")),
                                   subset = variable == 'AUC',
                   )
                 })
                 
                 output$fwhm <- renderPlot({
                   
                   if (isFALSE("FWHM" %in% input$variables)){return(NULL)}
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 FWHM module...")
                   on.exit(progress$close())
                   
                   lattice::xyplot(value ~ time | variable * Instrument,
                                   group = peptide,
                                   data = data(),
                                   type = 'b',
                                   #panel = .iqrPanel,
                                   auto.key = list(space = "right"),
                                   main = "MODULE",
                                   #ylim = quantile(data()$value[data()$variable == 'FWHM'], c(0.05, 0.95), na.rm = TRUE),
                                   subset = variable == 'FWHM'
                   )
                 })
                 
                 
                 return(data)
               }
  )
}