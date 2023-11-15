#R

configServer <- function(id){
  moduleServer(id,
               function(input, output, session) {
                 rootdirraw <- reactive({
                   cands <- c("/srv/www/htdocs/", "/scratch/DIAQC/qc/dump", "/Users/cp/Downloads/dump/")
                   for (d in cands){
                     if (dir.exists(d))return(d)
                   }
                   NULL
                 })
                 
                 rootdir <- reactive({
                   cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
                   for (d in cands){
                     if (dir.exists(d))return(d)
                   }
                   NULL
                 })
                 return(list(rootdir = rootdir(),
                             rootdirraw = rootdirraw()))
               }
  )
}

autoQC01UI <- function(id){
  ns <- NS(id)
  tagList(
    plotOutput(ns("apex")),
    plotOutput(ns("auc"))
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
                       col_types = cols(time = col_datetime(format = "%s")),
                       trim_ws = TRUE) -> rv
                   
                   rv[order(rv$time), ] -> rv
                   message(paste0("autoQC01 APEX wide nrow: ", nrow(rv)))
                   
                   rv$Instrument <- NA
                   rv |> .assignInstrument(coln = 'filename')
                 })
                 
                 data <- reactive({
                   shiny::req(autoQC01APEXwide())
                   
                   autoQC01APEXFilter <- autoQC01APEXwide()$Instrument %in% instrument() &
                     filterValues$timeMin < autoQC01APEXwide()$time & autoQC01APEXwide()$time < filterValues$timeMax
                 
                   autoQC01APEXwide()[autoQC01APEXFilter, ] |>
                     reshape2::melt(id.vars = c("filename", "time", "Instrument", "peptide")) 
                 })
                 
                 output$apex <- renderPlot({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 APEX module...")
                   on.exit(progress$close())
                   
                   lattice::xyplot(value ~ time | variable * Instrument,
                                   group = peptide,
                                   data = data(),
                                   type = 'b',
                                   auto.key = list(space = "bottom"),
                                   main = "MODULE",
                                   subset = variable == 'APEX'
                                   )
                 })
                 
                 output$auc <- renderPlot({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 AUC module...")
                   on.exit(progress$close())
                   
                   lattice::xyplot(value ~ time | variable * Instrument,
                                   group = peptide,
                                   data = data(),
                                   # type = 'b',
                                   panel = .iqrPanel,
                                   auto.key = list(space = "bottom"),
                                   main = "MODULE",
                                   ylim = quantile(data()$value, c(0.05, 0.95)),
                                   subset = variable == 'AUC'
                   )
                 })
                 
                 
                 return(data)
               }
  )
}