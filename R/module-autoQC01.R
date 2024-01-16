#R

#' @importFrom shinydashboard box 
#' @importFrom shiny fluidRow
#' @export
autoQC01UI <- function(id){
  ns <- NS(id)

  tagList(
    shinydashboard::box(title = "InstrumentEvents",
                        fluidRow(htmlOutput(ns("instrumentEventsOutput"))),
                        footer = "once enabled it shows the instrument events.",
                        status = "primary",
                        solidHeader = TRUE,
                        width = 12),
    # TODO(cp): make that reactive and enable status mtime.
    htmlOutput(NS(id, "AucApexFwhmUI")),
    htmlOutput(NS(id, "rawrrEnableUI")),
    #rawrrUI(NS(id, "rawrr01")),
    shinydashboard::box(title = "autoQC01 iRT peptide fit",
                        fluidRow(
                          column(3, offset = 0, htmlOutput(NS(id, "autoQC01Variable")))
                        ),
                        fluidRow(
                          column(12, offset = 0, plotOutput(NS(id, "autoQC01Plot")))
                        ),
                        status = "primary",
                        solidHeader = TRUE,
                        width = 12)
  )
}

#' autoQC01Server shiny module
#'
#' @param id 
#' @param filterValues time and instrument information
#'
#' @return data.frame 
#' @export
autoQC01Server <- function(id, filterValues, BFabric, inputfile){
  moduleServer(id,
               function(input, output, session) {
                 ns <- NS(id)
                 rootdirraw <- reactive({
                   config <- configServer("config2")
                   config$rootdirraw
                 })
                 
                 vals <- reactiveValues(fn = NA,
                                        mZ = .iRTmz(),
                                        hover = NA)
                 
                 rawrrServer("rawrr01", vals) #[names(.iRTmz()) %in% input$peptides])
                 ## observeEvent hover ================
                 observeEvent(input$hoverInfo, {
                   if(!is.null(input$hoverInfo)){
                     np <- nearPoints(data(), input$hoverInfo, xvar = "time", yvar = "value")
                     
                     cat("Hover (throttled):\n")
                     str(np$filename[1])
                     vals$hover <- np[1, ]
                   }
                 })
                 output$AucApexFwhmUI <- renderUI({
                   p <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                          "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                          "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                          "LFLQFGAQGSPFLK")
                   
                   v <- c("APEX", "AUC", "AUC.lg2", "FWHM")
                   
                   shinydashboard::box(title = paste0("AUC | APEX | FWHM  plots - mtime: ",
                                                      file.mtime(inputfile) |> strftime("%a %F %T")),
                     fluidRow(
                       column(9, offset = 0,
                              tagList(
                                plotOutput(NS(id, "auc"),
                                           click = NS(id, "plot_click"),
                                           hover = hoverOpts(NS(id, "hoverInfo"), delay = 200, nullOutside = TRUE),
                                           height = 600))),
                       column(3, offset = 0,
                              tagList(
                                selectInput(ns('peptides'), "peptides", multiple = TRUE, choices = p, selected = p[c(1, 6, 11)]),
                                selectInput(ns('variables'), "variables", multiple = TRUE, choices = v, selected = c("APEX", "AUC.lg2")),
                                tableOutput(NS(id, "hoverInfo")),
                                #sliderInput(ns('AUC.lg2.cutoff'), "AUC.lg2.cutoff", min = 9, max = 40, value = c(10, 40) , step = 1),
                              )
                       )
                     ), footer = "click to analyze selected raw file below.",
                     status = .status(inputfile),
                     solidHeader = TRUE,
                     width = 12
                   )
                 })
                 
                 output$rawrrEnableUI <- renderUI({
                   shiny::req(vals$fn)
                   if (file.exists(vals$fn)){
                     rawrrUI(NS(id, "rawrr01"))
                   }
                 })
                 ## autoQC01 ---------
                 autoQC01wide <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Reading autoQC01 lm data ...")
                   on.exit(progress$close())
                   config <- configServer("config1")
                   rv <- config$rootdir |>
                     file.path("autoQC01.csv") |>
                     readr::read_delim(
                       delim = ";",
                       escape_double = FALSE,           
                       col_types = readr::cols(time = col_datetime(format = "%s"),                       
                                               size = col_integer(),                       
                                               n = col_integer()),
                       trim_ws = TRUE)

                   rv$Instrument <- NA
                   rv |> .assignInstrument(coln = 'filename')
                 })
                 
                 autoQC01Long <- reactive({
                   shiny::req(filterValues$instrument,
                              autoQC01wide())
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Reshaping autoQC01 lm data ...")
                   on.exit(progress$close())
                   
                   autoQC01wideFilter <- autoQC01wide()$Instrument %in% filterValues$instrument & 
                     filterValues$timeMin < autoQC01wide()$time & autoQC01wide()$time < filterValues$timeMax  
                   
                   
                   rv <- autoQC01wide()[autoQC01wideFilter, ] |>
                     reshape2::melt(id.vars = c("md5", "filename", "time", "Instrument"))
                 })
                 
                 dataPeptideFit <- reactive({
                   shiny::req(autoQC01Long())
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Filtering autoQC01 data ...")
                   on.exit(progress$close())
                   
                   now <- Sys.time()
                   
                   autoQC01LongFilter <- autoQC01Long()$variable %in% input$autoQC01Variables 
                   
                   autoQC01Long()[autoQC01LongFilter, ]
                 })
                 
                 
                 autoQC01Variables <- reactive({
                   autoQC01Long()$variable |> unique()
                 })
                 
                 output$autoQC01Variable <- renderUI({
                   selectInput(NS(id,'autoQC01Variables'), 'Variables',
                               autoQC01Variables(),
                               multiple = TRUE,
                               selected = c('slope', 'r.squared', 'intercept'))
                 })
                 
                 #### autoQC01 lattice::xyplot -----------------
                 output$autoQC01Plot <- renderPlot({
                   shiny::req(dataPeptideFit())
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting autoQC01 lm ...")
                   on.exit(progress$close())
                   
                   if (filterValues$useBFabric){
                     .lattice(dataPeptideFit(),
                              useBfabric = filterValues$useBFabric,
                              bfabricInstrumentEvents = BFabric$bfabricInstrumentEventsFiltered()$time)
                   }else{
                     .lattice(dataPeptideFit())
                   }
                   
                 })
                 
                 ############################
                 
                 ## data AUC|APEX|FWHM ==================
                 autoQC01APEXwide <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   msg <- paste0("Reading ", inputfile |> basename(), " ...")
                   progress$set(message = msg, detail = "APEX | AUC | FWHM", value = 1/5)
                   on.exit(progress$close())
                   
                   inputfile |>
                     readr::read_delim(
                       delim = ";",
                       escape_double = FALSE,           
                       col_names =c('filename', 'time', 'peptide', 'AUC', 'APEX', 'FWHM'),
                       col_types = readr::cols(time = col_datetime(format = "%s")),
                       trim_ws = TRUE) -> rv
                   
                   rv$AUC.lg2 <- log(rv$AUC, 2)
                   
                   ## TODO(cp): perform the ordering before!
                   progress$set(detail = "ordering ...", value = 2/5)                  
                   rv[order(rv$time), ] -> rv
                   message(paste0("autoQC01 module APEX wide nrow: ", nrow(rv)))
                   
                   progress$set(detail = "skip checking if file.exists", value = 3/5)
                   #rv$file.exists <- file.path(config$rootdirraw, rv$filename) |> sapply(FUN = file.exists)
                   rv$file.exists <- TRUE
                   
                   progress$set(detail = "assigning instrumentss", value=  4/5)                   
                   rv$Instrument <- NA
                   rv |> .assignInstrument(coln = 'filename')
                 })

                 peptides_d <- reactive({
                   shiny::req(input$peptides)
                   input$peptides
                 }) |> debounce(2000)

                 data <- reactive({
                   # shiny::req(autoQC01APEXwide(), peptides_d())
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Reshaping autoQC01", detail = "Filtering", value = 1/3)
                   on.exit(progress$close())
                   
                   autoQC01APEXFilter <- autoQC01APEXwide()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin < autoQC01APEXwide()$time &
                     autoQC01APEXwide()$time < filterValues$timeMax &
                     autoQC01APEXwide()$peptide %in% peptides_d() 
                   
                   #&
                   #   as.integer(input$AUC.lg2.cutoff[1]) <= autoQC01APEXwide()$AUC.lg2 &
                   #  autoQC01APEXwide()$AUC.lg2 <= as.integer(input$AUC.lg2.cutoff[2])
                   
                   
                   
                   progress$set(detail = "reshape2::melt", value = 2/3)
                   autoQC01APEXwide()[autoQC01APEXFilter, ] |>
                     reshape2::melt(id.vars = c("filename", "file.exists", "time", "Instrument", "peptide")) -> rv
                   
                   
                   #progress$set(detail = "log AUC", value = 3/3)
                   #rv$value[rv$variable == "AUC"] <- log(rv$value[rv$variable == "AUC"], 2) 
                   
                   message(paste0("autoQC01 module APEX long nrow: ", nrow(rv)))
                   message(paste0(input$peptides, collapse = ', '))
                   message(input$variables)
                   #rv[, c("time", "variable", "value", "Instrument", "peptide", input$variable)]
                   rv
                 })
                 
                 ## -----------dblclick table ----------------
                 output$dblclick_info <- renderPrint({
                   #req(input$dblplot_click)
                   
                   cat("input$plot_dblclick:\n")
                   str(input$plot_dblclick)
                   np <- nearPoints(data(), input$dblplot_click, xvar = "time", yvar = "value")
                   message(paste0("plot_dblclick: ", np$filename, collapse = " | "))
                   #str(paste0(np$filename, collapse = " | "))
                 })
                 
                 output$hoverInfo <- renderUI({
                   shiny::req(vals$hover)
                   L <- tagList()
                   if (is.na(vals$hover$filename)){
                     L <- HTML("Hover over a point to see file information. <br>
                               Click on a point to trace peptides and gather raw file header information.")
                   }else{
                     L <- renderTable(t(vals$hover[, c('filename', 'time')]), rownames = T, colnames = F)
                   }
                   
                   shinydashboard::box(title = 'file information',
                                       L,
                                       status = .fileStatus(rootdirraw(), vals$hover$filename),
                                       solidHeader = TRUE,
                                       collapsible = FALSE,
                                       width = 12)
                 })
                 ## -----------click table ----------------
                 output$nearAuc <- renderTable({
                   req(input$hoverInfo)
                   np <- nearPoints(data(), input$hoverInfo, xvar = "time", yvar = "value")
                   message(paste0("hoverInfo: ", np$filename, collapse = " | "))
                   t(np[1,])
                 }, rownames = T, colnames = F)
                 
                 ## observeEvent =================
                 observeEvent(input$plot_click, {
                   np <- nearPoints(data(), input$plot_click, xvar = "time", yvar = "value")
                   
                   message(paste0("plot_click: ", np$filename, collapse = " | "))
                   
                   config <- configServer("config2")
                   vals$fn <- file.path(config$rootdirraw, np$filename[1])
                 })

                 observeEvent(peptides_d(), {
                   vals$mZ <- .iRTmz()[names(.iRTmz()) %in% peptides_d()] 
                 })

                 ## ----------- ggplot2::ggplot -----------------
                 output$auc <- renderPlot({
                   shiny::req(data(), peptides_d())
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Plotting peptide data ...",
                                detail  = paste0(input$variables, collapse = ' |'))
                   on.exit(progress$close())
                   if (nrow(data()) > 0){
                     data() |>
                       subset(variable %in% input$variables) |>
                       ggplot2::ggplot(ggplot2::aes(time, value)) +
                       ggplot2::geom_point(ggplot2::aes(color = peptide), alpha = 0.4) +
                       ggplot2::geom_line(ggplot2::aes(group = peptide, color = peptide), alpha = 0.4) +
                       ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1)  -> gp
                     if (filterValues$useBFabric){
                       gp + ggplot2::geom_vline(xintercept = BFabric$bfabricInstrumentEventsFiltered()$time, linetype="dashed",
                                                color = "red", size = 1) -> gp
                     }
                     gp
                   }else{.missing()}
                 }, res = 96)

                 output$instrumentEventsOutput <- renderUI({
                   if (filterValues$useBFabric){
                     shiny::req(BFabric$bfabricInstrumentEventsFiltered())
                     DT::renderDataTable({ BFabric$bfabricInstrumentEventsFiltered() })
                   }
                 })
                 return(autoQC01wide)
               }
  )
}