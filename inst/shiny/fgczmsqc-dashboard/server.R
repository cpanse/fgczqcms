#R
## Christian Panse <cp@fgcz.ethz.ch> 2023-11-09
## server
library(shinylogs)

# requirements ===========
stopifnot(require(readr),
          require(reshape2),
          require(shinydashboard),
          require(lattice),
          require(rawDiag),
          require(rawrr))

source('helpers-ggplot2.R')
source('helpers-readr.R')
source('helpers.R')
source('module-autoQC01.R')
source('module-autoQC03.R')
source('module-bfabricInstrumentEvent.R')
source('module-config.R')
source('module-rawrr.R')

# define server logic ============
function(input, output, session) {
  track_usage(storage_mode = store_json(path = "logs-qc/"))
  
  #  reactives =============
  
  ## >> reactiveValues defined here ##################
  vals <- reactiveValues(timeRangeInSecs = 14 * 3600 * 24,
                         timeMin = (Sys.time() - (7 * 3600 * 24)),
                         timeMax = Sys.time(),
                         heightPx = 400,
                         addSmoothing = FALSE,
                         instrument = NULL,
                         useBFabric = FALSE,
                         peptide = c('LGGNEQVTR', 'YILAGVENSK', 'DGLDAASYYAPVR',
                                     'GTFIIDPAAVIR'))
  
  
  rootdir <- reactive({
    cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  
  ### initialize modules =============
  BFabric <- bfabricInstrumentEventServer("bfabric01", filterValues = vals)
  
  autoQC01 <- autoQC01Server("autoQC01",
                             filterValues = vals,
                             BFabric = BFabric,
                             inputfile = file.path(rootdir(),
                                                   'autoQC01-fit-apex-auc-fwhm.txt'))
  
  autoQC03DDA <- autoQC03Server("autoQC03-DDA",
                                filterValues = vals,
                                BFabric = BFabric,
                                inputfile = file.path(rootdir(), "comet.RData"),
                                readFUN = .readComet,
                                ggplot2FUN = .ggplotAutoQC03,
                                title = "DDA (Data-dependent acquisition)",
                                footer = "Running Comet (refer to https://github.com/UWPR/Comet) and utilizing UP000005640 FASTA as input generates the graphs.")
  
  autoQC03DIA <- autoQC03Server("autoQC03-DIA",
                                filterValues = vals,
                                BFabric = BFabric,
                                inputfile = file.path(rootdir(), "autoQC03-diann.txt"),
                                readFUN = .readDIANN,
                                ggplot2FUN = .ggplotAutoQC03,
                                title = "DIA (Data-independent acquisition)",
                                footer = "We compose the graphs by utilizing DIA-NN (check it out at https://github.com/vdemichev/DiaNN) and incorporating UP000005640 FASTA as input. We convert Orbitrap raw  files prior to mzML, employing ProteoWizard through  Docker and wine.")
  
  
  
  autoQC01alpha <- autoQC03Server("__autoQC01__",
                                  filterValues = vals,
                                   BFabric = BFabric,
                                   inputfile = file.path(rootdir(), 'autoQC01-fit-apex-auc-fwhm.txt'),
                                   readFUN = .readAutoQC01,
                                   ggplot2FUN = .ggplotAutoQC01,
                                   title = "autoQC01 - Biognosys iRT peptides runs",
                                   footer = "graphs pre-computed AUC | APEX | FWHM  values by utilizing rawrr (https://bioconductor.org/packages/rawrr/). Of note, only data acquired by Orbitraps are presented.")
  
  output$autoQC01 <- renderUI({
    autoQC01UI("autoQC01")
  })
  
  rawFileHeader <- reactive({
    shiny::req(input$file)
    
    rootdirraw() |>
      file.path(input$file) |>
      rawrr::readFileHeader()
  })
  
  iRTmz <- reactive({
    iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384,
               683.8282, 683.8541, 699.3388, 726.8361, 776.9301)
    
    names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                      "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                      "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                      "LFLQFGAQGSPFLK")
    
    iRTmz
  })

  ####### instruments -----------
  instruments <- reactive({
    names(.getInstruments())
  })
  
  
  rootdirraw <- reactive({
    cands <- c("/srv/www/htdocs/", "/scratch/DIAQC/qc/dump", "/Users/cp/Downloads/dump/")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  observeEvent({ input$addSmoothing }, {
    vals$addSmoothing <- input$addSmoothing
  })
  
  observeEvent({ input$heightPx }, {
    vals$heightPx <- as.integer(input$heightPx)
  })

  observeEvent({input$peptide}, {
    vals$peptide <- input$peptide
  })
  observeEvent({ input$useBfabric }, {
    vals$useBFabric <- input$useBfabric
  })
  
  observeEvent({ input$instrument }, {
    vals$instrument <- input$instrument
  })
  observeEvent({ input$timeRange }, {
    #shiny::req(input$timeRange)

    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    if (timeDiff < as.integer(input$timeRange) * 3600 * 24){
      vals$timeRangeInSecs <- as.integer(input$timeRange) * 3600 * 24
      
      if (vals$timeMin > Sys.time() - 21*24*3600){
        vals$timeMin <- Sys.time() - 21*24*3600
      }
      
      if (vals$timeMax < Sys.time()){
        vals$timeMax <- Sys.time()
      }
      
    }else{
      warning("timeDiff is higher than timeRange")
    }
  })
  
  observe({
    dt <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    
    msg <- paste0("input timeRange: ", vals$timeRangeInSecs,
                  " | timeMin: ", vals$timeMin,
                  " | timeMax: ", vals$timeMax,
                  " | timeDiff: ", dt)
    
    message(msg)
  })
  ## autoQC03 debounce --------
  ## TODO(cp): rename it
  autoQC03TimeRange_d <- reactive({
    shiny::req(input$cometTimeRange)
    input$cometTimeRange
  }) |> debounce(500)

  observeEvent({autoQC03TimeRange_d()}, {
    #shiny::req(input$cometTimeRange)
    
    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    
    if (timeDiff < vals$timeRangeInSecs){
      vals$timeMin <- autoQC03TimeRange_d()[1]
      vals$timeMax <- autoQC03TimeRange_d()[2]
      
    }else{
      warning("timeDiff is higher than timeRange")
    }
  })
  ## autoQC01 debounce --------
  autoQC01TimeRange_d <- reactive({
    shiny::req(input$autoQC01TimeRange)
    input$autoQC01TimeRange
  }) |> debounce(500)
  observeEvent({autoQC01TimeRange_d()}, {
    # shiny::req(input$autoQC01TimeRange)
   
    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    
    if (timeDiff < vals$timeRangeInSecs){
      vals$timeMin <- autoQC01TimeRange_d()[1]
      vals$timeMax <- autoQC01TimeRange_d()[2]
    }else{
      warning("timeDiff is higher than timeRange")
    }
  })
  #  renderUIs =============
  ## TimeSliders ---------------
  output$autoQC01TimeSlider <- renderUI({
    now <- Sys.time()
    if (vals$timeMin < (now - (vals$timeRangeInSecs))){
      vals$timeMin <- (now - (vals$timeRangeInSecs)) + 1
    }
    
    shinydashboard::box(
      sliderInput("autoQC01TimeRange", "Observation range:",
                  min = (now - (vals$timeRangeInSecs)),
                  max = now,
                  value = c(vals$timeMin, vals$timeMax),
                  timeFormat = "%F",
                  step = 3600 * 24,
                  width = "95%"),
      footer = "choose time range of autoQC01 files.", width = 12)
  })
  
  ## TODO(cp): rename to autoQC03TimeSlider
  output$cometTimeSlider <- renderUI({
    now <- Sys.time()
    
    if (vals$timeMin < (now - (vals$timeRangeInSecs))){
      vals$timeMin <- (now - (vals$timeRangeInSecs)) + 3600
    }
    
    shinydashboard::box(
      sliderInput("cometTimeRange", "Observation range:", 
                  min = (now - (vals$timeRangeInSecs)),
                  max = now,
                  value = c(vals$timeMin, vals$timeMax),
                  timeFormat = "%F",
                  step = 3600 * 24,
                  width = "95%"),
      footer = "choose time range of autoQC03 data.", width = 12)
  })
  
  
  output$diannTimeSlider <- renderUI({
    mintime <- min(diannLong()$time)
    now <- Sys.time()
    
    maxtime <- (1 + difftime(now, mintime, units = 'days') |>
                  round() |>
                  as.integer())
    
    sliderInput("diannDays", "Observation range in days:", min = 0,
                max = maxtime,
                value = c(0, min(28, maxtime)), width = "100%")
  })
  
  
#  output$instrument <- renderUI({
#   # fluidRow(
#      selectInput('instrument', 'instruments',
#                  instruments(),
#                  multiple = FALSE,
#                  selected = instruments()[1])
#   # )
#  })
  
  output$useBfabric <- renderUI({
    if (require(bfabricShiny)){
      L <- checkboxInput('useBfabric',
                                  'show B-Fabric Instrument Events',
                                  value = FALSE)
      return(L)
    }
    NULL
  })
  
  ##### TIC -------------------------------
  ticData <- reactive({
    shiny::req(input$ticfile)
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading TIC data ...")
    on.exit(progress$close())
    
    count <- 0
    
    rv <- file.path(rootdirraw(), input$ticfile) |>
      lapply(FUN = function(f){
        count <- count + 1
        progress$set(message = "Reading TIC from file", detail = basename(f), value = count/length(input$ticfile))
        rawrr::readChromatogram(f, type='tic') })
    
    rv
  })
  
  
  ## Composing summary data ============
  dataSummary <- reactive({
    shiny::req(autoQC03DDA(), autoQC03DIA(), autoQC01())
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = paste0("Composing summary data ..."))
    on.exit(progress$close())
    
    d1 <- data.frame(time = autoQC01()$time,
                     size = autoQC01()$size,
                     Instrument = autoQC01()$Instrument,
                     method = "Biognosys iRT (autoQC01)")
    
    autoQC03DDA() |> 
      subset(variable == "size") -> dd2
    
    d2 <- data.frame(time = dd2$time,
                     size = dd2$value,
                     Instrument = dd2$Instrument,
                     method = "DDA (comet)")
    
    autoQC03DIA() |> 
      subset(variable == "Size") -> dd3
    
    d3 <- data.frame(time = dd3$time,
                     size = dd3$value,
                     Instrument = dd3$Instrument,
                     method = "DIA (DIA-NN)")
    
    
    list(d1, d2, d3) |>
      Reduce(f = rbind)
  })
  
  superpose.symbol <- reactive({
    t <- trellis.par.get("superpose.symbol")
    cm <- c("#94C6FF", "#FFBBA9", "#76E3B8", "#FFD6AD", 
            "#BCE1FF", "#FFF691", "#FFC1E1")
    t$col <- cm
    t$fill <- cm
    
    return(t)
  })
  
  ## plotSummaryLCMSruns  --------
  output$plotSummaryLCMSruns  <- renderPlot({
    shiny::req(dataSummary())
    trellis.par.set("superpose.symbol", superpose.symbol())
    lattice::dotplot(Instrument ~ time | method,
                     groups = method,
                     data = dataSummary(),
                     alpha = 0.2,
                     cex = 2.4,
                     pch = 22,
                     layout = c(1, 3)
                     
                     
    )
  })
  
  
  output$plotSummaryCumsum  <- renderPlot({
    shiny::req(dataSummary())
    
    trellis.par.set("superpose.symbol", superpose.symbol())
    
    barchart(size/1024^3 ~ format(time, "%Y-%m") | method,
             data = dataSummary(),
             layout = c(1, 3),
             scales = list(x = list(rot = 45)))
  })
  
  ## printSummary --------
  output$summaryFrequency <- renderTable({
    (table(dataSummary()$method, dataSummary()$Instrument)) |> as.data.frame()
  }, rownames = FALSE, colnames = FALSE)
  
  ## plot TICs --------
  output$plotTIC <- renderPlot({
    shiny::req(ticData())
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Plotting TIC data ...")
    on.exit(progress$close())
    
    
    par(mfrow = c(length(input$ticfile), 1))
    
    ticData() |> lapply(function(x){plot(x)})
  }, height = function(){300 * length(input$ticfile)})
  
  #### render rawFileHeader  ----
  output$rawFileHeader <- renderPrint({
    shiny::req(rawFileHeader())
    capture.output(rawFileHeader())
  })
  
  output$cometBfabricInstrumentEventsOutput <- renderUI({
    shiny::req(input$useBfabric)
    shiny::req(BFabric$bfabricInstrumentEventsFiltered())
    
    DT::renderDataTable({ BFabric$bfabricInstrumentEventsFiltered() })
    
  })

  .determineLastEntry <- function(x){
    stopifnot('time' %in% colnames(x),
              'Instrument' %in% colnames(x))
    x[, c('time', 'Instrument')] |>
      split(f = x$Instrument) |>
      lapply(FUN = function(o){
        o.max <- max(o$time)
        idx <- which(o$time == o.max)
        o[idx[1], c('time', 'Instrument')]
      }) |>
      Reduce(f = rbind) -> x
    format(x$time, '%Y-%m-%d %H:%M')  -> x$time
    x[rev(order(x$time)), ]
  }

  output$lastEntryAutoQC01 <- renderTable({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "determine last entry", detail = "autoQC01")
    on.exit(progress$close())
    
    .determineLastEntry(autoQC01alpha())
  })
  
  output$lastEntryAutoQC03dda <- renderTable({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "determine last entry", detail = "autoQC03 dda")
    on.exit(progress$close())
    
    .determineLastEntry(autoQC03DDA())
  })
  
  output$lastEntryAutoQC03dia <- renderTable({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "determine last entry", detail = "autoQC03 dia")
    on.exit(progress$close())
    
    .determineLastEntry(autoQC03DIA())
  })
  
  # Values from cdata returned as text
  output$clientDataText <- renderText({
    cdata <- session$clientData
    cnames <- names(cdata)
    
    allvalues <- lapply(cnames, function(name) {
      paste(name, cdata[[name]], sep = " = ")
    })
    paste(allvalues, collapse = "\n")
  })
  
  output$sessionInfo <- renderPrint({
    sessionInfo() |> capture.output()
  })
}
