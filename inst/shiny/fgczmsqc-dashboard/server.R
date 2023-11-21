#R
## Christian Panse <cp@fgcz.ethz.ch> 2023-11-09
## server

# requirements ===========
stopifnot(require(readr),
          require(reshape2),
          require(shinydashboard),
          require(lattice),
          require(rawrr))

source('helpers.R')
source('module-config.R')
source('module-autoQC01.R')
source('module-autoQC03.R')
source('module-bfabricInstrumentEvent.R')
source('module-rawrr.R')

# define server logic ============
function(input, output, session) {
  #  reactives =============
  
  ## >> reactiveValues defined here ##################
  vals <- reactiveValues(timeRangeInSecs = 14 * 3600 * 24,
                         timeMin = (Sys.time() - (7 * 3600 * 24)),
                         timeMax = Sys.time(),
                         instrument = NULL,
                         useBFabric = FALSE)
  
  
  rootdir <- reactive({
    cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  
  ### initialize modules =============
  BFabric <- bfabricInstrumentEventServer("bfabric01", filterValues = vals)
  
  autoQC01 <- autoQC01Server("autoQC01", filterValues = vals, BFabric = BFabric)
  
  autoQC03DDA <- autoQC03Server("autoQC03-DDA", filterValues = vals,
                                BFabric = BFabric,
                                inputfile = file.path(rootdir(), "comet.RData"),
                                readFUN = .readComet, title="DDA")
  
  autoQC03DIA <- autoQC03Server("autoQC03-DIA", filterValues = vals,
                                BFabric = BFabric,
                                inputfile = file.path(rootdir(), "autoQC03-diann.txt"),
                                readFUN = .readDIANN, title="DIA")
  
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
  
  #### iRTprofileRawDDA --------------
  iRTprofileRawDDA <- reactive({
    shiny::req(input$file)
    shiny::req(input$Ms1ppmError)
    shiny::req(iRTmz())
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading iRT peptide profiles ...")
    on.exit(progress$close())
    
    file.path(rootdirraw(), input$file) |>
      rawrr::readChromatogram(mass = iRTmz(),
                            tol = as.integer(input$Ms1ppmError),
                            type = "xic",
                            filter = "ms")
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
  
  files <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Composing file list ...")
    on.exit(progress$close())
    
    c( #(autoQC01Data()$filename |> unique()),
      (diannData()$File.Name |> unique()),
      (cometData()$File.Name |> unique())) |>
      unique() |>
      lapply(function(f){
        if(file.exists(file.path(rootdirraw(), f))){return(f)}
        else{
          # msg <- paste0("Can not find file ", f)
          # message(msg)
        }
        NULL
      }) |>
      unlist() 
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
  
  observeEvent({input$cometTimeRange}, {
    #shiny::req(input$cometTimeRange)
    
    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    
    if (timeDiff < vals$timeRangeInSecs){
      vals$timeMin <- input$cometTimeRange[1]
      vals$timeMax <- input$cometTimeRange[2]
      
    }else{
      warning("timeDiff is higher than timeRange")
    }
  })
  
  observeEvent({input$autoQC01TimeRange}, {
   # shiny::req(input$autoQC01TimeRange)
    
    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    
    if (timeDiff < vals$timeRangeInSecs){
      vals$timeMin <- input$autoQC01TimeRange[1]
      vals$timeMax <- input$autoQC01TimeRange[2]
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
                  width = "95%"), footer = "choose time range of autoQC01 files.", width = 12)
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
  
  
  output$instrument <- renderUI({
    fluidRow(selectInput('instrument', 'Instruments',
                              instruments(),
                              multiple = FALSE,
                              selected = instruments()[1]))
  })
  
  output$useBfabric <- renderUI({
    if (require(bfabricShiny)){
      L <- fluidRow(checkboxInput('useBfabric',
                                 'show B-Fabric Instrument Events',
                                 value = FALSE))
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
  
  ## plotSummary --------
  output$plotSummary  <- renderPlot({
    shiny::req(dataSummary())
   
    trellis.par.set("superpose.symbol", superpose.symbol())
    
    lattice::dotplot(Instrument ~ time | method,
                     groups = method,
                     data = dataSummary(),
                     alpha = 0.2,
                     cex = 2.4,
                     pch = 22,
                     layout = c(1, 3),
                     main = "instrument events")
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
  output$summary <- renderPrint({
    capture.output( table(dataSummary()$method, dataSummary()$Instrument))
  })

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
  
  #### render sessionInfo ----
  output$sessionInfo <- renderPrint({
    (.getInstruments()[input$instrument] |> unlist() |> paste(collapse = ";") |> message())
    
    capture.output(sessionInfo())
  })
}
