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
source('module-bfabricInstrumentEvent.R')

# define server logic ============
function(input, output, session) {
  #  reactives =============
  
  ## >> reactiveValues defined here ##################
  vals <- reactiveValues(timeRangeInSecs = 14 * 3600 * 24,
                         timeMin = (Sys.time() - (7 * 3600 * 24)),
                         timeMax = Sys.time(),
                         instrument = NULL)
  
  
  rootdir <- reactive({
    cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  
  autoQC01 <- autoQC01Server("autoQC01", filterValues = vals)
  
  BFabric <- bfabricInstrumentEventServer("bfabric01", filterValues = vals)
  
  output$autoQC01 <- renderUI({
    autoQC01UI("autoQC01")
  })
  
  output$bfabric <- renderUI({
    #shiny::req(input$useBFabric)
    shiny::req(BFabric$bfabricInstrumentEventsFiltered())
    bfabricInstrumentEventUI("bfabric01")
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
  
  #### iRTprofileRawDIA --------------
  iRTprofileRawDIA <- reactive({
    shiny::req(input$file)
    shiny::req(input$iRTpeptide)
    shiny::req(input$scanType)
    shiny::req(input$Ms2ppmError)
    
    progress <- shiny::Progress$new(session = session)
    msg <- paste0("Reading Ms2 profiles for ", input$iRTpeptide, " ...")
    progress$set(message = msg)
    on.exit(progress$close())
    
    yIonSeries <- (input$iRTpeptide |> 
      protViz::fragmentIon())[[1]]['y'] |>
        unlist() 
    
    print(yIonSeries[seq(1, nchar(input$iRTpeptide) - 1)])
    
    file.path(rootdirraw(), input$file) |>
     rawrr::readChromatogram(mass = yIonSeries[seq(1, nchar(input$iRTpeptide) - 1)],
                              tol = as.integer(input$Ms2ppmError),
                              type = "xic",
                              filter = input$scanType) 
  })
  
  
  scanType <- reactive({
    shiny::req(input$file)
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading scanTypes of raw file ...")
    on.exit(progress$close())
    
    file.path(rootdirraw(), input$file) |>
      rawrr::readIndex() -> idx
    
    idx$scanType |> 
      unique() |>
      sort() 
  })
  
  rtFittedAPEX <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Fitting iRT APEXs ...")
    on.exit(progress$close())
    
    iRTprofileRawDDA() |>
      rawrr:::pickPeak.rawrrChromatogram() |>
      .fitPeak.rawrrChromatogram(delta = 0.5, n = 400)
  })
  
  
  ###### comet --------
  cometWide <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading comet data ...")
    on.exit(progress$close())
    e <- new.env()
    fn <- rootdir() |>
      file.path("comet.RData") |>
      load(envir=e)
    
    e$comet$assignmentRate <- round (100 * e$comet$nConfidentPSM / e$comet$nPSM)
    
    ## rename columns
    colnames(e$comet)[colnames(e$comet) == "filename.y"] <- "File.Name"
    colnames(e$comet)[colnames(e$comet) == "Time"] <- "time"
    colnames(e$comet)[colnames(e$comet) == "instrument"] <- "Instrument"
    colnames(e$comet)[colnames(e$comet) == "nConfidentPeptide"] <- "nConfidentPeptides"
    
    cc <- c('md5', 'File.Name', 'time', 'Instrument',  'scanType', 'psmFdrCutoff',
            'nDecoyPSM', 'nConfidentPSM', 'nDecoyPeptide', 'nConfidentPeptides',
            'nDecoyProteins', 'nConfidentProteins', 'fdrPSM', 'fdrPeptide',
            'fdrProtein',  'size',   'nMS2', 'TIC')
    e$comet$time <- as.POSIXct(e$comet$time )
    e$comet[, cc]
  })
  
  cometLong <- reactive({
    shiny::req(cometWide())
    
    
    rv <- cometWide() |>
      reshape2::melt(id.vars = c("md5", "File.Name", "time", "Instrument",
                                 "scanType"))
    
    rv[grepl(input$regex, rv$File.Name), ]
  })
  
  cometData <- reactive({
    shiny::req(input$instrument)
    shiny::req(cometLong())
   
    now <- Sys.time()
    
    
    msg <- paste0("cometTimes:", sum(vals$timeMin < cometLong()$time), " | ", sum(cometLong()$time < vals$timeMax))
    message(msg) 
    
    cometFilter <- cometLong()$Instrument %in% input$instrument &
      cometLong()$variable %in% input$cometVariables &
      vals$timeMin < cometLong()$time & cometLong()$time < vals$timeMax
  
     message(paste0("comet filter length: ", sum(cometFilter)))
    
    cometLong()[cometFilter, ]
  }) 
  

  
  cometVariables <- reactive({
    cometLong()$variable |> unique()
  })
  
 
  
  
  ########### diann -----------
  diannWide <- reactive({
    S <- rootdir() |> 
      file.path("output.txt") |>
      readr::read_delim(
        delim = ";",
        escape_double = FALSE,
        col_types = cols(Time = col_datetime(format = "%s"),
                         Size = col_integer(),
                         Precursors.Identified = col_integer(),
                         Proteins.Identified = col_integer()), 
        trim_ws = TRUE) 
    
    S <- S[grepl(input$regex, S$File.Name), ]
    
    colnames(S)[colnames(S) == "Time"] <- "time"
    
    S$Instrument <- NA
    S |> .assignInstrument()
  })
  
  diannLong <- reactive({
    shiny::req(diannWide())
    diannWide() |> reshape2::melt(id.vars = c("Md5", "File.Name", "time", "Instrument"))
  })
  
  
  ####### instruments -----------
  instruments <- reactive({
    names(.getInstruments())
  })
  
  ####### variables -----------
  diannVariables <- reactive({
    shiny::req(diannLong())
    diannLong()$variable |> unique()
  })
   
  # TODO(cp): rename to diannData <-
  diannData <- reactive({
    # shiny::req(input$variables)
    shiny::req(input$instrument)
    shiny::req(diannLong())
    
    now <- Sys.time()
    diannLong()[diannLong()$Instrument %in% input$instrument &
             diannLong()$variable %in% input$diannVariables &
             input$diannDays[1] <=  difftime(now, diannLong()$time, units = "days") &
             difftime(now, diannLong()$time, units = "days") < input$diannDays[2], ]
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
  

  observeEvent({ input$instrument }, {
    vals$instrument <- input$instrument
  })
  observeEvent({ input$timeRange }, {
    #shiny::req(input$timeRange)
    
    timeDiff <- difftime(vals$timeMax, vals$timeMin, units = "secs")
    if (timeDiff < as.integer(input$timeRange) * 3600 * 24){
      vals$timeRangeInSecs <- as.integer(input$timeRange) * 3600 * 24
      
      if (vals$timeMin > Sys.time() - 7*24*3600){
        vals$timeMin <- Sys.time() - 7*24*3600
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
    sliderInput("autoQC01TimeRange", "Observation range:",
                min = (now - (vals$timeRangeInSecs)),
                max = now,
                value = c(vals$timeMin, vals$timeMax),
                timeFormat = "%F",
                step = 3600 * 24,
                width = "95%")
  })
  
  output$cometTimeSlider <- renderUI({
    now <- Sys.time()
    
    if (vals$timeMin < (now - (vals$timeRangeInSecs))){
      vals$timeMin <- (now - (vals$timeRangeInSecs)) + 3600
    }
       
    sliderInput("cometTimeRange", "Observation range:", 
                min = (now - (vals$timeRangeInSecs)),
                max = now,
                value = c(vals$timeMin, vals$timeMax),
                timeFormat = "%F",
                step = 3600 *24,
                width = "95%")
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


  
  output$scanTypeUI <- renderUI({
    shiny::req(scanType())
    list(radioButtons("Ms2ppmError", "Ms2 ppmError",
                      choices = c(10, 20, 30, 50, 100),
                      selected = 30,
                      inline = TRUE,
                      width = NULL),
         selectInput('scanType', 'scanType',
                      scanType(),
                      multiple = FALSE,
                      selected = scanType()[1]),
          selectInput('iRTpeptide', 'iRTpeptide fragment ions',
                      names(iRTmz()),
                      multiple = FALSE,
                      selected = names(iRTmz())[1]),
          sliderInput("rtSlider", "rtSlider", min = 0, max = 120,
                      #max = 1 + ((iRTprofileRawDDA()[[1]][['times']]) |> max() |> round()),
                      value = c(26, 29), width = "100%"))
  })
  
  output$diannVariable <- renderUI({
    defaulVariables <- c('Precursors.Identified', 'Proteins.Identified', 'FWHM.RT')
    selectInput('diannVariables', 'Variables',
                diannVariables(),
                multiple = TRUE,
                selected = defaulVariables)
  })
  
  output$ticFileInput <- renderUI({
    shiny::req(files())
   
      L <- tagList(selectInput('ticfile', 'TIC file candidates',
                            files(),
                            multiple = TRUE,
                            selected = files()[1]),
                hr())
    
    return(L)
  })
  
  output$fileInput <- renderUI({
    shiny::req(files())
    
    L <- tagList(fluidRow(box(selectInput('file', 'Files',
                                          files(),
                                          multiple = FALSE,
                                          selected = files()[1]),width="100%"), width = "100%"),
                 fluidRow(checkboxInput("showRawFileHeader", "show raw File Header", value = FALSE)),
                 fluidRow(checkboxInput("showIrtMS1Profile", "show iRT Ms Profile", value = FALSE)),
                 fluidRow(checkboxInput("showIrtMS2Profile", "show iRT Ms2 Profile" , value = FALSE)))
    return(L)
    
  })
  
  ### render fileOutput -----------
  output$fileOutput <- renderUI({
    shiny::req(files())
    L <- tagList()
    
    if (input$showRawFileHeader){
      L <- tagList(fluidRow(
        fluidRow(box(verbatimTextOutput("rawFileHeader"), width = 800))
      ))
    }
    
    if (input$showIrtMS1Profile){
      L <- append(L, tagList(
        fluidRow(radioButtons("Ms1ppmError", "Ms1 ppmError",
                              choices = c(5, 10, 20, 30, 50, 100),
                              selected = 10,
                              inline = TRUE,
                              width = NULL)),
        fluidRow(box(plotOutput("plotiRTDDAChromatograms"), width = "100%")),
        fluidRow(box(plotOutput("plotDDAiRTfits"))),
        fluidRow(box(plotOutput("plotDDAiRTprofiles"), width = "100%"))
      ))
    }
    
    if (input$showIrtMS2Profile){ 
      L <- append(L, tagList(
        fluidRow(shiny::htmlOutput("scanTypeUI")),
        fluidRow(box(plotOutput("plotDIAiRTprofiles"), width = "100%"))
      ))
    }
    
    return(L)
  })
  #### renderPlots DIA-NN -------------
  output$tableDIANN <-  DT::renderDataTable({ diannData() })
  output$tableComet <-  DT::renderDataTable({ cometData()  })
  
  output$plotiRTDDAChromatograms <- renderPlot({
    shiny::req(input$showIrtMS1Profile)
    shiny::req(iRTprofileRawDDA())
    
    if(input$showIrtMS1Profile){
      iRTprofileRawDDA() |>
        plot(main = input$file)
    }
  })

  #### plotDIAiRTprofiles  ------------
  output$plotDIAiRTprofiles <- renderPlot({
    shiny::req(iRTprofileRawDIA())
    
    message("calling iRTprofileRawDIA()  ...")
    iRTprofileRawDIA()  |>
      .plotChromatogramSet(xlim=input$rtSlider)
  })
  
  .fwhm <- function(x, y){
    ymax <- max(y)
    halfmax <- ymax / 2
    
    idxMax <- which(y == ymax)[1]
    
    for (i in 1:length(y)){
      if (y[i] >= halfmax){
        break
      }
    }
    
    fwhm <- 2 * (x[idxMax] - x[i])
    
    return(list(x1= x[i], y1 = y[i], 
                idxMax = idxMax,
                fwhm = fwhm))
  }
  
  output$plotDDAiRTprofiles <- renderPlot({
    shiny::req(input$showIrtMS1Profile)
    shiny::req(rtFittedAPEX())
    
    par(mfrow = c(3, 4), mar = c(6, 4, 4, 1))
    rtFittedAPEX <- rtFittedAPEX() |>
      lapply(function(x){
        AUC <- sum(diff(x$xx) * (head(x$yp, -1) + tail(x$yp,  -1))) / 2
        APEX <- x$xx[which.max(x$yp)[1]]
        # TODO(cp): determine Full width at half maximum (FWHM)
        FWHM <- .fwhm(x$xx, x$yp)
        r.squared <- x$r.squared[1]
        
        if (FWHM$fwhm > 0){
          plot(x$times, x$intensities,
               type='p',
               sub = sprintf("AUC: %.1e | APEX: %.1f | FWHM: %.1e", AUC, APEX, FWHM$fwhm),
               ylim = range(c(x$intensities, x$yp)),
               xlim = range(c(x$xx, x$times)),
               #xlim = range(APEX - 2 * FWHM$fwhm, APEX + 2 * FWHM$fwhm),
               main = paste(names(iRTmz())[which(x$mass == iRTmz())], x$mass));
          lines(x$xx, x$yp, col='red');
          segments(FWHM$x1, FWHM$y1, FWHM$x1 + FWHM$fwhm, FWHM$y1, col = 'green')
          abline(v = APEX, col = 'blue')
          legend("topleft", legend = c(sprintf("R^2: %.1e", r.squared)))
        }else{
          plot(x$times, x$intensities,main = paste(names(iRTmz())[which(x$mass == iRTmz())], x$mass))
          lines(x$xx, x$yp, col='violet');
        }
        x
      })
  }, height = 600)
  
  output$plotDDAiRTfits <- renderPlot({
    shiny::req(input$showIrtMS1Profile)
    shiny::req(rtFittedAPEX())
    
    iRTscore <- c(-24.92, 19.79, 70.52, 87.23, 0, 28.71, 12.39, 33.38, 42.26, 54.62, 100)
    
    rtFittedAPEX <- rtFittedAPEX() |>
      sapply(function(x){x$xx[which.max(x$yp)[1]]})
    
    fit <- lm(rtFittedAPEX ~ iRTscore)
    par(mfrow = c(1, 1), mar=c(5, 5, 4, 1))
    plot(rtFittedAPEX ~ iRTscore,
         ylab = 'Retention time [min]',
         xlab = "iRT score",
         pch=16, frame.plot = FALSE)
    abline(fit, col = 'grey')
    abline(v = 0, col = "grey", lty = 2)
    legend("topleft", legend = paste("Regression line: ", "rt =",
                                     format(coef(fit)[1], digits = 4), " + ",
                                     format(coef(fit)[2], digits = 2), "score",
                                     "\nR2: ", format(summary(fit)$r.squared, digits = 4)),
           bty = "n", cex = 0.75)
    text(iRTscore, rtFittedAPEX, names(iRTmz()), pos = 4,cex = 0.5)
  })
  
  
  output$iRTpeptides <- renderUI({
    if(input$showIrtMS2Profile){
      selectInput('iRTpeptides', 'iRTpeptides',
                  files(),
                  multiple = FALSE,
                  selected = iRTpeptides()[1])
    }
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
  
 
  
  #### cometVariable ------------
 
  
  
  output$cometVariable <- renderUI({
    
    defaulVariables <- c('nConfidentProteins', 'nConfidentPeptides', 'nMS2')
    
    selectInput('cometVariables', 'Variables',
                cometVariables(),
                multiple = TRUE,
                selected = defaulVariables)
    
  })
  
  
  #### comet lattice::xyplot -----------------
  output$cometPlot <- renderPlot({
    shiny::req(cometData())
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Plotting comet data ...")
    on.exit(progress$close())
    
    .lattice(cometData(),
             useBfabric = input$useBfabric,
             bfabricInstrumentEvents = BFabric$bfabricInstrumentEventsFiltered()$time,
             group = scanType)
  
  })# height = function(){400 * length(cometData()$instrument |> unique())})
  
  #### DIA-NN lattice::xyplot -----------------
  output$diannPlot <- renderPlot({
    shiny::req(diannData())
    
    .lattice(diannData(), input$useBfabric, BFabric$bfabricInstrumentEventsFiltered()$time)
    
  })#, height = function(){400 * length(diannData()$Instrument |> unique())})
 
  
  ## Summary ============
  summaryData <- reactive({
    #d1 <- data.frame(time = autoQC01wide()$time,
    #                 size = autoQC01wide()$size,
    #                 Instrument = autoQC01wide()$Instrument,
    #                 method = "Biognosys iRT (autoQC01)")
    
    d2 <- data.frame(time = cometWide()$time,
                     size = cometWide()$size,
                     Instrument = cometWide()$Instrument,
                     method = "DDA (comet)")
    
    d3 <- data.frame(time = diannWide()$time,
                     size = diannWide()$Size,
                     Instrument = diannWide()$Instrument,
                     method = "DIA (DIA-NN)")
    
    
    (list(d2, d3) |>
        Reduce(f = rbind))
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
    shiny::req(summaryData())
   
    trellis.par.set("superpose.symbol", superpose.symbol())
    
    lattice::dotplot(Instrument ~ time | method,
                     groups = method,
                     data = summaryData(),
                     alpha = 0.2,
                     cex = 2.4,
                     pch = 22,
                     layout = c(1, 3),
                     main = "instrument events")
  })
  
  output$plotSummaryCumsum  <- renderPlot({
    shiny::req(summaryData())
    
    trellis.par.set("superpose.symbol", superpose.symbol())
    
    
    barchart(size/1024^3 ~ format(time, "%Y-%m") | method, data = summaryData(),
             layout = c(1,3),
             scales = list(x = list(rot=45)))
  })
  

  
  ## printSummary --------
  output$summary <- renderPrint({
    capture.output( table(summaryData()$method, summaryData()$Instrument))
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
  
  output$autoQC01BfabricInstrumentEventsOutput <- renderUI({
    shiny::req(input$useBfabric)
    shiny::req(BFabric$bfabricInstrumentEventsFiltered())
    
    DT::renderDataTable({ BFabric$bfabricInstrumentEventsFiltered() })
    
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
