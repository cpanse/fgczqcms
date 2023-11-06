#R
## server

# requirements ===========
stopifnot(require(readr),
          require(reshape2),
          require(shinydashboard),
          require(lattice),
          require(rawrr))

#stopifnot(require(bfabricShiny))


.plotChromatogramSet <- function (x, diagnostic = FALSE, ...) 
{
  stopifnot(attr(x, "class") == "rawrrChromatogramSet")
  if (attr(x, "type") == "xic") {
    plot(0, 0, type = "n", frame.plot = FALSE, xlab = "Retention Time [min]", 
         ylab = "Intensities", ylim = range(unlist(lapply(x, function(o) {
           o$intensities
         }))), ...)
    cm <- hcl.colors(length(x), "Set 2")
    mapply(function(o, co) {
      lines(o$times, o$intensities, col = co)
    }, x, cm)
    legend("topleft", as.character(sapply(x, function(o) {
      o$mass
    })), col = cm, pch = 16, title = "target mass [m/z]", 
    bty = "n", cex = 0.75)
    if (diagnostic) {
      legend("topright", legend = paste(c("File: ", "Filter: ", 
                                          "Type: ", "Tolerance: "), c(basename(attr(x, 
                                                                                    "file")), attr(x, "filter"), attr(x, "type"), 
                                                                      attr(x, "tol"))), bty = "n", cex = 0.75, text.col = "black")
    }
  }
  invisible(x)
}


# utils ===========

.qcpanel <- function(x, y, ...){
  lattice::panel.abline(h = quantile(y, c(0.05, 0.25, 0.5, 0.75, 0.9)),
               col='grey', lwd=c(1,2,3,2,1))
  lattice::panel.xyplot(x, y, ...)
}

.iqrPanel <- function(x, y, ...){
  q <- quantile(y, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  iqr <- q[3] - q[1]
  lattice::panel.abline(h = c(q[3] + 1.5 * iqr, q[1] - 1.5 * iqr),
                        col = 'lightgrey', lwd = 0.75)
  lattice::panel.abline(h = q,
                        col = c('grey', 'green', 'grey'), lwd = c(0.75, 2, 0.75))
  
  filter <- ((q[1] - 1.5 * iqr) < y & y < (q[3] + 1.5 * iqr))
  idx <- order(x[filter])
  
  lattice::panel.xyplot(x[filter][idx], y[filter][idx], ..., type = 'b', pch=16)
  lattice::panel.xyplot(x[!filter], y[!filter], ..., type = 'p', col='lightgrey')
}


.assignInstrument <- function(x){
  for (i in c('QEXACTIVE_1', 'LUMOS_1', 'LUMOS_2', 'EXPLORIS_1', 'EXPLORIS_2', 'FUSION_2')){
    idx <- grepl(i, x$File.Name) 
    x$Instrument[which(idx)] <- i
  }
  x
}

.tic <- function(f){
  message(f)
  f |>
    rawrr::readChromatogram(type='tic') |>
    plot()
}

.missing <- function(){
  plot(0,0, xlab = '', ylab = '', type = 'n', axes = FALSE)
  text(0, 0, "missing\ndata or plot", cex=5)
}

# define server logic ============
function(input, output, session) {
  #  reactives =============
  
  iRTmz <- reactive({
    iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384, 683.8282,
               683.8541, 699.3388, 726.8361, 776.9301)
    
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
      rawrr:::fitPeak.rawrrChromatogram()
  })
  
  
  ###### comet --------
  comet <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading comet data ...")
    on.exit(progress$close())
    e <- new.env()
    fn <- rootdir() |> file.path("comet.RData")  |> load(envir=e)
    
    e$comet$assignmentRate <- round (100 * e$comet$nConfidentPSM / e$comet$nPSM)
    e$comet$Time <- as.POSIXct(e$comet$time )
    
    cc <- c('md5', 'filename.y', 'Time', 'instrument',  'scanType', 'psmFdrCutoff',
            'nDecoyPSM', 'nConfidentPSM', 'nDecoyPeptide', 'nConfidentPeptide',
            'nDecoyProteins', 'nConfidentProteins', 'fdrPSM', 'fdrPeptide',
            'fdrProtein',  'size',   'nMS2', 'TIC')
    
    
    e$comet[grepl(input$regex, e$comet$filename.y), cc]
  })
  
  cometFilter <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Filtering comet data acc. time and instrument ...")
    on.exit(progress$close())
    
    now <- Sys.time()
    
    filter <- input$cometDays[1] <= difftime(now, comet()$Time, units = "days")  &
      difftime(now, comet()$Time, units = "days") < input$cometDays[2] &
      comet()$instrument %in% input$instrument 
    
    filter
  })
  
  cometLong <- reactive({
    rv <- comet() |>
      reshape2::melt(id.vars = c("md5", "filename.y", "Time", "instrument",
                                 "scanType"))
  })
  
  cometVariables <- reactive({
    cometLong()$variable |> unique()
  })
  
  cometData <- reactive({
    shiny::req( input$instrument)
    #shiny::req( input$cometVariables )
    shiny::req(cometLong())
    
    now <- Sys.time()
    cometLong()[cometLong()$instrument %in% input$instrument &
                  cometLong()$variable %in% input$cometVariables &
                  input$cometDays[1] <=  difftime(now, cometLong()$Time, units = "days") &
                  difftime(now, cometLong()$Time, units = "days") < input$cometDays[2], ]
  }) 
  
  
  ########### diann -----------
  wide <- reactive({
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
    
    S$Instrument <- NA
    S |> .assignInstrument()
  })
  
  long <- reactive({
    wide() |> reshape2::melt(id.vars = c("Md5", "File.Name", "Time", "Instrument"))
  })
  
  
  ####### instruments -----------
  instruments <- reactive({
    c(wide()$Instrument |> unique(),
      comet()$instrument |> unique()) |> 
      unique()
  })
  
  ####### variables -----------
  variables <- reactive({
    long()$variable |> unique()
  })
   
  # TODO(cp): rename to diannData <-
  diannData <- reactive({
    # shiny::req(input$variables)
    shiny::req(input$instrument)
    
    now <- Sys.time()
    long()[long()$Instrument %in% input$instrument &
             long()$variable %in% input$variables &
             input$diannDays[1] <=  difftime(now, long()$Time, units = "days") &
             difftime(now, long()$Time, units = "days") < input$diannDays[2], ]
  }) 
  
  rootdir <- reactive({
    cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  rootdirraw <- reactive({
    cands <- c("/srv/www/htdocs/", "/scratch/DIAQC/qc/dump", "/Users/cp/Downloads/dump/")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
  })
  
  files <- reactive({
    shiny::req(diannData())
    shiny::req(cometData())

    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Composing file list ...")
    on.exit(progress$close())
    
    c(diannData()$File.Name[diannData()$Instrument %in% input$instrument],
      cometData()$filename.y[comet()$instrument %in% input$instrument]) |>
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
  
  #  renderUIs =============
  output$diannTimeSlider <- renderUI({
    mintime <- min(long()$Time)
    now <- Sys.time()
    
    maxtime <- (1 + difftime(now, mintime, units = 'days') |>
              round() |>
              as.integer())
    
    sliderInput("diannDays", "Observation range in days:", min = 0,
                max = maxtime,
                value = c(0, min(28, maxtime)), width = "100%")
  })
  
  
  output$instrument <- renderUI({
    selectInput('instrument', 'Instruments',
                instruments(),
                multiple = TRUE,
                selected = instruments()[1])
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
          selectInput('iRTpeptide', 'iRTpeptide',
                      names(iRTmz()),
                      multiple = FALSE,
                      selected = names(iRTmz())[1]),
          sliderInput("rtSlider", "rtSlider", min = 0, max = 120,
                      #max = 1 + ((iRTprofileRawDDA()[[1]][['times']]) |> max() |> round()),
                      value = c(26, 29), width = "100%"))
  })
  
  output$variable <- renderUI({
    defaulVariables <- c('Precursors.Identified', 'Proteins.Identified', 'FWHM.RT')
    selectInput('variables', 'Variables',
                variables(),
                multiple = TRUE,
                selected = defaulVariables)
  })
  
  output$ticFileInput <- renderUI({
    shiny::req(files())
    
    L <- list(selectInput('ticfile', 'TIC file candidates',
                          files(),
                          multiple = TRUE,
                          selected = files()[1]),
              hr())
    
    
    return(L)
  })
  
  output$fileInput <- renderUI({
    shiny::req(files())
    
    L <- list(selectInput('file', 'Files',
                          files(),
                          multiple = FALSE,
                          selected = files()[1]),
              hr(),
              checkboxInput("showRawFileHeader", "show Raw File Header", value = FALSE, width = NULL),
              checkboxInput("showIrtMS1Profile", "show iRT Ms Profile", value = FALSE, width = NULL),
              checkboxInput("showIrtMS2Profile", "show iRT Ms2 Profile" , value = FALSE, width = NULL),
              hr())
  
  
   return(L)
  })
  
  ### render fileOutput -----------
  output$fileOutput <- renderUI({
    L <- list()
    
    if (input$showRawFileHeader){
      L <- append(L, fluidRow(
        box(verbatimTextOutput("rawFileHeader"), width = 800)
      ))
    }
    
    if (input$showIrtMS1Profile){
      L <- append(L, list(
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
      L <- append(L, list(
        fluidRow(shiny::htmlOutput("scanTypeUI")),
        fluidRow(box(plotOutput("plotDIAiRTprofiles"), width = "100%"))
      ))
    }
    
    return(L)
  })
  #### renderPlots DIA-NN -------------
  output$tableDIANN <- renderDataTable({ diannData() })
  
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
  
  output$plotDDAiRTprofiles <- renderPlot({
    shiny::req(input$showIrtMS1Profile)
    shiny::req(rtFittedAPEX())
    
    par(mfrow = c(2, 6), mar = c(4, 4, 4, 1))
    rtFittedAPEX <- rtFittedAPEX() |>
      lapply(function(x){
        plot(x$times, x$intensities,
             type='p',
             ylim = range(c(x$intensities,x$yp)),
             main = paste(names(iRTmz())[which(x$mass == iRTmz())], x$mass));
        lines(x$xx, x$yp, col='red'); x})
    
  })
  
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
    
    print(input$ticfile)
    
    rv <- file.path(rootdirraw(), input$ticfile) |>
      lapply(FUN = function(f){
        progress$set(message = "Reading", detail = f)
        rawrr::readChromatogram(f, type='tic') })
    
    rv
  })
  
  output$plotTIC <- renderPlot({
  
    shiny::req(ticData())
    
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Plotting TIC data ...")
    on.exit(progress$close())
    
    par(mfrow = c(length(input$ticfile), 1))
    
    ticData() |> lapply(function(x){plot(x)})
    }, height = function(){300 * length(input$ticfile)})
  
  #### DIA-NN lattice::xyplot -----------------
  output$diannPlot <- renderPlot({
    shiny::req(diannData())

    lattice::xyplot(value ~ Time | variable * Instrument,
                    group = Instrument,
                    data = diannData(),
                    scales = list(y = list(relation = "free")),
                    panel = .iqrPanel,
                    sub = "Interquantile range (IQR): inbetween grey lines; median green; outliers: lightgrey.",
                    auto.key = list(space = "bottom"))
    
  },
  height = function(){400 * length(diannData()$Instrument |> unique())})
  
  #### cometVariable ------------
  output$cometVariable <- renderUI({
    
    defaulVariables <- c('nConfidentProteins', 'nConfidentPeptide', 'nMS2')
    
    selectInput('cometVariables', 'Variables',
                cometVariables(),
                multiple = TRUE,
                selected = defaulVariables)
    
  })
  
  #### comet time slider -----------------
  output$cometTimeSlider <- renderUI({
    mintime <- min(cometLong()$Time)
    now <- Sys.time()
    
    maxtime <- (1 + difftime(now, mintime, units = 'days') |>
                  round() |>
                  as.integer())
    
    sliderInput("cometDays", "Observation range in days:", min = 0,
                max = maxtime,
                value = c(0, min(maxtime, 28)), width = "100%")
  })
  #### comet lattice::xyplot -----------------
  
  imgHeight <- reactive({input$height |> as.integer()})
  
  output$cometPlot <- renderPlot({
    shiny::req(cometData())

    if(cometData() |> nrow() > 0){
      lattice::xyplot(value ~ Time | variable * instrument,
                      group = scanType,
                      data = cometData(),
                      scales = list(y = list(relation = "free")),
                      panel = .iqrPanel,
                      sub = "Interquantile range (IQR): inbetween grey lines; median green; outliers: lightgray.",
                      auto.key = list(space = "bottom"))
    }
    else{.missing()}
  },
  height = function(){400 * length(cometData()$instrument |> unique())})

  rawFileHeader <- reactive({
    shiny::req(input$file)
    
    rootdirraw() |>
      file.path(input$file) |>
      rawrr::readFileHeader()
    })
  
  #### print rawFileHeader  ----
  output$rawFileHeader <- renderPrint({
    shiny::req(rawFileHeader())
    
    capture.output(rawFileHeader())
  })
  
  
  #### sessionInfo ----
  output$sessionInfo <- renderPrint({
    capture.output(sessionInfo())
  })
  
  
}
