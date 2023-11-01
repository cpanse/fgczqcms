#R
## server

# requirements ===========
stopifnot(require(readr), require(reshape2), require(shinydashboard))

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
    rawrr::readChromatogram( type='tic') |>
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
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading iRT peptide profiles ...")
    on.exit(progress$close())
    
    file.path(rootdirraw(), input$file) |>
      rawrr::readChromatogram(mass = iRTmz(),
                              tol = as.integer(input$ppmError),
                              type = "xic",
                              filter = "ms") 
  })

  #### iRTprofileRawDIA --------------
  iRTprofileRawDIA <- reactive({
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Reading Ms2 profiles ...")
    on.exit(progress$close())
    
    yIonSeries <- ("ELVIS" |> 
        protViz::fragmentIon())[[1]]['y'] |>
        unlist() 
    
    
    yIonSeries <- (input$iRTpeptide |> 
      protViz::fragmentIon())[[1]]['y'] |>
        unlist() 
    
    print(yIonSeries[seq(1, nchar(input$iRTpeptide) - 1)])
    
    file.path(rootdirraw(), input$file) |>
     rawrr::readChromatogram(mass = yIonSeries[seq(1, nchar(input$iRTpeptide) - 1)],
                              tol = as.integer(input$ppmError),
                              type = "xic",
                              filter = input$scanType) 
  })
  
  
  scanType <- reactive({
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
    fn <- rootdir() |> file.path("comet-20231030.RData")  |> load(envir=e)
    
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
  
  data <- reactive({
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
    progress <- shiny::Progress$new(session = session)
    progress$set(message = "Composing file list ...")
    on.exit(progress$close())
    
    c(data()$File.Name[data()$Instrument %in% input$instrument],
      comet()$filename.y[comet()$instrument %in% input$instrument]) |>
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
                value = c(0, min(28, maxtime)), width = 1000)
  })
  
  
  output$instrument <- renderUI({
    selectInput('instrument', 'Instruments',
                instruments(),
                multiple = TRUE,
                selected = instruments()[1])
  })
  
  output$scanType <- renderUI({
    if(input$plotDiannMs2){
      list( selectInput('scanType', 'scanType',
                        scanType(),
                        multiple = FALSE,
                        selected = scanType()[1]),
            selectInput('iRTpeptide', 'iRTpeptide',
                        names(iRTmz()),
                        multiple = FALSE,
                        selected = names(iRTmz())[1]),
            sliderInput("rtSlider", "rtSlider", min = 0,
                        max = 1 + ((iRTprofileRawDDA()[[1]][['times']]) |> max() |> round()),
                        value = c(26,29), width = 1000))
    }else{
      shiny::renderText("set off")
    }
  })
  
  output$variable <- renderUI({
    defaulVariables <- c('Precursors.Identified', 'Proteins.Identified')
    selectInput('variables', 'Variables',
                variables(),
                multiple = TRUE,
                selected = defaulVariables)
  })
  
  output$file <- renderUI({
    selectInput('file', 'Files',
                files(),
                multiple = FALSE,
                selected = files()[1])
  })
  
  #### renderPlots DIA-NN -------------
  output$tableDIANN <- renderDataTable({ data() })
  
  output$tableComet <-  DT::renderDataTable({ cometData()  })
  
  output$plotiRTDDAChromatograms <- renderPlot({
    iRTprofileRawDDA() |>
      plot(main = input$file)
  })
  
  
  #### plotDIAiRTprofiles  ------------
  output$plotDIAiRTprofiles <- renderPlot({
    if(input$plotDiannMs2){
      
      message("calling iRTprofileRawDIA()  ...")
      
      
      iRTprofileRawDIA()  |>
        .plotChromatogramSet(, xlim=input$rtSlider)
      
    }else{
      shiny::renderText("set off")
    }
  })
  
  output$plotDDAiRTprofiles <- renderPlot({
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
    if(input$plotDiannMs2){
      selectInput('iRTpeptides', 'iRTpeptides',
                  files(),
                  multiple = FALSE,
                  selected = iRTpeptides()[1])
    }
  })
  

  output$plotTIC <- renderPlot({ .tic(file.path(rootdirraw(), input$file)) })
  
  #### DIA-NN lattice::xyplot -----------------
  output$diannPlot <- renderPlot({
    if(data() |> nrow() > 0){
      lattice::xyplot(value ~ Time | variable * Instrument,
                      group = Instrument,
                      data = data(),
                      scales = 'free',
                      type = 'b',
                      #layout = c(1, length(input$variables)),
                      auto.key = list(space = "bottom"))}
    else{.missing()}
  })
  
  #### cometVariable ------------
  output$cometVariable <- renderUI({
    
    defaulVariables <- c('nMS2')
    
    selectInput('cometVariables', 'Variables',
                cometVariables(),
                multiple = FALSE,
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
                value = c(0, min(maxtime, 28)), width = 1000)
  })
  #### comet lattice::xyplot -----------------
  output$cometPlot <- renderPlot({
    if(cometData() |> nrow() > 0){
      lattice::xyplot(value ~ Time | variable * instrument,
                      group = scanType,
                      data = cometData(),
                      scales = 'free',
                      type = 'b',
                      pch = 16,
                      #layout = c(1, length(input$variables)),
                      auto.key = list(space = "bottom"))
    }
    else{.missing()}
  })
  
  
}
