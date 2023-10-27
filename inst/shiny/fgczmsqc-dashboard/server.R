#R
## server

# requirements ===========
stopifnot(require(readr), require(reshape2), require(shinydashboard))


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
  rawrr::readChromatogram(f, type='tic') |>
    plot()
}

.missing <- function(){
  plot(0,0, xlab = '', ylab = '', type = 'n', axes = FALSE)
  text(0, 0, "missing\ndata or plot", cex=5)
}

# define server logic ============
function(input, output, session) {
  #  reactives =============
  iRTprofileRawDDA <- reactive({
    iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384, 683.8282,
               683.8541, 699.3388, 726.8361, 776.9301)
    
    names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                      "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                      "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                      "LFLQFGAQGSPFLK")
    
    rawrr::readChromatogram(input$file, mass = iRTmz, tol = as.integer(input$ppmError), type = "xic", filter = "ms") 
  })
  
  iRTprofileRawDIA <- reactive({
    #plot(0, 0, sub=f, type = 'n', xlab='', ylab=''); text(0,0, "t.b.d", cex=5)
    iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384, 683.8282,
               683.8541, 699.3388, 726.8361, 776.9301)
    
    names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                      "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                      "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                      "LFLQFGAQGSPFLK")
    
    yIonSeries <- names(iRTmz)[9] |>
      lapply(function(x){
        y <- protViz::fragmentIon(x)[[1]]$y
        y[seq(1, length(y) - 1)]
      }) |> unlist()
    
    input$file |>
      rawrr::readChromatogram(mass = yIonSeries, tol = 15, type = "xic", filter = "ms2") 
  })
  
 
  wide <- reactive({
    S <- rootdir() |> 
      file.path("output.txt") |> readr::read_delim(
                      delim = ";",
                      escape_double = FALSE,
                      col_types = cols(Time = col_datetime(format = "%s"),
                                       Size = col_integer(),
                                       Precursors.Identified = col_integer(),
                                       Proteins.Identified = col_integer()), 
                      trim_ws = TRUE) 
    
   
    
    S$Instrument <- NA
    S |> .assignInstrument()
  })
  
  long <- reactive({
    wide() |> reshape2::melt(id.vars = c("Md5", "File.Name", "Time", "Instrument"))
  })
  
  instruments <- reactive({
    wide()$Instrument |> unique()
  })
  
  variables <- reactive({
    long()$variable |> unique()
  })
  
  data <- reactive({
    now <- Sys.time()
    long()[long()$Instrument %in% input$instrument &
             long()$variable %in% input$variables &
             input$days[1] <=  difftime(now, long()$Time, units = "days") &
             difftime(now, long()$Time, units = "days") < input$days[2], ]
  }) 
  
  rootdir <- reactive({
    cands <- c("/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
    for (d in cands){
      if (dir.exists(d))return(d)
    }
    NULL
    })
  
  files <- reactive({
    data()$File.Name |>
      unique() |>
      lapply(function(f){file.path(rootdir(), f)}) |>
      lapply(function(f){
        if(file.exists(f)){return(f)}
        else{
          msg <- paste0("Can not find file ", f)
          message(msg)
        }
        NULL
        }) |>
      unlist()
  })
  
  #  renderUIs =============
  output$instrument <- renderUI({
    selectInput('instrument', 'Instruments',
                instruments(),
                multiple = TRUE,
                selected = instruments())
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
  
  # renderPlots ==========
  output$plot1 <- renderPlot({
    if(data() |> nrow() > 0){
      lattice::xyplot(value ~ Time | variable,
                      group = Instrument,
                      data = data(),
                      scales = 'free',
                      type = 'b',
                      layout = c(1, length(input$variables)),
                      auto.key = list(space = "bottom"))}
    else{.missing()}
    
  })
  
  output$table <- renderDataTable({ data() })
  
  output$plotiRTDDAChromatograms <- renderPlot({
    iRTprofileRawDDA() |> plot(main=gsub(rootdir(), "", input$file))
    
  })
  
  output$plotiRTprofiles <- renderPlot({
    par(mfrow = c(2, 6), mar=c(4,4,4,1))
    rtFittedAPEX <- iRTprofileRawDDA() |>
      rawrr:::pickPeak.rawrrChromatogram() |>
      rawrr:::fitPeak.rawrrChromatogram() |>
      lapply(function(x){
        plot(x$times, x$intensities, type='p',
             ylim=range(c(x$intensities,x$yp)),
             main=x$mass); lines(x$xx, x$yp,
                                 col='red'); x}) |>
      sapply(function(x){x$xx[which.max(x$yp)[1]]})
  })
  
  output$plotiRTfits <- renderPlot({
    iRTscore <- c(-24.92, 19.79, 70.52, 87.23, 0, 28.71, 12.39, 33.38, 42.26, 54.62, 100)
    
    rtFittedAPEX <- iRTprofileRawDDA() |>
      rawrr:::pickPeak.rawrrChromatogram() |>
      rawrr:::fitPeak.rawrrChromatogram() |>
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
    #text(iRTscore, rt, iRTmz, pos=1,cex=0.5)
  })
  
  output$plotiRTfits2 <- renderPlot({plot(0,0)})

}
