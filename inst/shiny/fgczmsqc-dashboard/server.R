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


# define server logic ============
function(input, output, session) {
  
  #  reactives =============
  wide <- reactive({
    S <- readr::read_delim("~/Downloads/output.txt", 
                      delim = ";", escape_double = FALSE, col_types = cols(Time = col_datetime(format = "%s"), 
                                                                           Size = col_integer(), Precursors.Identified = col_integer(), 
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
  
  # renderPlots ==========
  output$plot1 <- renderPlot({
    lattice::xyplot(value ~ Time | variable,
                    group = Instrument,
                    data = data(),
                    scales = 'free',
                    type = 'b',
                    layout = c(1, length(input$variables)),
                    auto.key = list(space = "bottom"))
  })
}
