#R


stopifnot(require(rawrr))

source("helpers-rawrr.R")

rawrrUI <- function(id){
  tagList(
    shinydashboard::box(
      column(12, offset = 0,
             fluidRow(
               column(6, offset = 0,
                      fluidRow(plotOutput(NS(id, "plotChromatograms"))),
                      fluidRow(sliderInput(NS(id, "rtwindow"), "rt window",
                                           min = 0.0, max = 40, value = c(10, 29), step = 0.5,
                                           width = "100%"))),
               column(6, offset = 0,
                      fluidRow(plotOutput(NS(id, "plotTools"))),
                      column(5, offset = 1,
                             fluidRow(radioButtons(NS(id, "ppmError"), "ppm error",
                                                   choices = c(5, 10, 20, 30, 50, 100),
                                                   selected = 10,
                                                   inline = TRUE,
                                                   width = 300),
                             )
                      )
               )
             ),
             fluidRow(
               plotOutput(NS(id, "plotProfiles"))
             )
      ), title = "rawrr module",
      footer = "The graphs trace the AUC | APEX | FWHM value computation 
      by extracting and
      fitting chromatographic profiles of each given m/z value. The rawrr
      package extracts the ion chromatograms from the selected Orbitrap raw
      file.",
      collapsible = TRUE,
      collapsed = FALSE,
      status = "primary",
      solidHeader = TRUE,
      width = 12)
  )
  
}

#' Extracts and plots chromatographic profile of a given m/z vector using the 
#' rawrr Rpkg from Bioconductor.
#'
#' @param id shiny namespace id
#' @param vals a reactive variable containing fn, the rawfile, and a numeric
#' vector of mZ values. The names of the vector, e.g., AA sequeces, are used as
#' labels.#'
#' @return nothing
rawrrServer <- function(id, vals){
  
  moduleServer(id,
               function(input, output, session) {
                 
                 observeEvent(vals$fn, {
                   message(paste0("fn has changed to ", vals$fn))
                 })
                 
                 rawfile <- reactive({ vals$fn })

                 peptideProfile <- reactive({
                   # invalidateLater(2000, session)
                   reactiveTimer(intervalMs = 2000, session)
                   #shiny::req(rawfile())
                   shiny::req(rawfile(), input$ppmError, vals$mZ)
                   #shiny::req(vals$mZ)
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Reading mZ profiles ...", detail = 'using rawrr')
                   on.exit(progress$close())
                   
                   if (file.exists(rawfile())){
                     message("Reading raw file: ", vals$fn)
                     rv <- rawrr::readChromatogram(rawfile = rawfile(),
                                                   mass = vals$mZ,
                                                   tol = as.integer(input$ppmError),
                                                   type = "xic",
                                                   filter = "ms")
                     
                     message(paste0("length = ", length(rv)))
                     return(rv)
                   }
                 })
                 
                 output$plotChromatograms <- renderPlot({
                   shiny::req(peptideProfile())
                   
                   if (length(peptideProfile()) > 0){
                     progress <- shiny::Progress$new(session = session)
                     progress$set(message = "Plotting chromatograms ...")
                     on.exit(progress$close())
                     
                     message("Plotting chromatograms ...")
                     
                     peptideProfile() |> 
                       .rawrrChromatogramSet(main = basename(vals$fn),
                                             xlim = as.numeric(input$rtwindow))
                   }else{
                     # TODO(cp): why it never goes here???
                     # TODO(cp): put the rawrr logo here.
                     message("Plotting .missing() ...")
                     .missing()
                   }
                 })
                 
                 output$plotTools <- renderPlot({
                   op <- par(mfrow = c(2, 2), mar=c(4,4,4,1))
                   .rawrr_logo()
                   
                   peptideProfile() |>
                     lapply(function(x)diff(x$times)) |>
                     unlist() |> hist(xlab = "diff(times)", main = "Histogram")
                   
                   
                 })
                 
                 output$plotProfiles <- renderPlot({
                   if (length(peptideProfile()) > 0){
                     progress <- shiny::Progress$new(session = session)
                     progress$set(message = "Plotting peptide profiles ...",
                                  detail = 'pick- and fit-ing peaks ...', value = 2/5)
                     on.exit(progress$close())
                     
                     
                     if(length(vals$mZ) < 7){
                       op <- par(mfrow = c(1, length(vals$mZ)))
                     }else{
                       op <- par(mfrow = c(2, 6))
                     }
                     

                     peptideProfile() |>
                       .pickPeak.rawrrChromatogram() |>
                       .fitPeak.rawrrChromatogram(delta = 0.5, n = 200) |>
                       lapply(FUN = .plotGaussianPeakProfile)
                   }else{
                     .missing()
                   }
                 })
               }
  )
}

