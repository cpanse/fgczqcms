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
                      fluidRow(plotOutput(NS(id, "plotProfiles"))),
                      column(5, offset = 1,
                             fluidRow(radioButtons(NS(id, "ppmError"), "ppm error",
                                                   choices = c(5, 10, 20, 30, 50, 100),
                                                   selected = 10,
                                                   inline = TRUE,
                                                   width = 300),
                             )
                      )
               )
             )
      ), title = "rawrr module",
      footer = "extracts and plots chromatographic profiles of a given m/z list by reading the Orbitrap raw file.",
      collapsible = TRUE, width = 12)
  )
}

#' extracts and plots chromatographic profile of a given m/z list
#'
#' @param id shiny namespace id
#' @param rawfile the complete path to the raw file
#' @param mZ a vector of m/z values
#'
#' @return plots
rawrrServer <- function(id, vals){
  moduleServer(id,
               function(input, output, session) {
                 
                 observeEvent(vals$fn, {
                   message(paste0("fn has changed to ", vals$fn))
                 })
                 
                 rawfile <- reactive({vals$fn})
                 
                 peptideProfile <- reactive({
                   shiny::req(rawfile())
                   #shiny::req(rawfile(), input$ppmError, vals$mZ)
                   #shiny::req(vals$mZ)
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Reading mZ profiles ...")
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
                   #  return(rv)
                   #}else{
                   #  return(NULL)
                   #}
                 })
                 
                 output$plotChromatograms <- renderPlot({
                   if (length(peptideProfile()) >0){
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
                 
                 output$plotProfiles <- renderPlot({
                   if (length(peptideProfile()) > 0){
                     progress <- shiny::Progress$new(session = session)
                     progress$set(message = "Plotting peptide profiles ...",
                                  detail = 'pick- and fit-ing peaks ...', value = 2/5)
                     on.exit(progress$close())
                     
                     if (length(vals$mZ) > 9){
                       op <- par(mfrow = c(3, 4))
                     }else if (length(vals$mZ) > 4){
                       op <- par(mfrow = c(3, 3))
                     }else if (length(vals$mZ) > 2){
                       op <- par(mfrow = c(2, 2))
                     }else{
                       op <- par(mfrow = c(1, length(vals$mZ)))
                     }
                     
                     
                     ## TODO(cp): put that into the helper file
                     peptideProfile() |>
                       .pickPeak.rawrrChromatogram() |>
                       .fitPeak.rawrrChromatogram(delta = 0.5, n = 200) |>
                       lapply(function(x){
                         AUC <- sum(diff(x$xx) * (head(x$yp, -1) + tail(x$yp,  -1))) / 2
                         APEX <- x$xx[which.max(x$yp)[1]]
                         peptide <- names(vals$mZ)[which(x$mass == vals$mZ)]
                         progress$set(detail = paste0("Render ", peptide), value = 4/5)
                         
                         FWHM <- .fwhm(x$xx, x$yp)
                         r.squared <- x$r.squared[1]
                         #df <- data.frame(filename = file, time = time, peptide = peptide, auc = NA, apex = NA, FWHM = NA)
                         
                         if (.isFWHM(FWHM, x$times)){
                           plot(x$times, x$intensities,
                                type='p',
                                sub = sprintf("AUC: %.1e | APEX: %.1f | FWHM: %.1e", AUC, APEX, FWHM$fwhm),
                                ylim = range(c(x$intensities, x$yp)),
                                xlim = range(APEX - 2 * FWHM$fwhm, APEX + 2 * FWHM$fwhm),
                                main = paste0(peptide, x$mass, collapse = " | "));
                           # legend("topleft", legend = c(sprintf("R^2: %.1e", r.squared)), cex = 0.5)
                           lines(x$xx, x$yp, col='red');
                           segments(FWHM$x1, FWHM$y1, FWHM$x1 + FWHM$fwhm, FWHM$y1, col = 'green')
                           abline(v = APEX, col = 'blue')
                           
                           #msg <- paste0(file, peptide, AUC, APEX, FWHM$fwhm, sep='\t|\t')
                           #df <- data.frame(filename = file, time = time, peptide=peptide, auc=AUC, apex=APEX, FWHM=FWHM$fwhm)
                         }else{
                           plot(x$times, x$intensities, main = paste(peptide, x$mass), sub = 'fitting failed!')
                         }
                       }) # lapply
                   }else{
                     .missing()
                   }
                 })
               }
  )
}

