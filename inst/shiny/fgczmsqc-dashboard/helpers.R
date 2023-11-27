#R
## helper functions for fgcz ms qc shiny dashboard


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

.fitChromatographicPeak <- function (x, y, delta = 0.25, n = 200) 
{
  peak <- data.frame(logy = log(y + 1), x = x)
  x.mean <- mean(peak$x)
  peak$xc <- peak$x - x.mean
  weights <- y^2
  fit <- lm(logy ~ xc + I(xc^2), data = peak, weights = weights)
  x0 <- -fit$coefficients[2]/(2 * fit$coefficients[3])
  
  ## predict
  xx <- with(peak, seq(min(xc) - delta, max(xc) + delta, length = n))
  yp <- exp(predict(fit, data.frame(xc = xx)))
  data.frame(xx = xx + x.mean, yp = yp, r.squared = summary(fit)$r.squared)
}

.fitPeak.rawrrChromatogram <- function (x, delta = 0.25, n = 200) 
{
  lapply(x, function(y) {
    fittedPeak <- .fitChromatographicPeak(y$times, y$intensities, delta = delta, n)
    rv <- y
    rv$xx <- fittedPeak$xx
    rv$yp <- fittedPeak$yp
    rv$r.squared <- fittedPeak$r.squared
    rv
  })
}

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
  # lattice::panel.abline(...)
  lattice::panel.xyplot(x[filter][idx], y[filter][idx], ..., type = 'b', pch=16)
  lattice::panel.xyplot(x[!filter], y[!filter], ..., type = 'p', col='lightgrey')
}

#' 
.lattice <- function(S, useBfabric = FALSE, bfabricInstrumentEvents, ...) {
  if (nrow(S) > 0){
    t <- trellis.par.get("superpose.symbol")
    cm <- c("#94C6FF", "#FFBBA9", "#76E3B8", "#FFD6AD", 
            "#BCE1FF", "#FFF691", "#FFC1E1")
    t$col <- cm
    t$fill <- cm
    trellis.par.set("superpose.symbol", t)
    
    lattice::xyplot(value ~ time | variable * Instrument,
                    data = S,
                    scales = list(y = list(relation = "free")),
                    panel = function(x, y, ...){
                      .iqrPanel(x, y, ...)
                      try(if (useBfabric){
                        lattice::panel.abline(v = bfabricInstrumentEvents, col = '#FF1111')
                      }, TRUE)
                    },
                    sub = "Interquantile range (IQR): inbetween grey lines; median green; outliers: lightgrey.",
                    auto.key = list(space = "bottom"), ...)
  }else{
    .missing()
  }
}

#' @return a list of B-Fabric instrumentids
.getInstruments <- function(){
  list(EXPLORIS_2 = 335,
       EXPLORIS_1 = 253,
       LUMOS_2 = 252,
       TIMSTOF_1 = 218,
       LUMOS_1 = 214,
       QEXACTIVE_1 = 93,
       FUSION_2 = 73,
       QEXACTIVEHFX_1 = -1,
       QEXACTIVEHF_1 = -1,
       QEXACTIVEHF_2 = -1,
       QEXACTIVEHF_4 = -1,
       QEXACTIVE_2 = -1,
       FUSION_1 = -1,
       VELOS_1 = -1)
}


#' Assign instruments of a given FGCZ file path
#'
#' @param x a data frame with a column 'File.Name'
#' @param coln the column name of the file name
#'
#' @return a data.frame with assigned instruments
.assignInstrument <- function(x, coln = 'File.Name'){
  stopifnot('Instrument' %in% colnames(x))
  
  for (i in names(.getInstruments())){
    idx <- grepl(i, x[[coln]]) 
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

.testBfabric <- function(){
  bfabricShiny::readPages(login,
                          webservicepassword,
                          endpoint = 'user',
                          query = list(login = 'cpanse'),
                          posturl = bfabricposturl) -> rv
  
  stopifnot(rv[[1]]$login == 'cpanse')
  rv[[1]]$login == 'cpanse'
}


.iRTmz <- function(){
  iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384,
             683.8282, 683.8541, 699.3388, 726.8361, 776.9301)
  
  names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                    "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                    "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                    "LFLQFGAQGSPFLK")
  
  iRTmz
}

.iRTscores <- function(){
  iRTscore <- c(-24.92, 19.79, 70.52, 87.23, 0, 28.71, 12.39, 33.38, 42.26, 54.62, 100)
  names(iRTscore) <- names(.iRTmz())
  iRTscore
}


.readDIANN <- function(filename){
  filename |>
    readr::read_delim(
      delim = ";",
      escape_double = FALSE,
      col_types = cols(Time = col_datetime(format = "%s"),
                       Size = col_integer(),
                       Precursors.Identified = col_integer(),
                       Proteins.Identified = col_integer()), 
      trim_ws = TRUE) -> S
  
  #S <- S[grepl(input$regex, S$File.Name), ]
  
  colnames(S)[colnames(S) == "Time"] <- "time"
  
  #idx <- order(S$time)
  S$Instrument <- NA
  S |>
    .assignInstrument() |>
    reshape2::melt(id.vars = c("Md5", "File.Name", "time", "Instrument"))
}



#' read comet qc data from RData file
#'
#' @param filename 
#'
#' @return a data.frame object
#'
#' @examples
#' .readComet('~/Downloads/dump/comet.RData') -> S
.readComet <- function(filename){
  e <- new.env()
  filename |>
    load(envir = e)
  
  e$comet$assignmentRate <- round (100 * e$comet$nConfidentPSM / e$comet$nPSM)
  
  ## rename columns
  colnames(e$comet)[colnames(e$comet) == "filename.y"] <- "File.Name"
  colnames(e$comet)[colnames(e$comet) == "Time"] <- "time"
  colnames(e$comet)[colnames(e$comet) == "md5"] <- "Md5"
  colnames(e$comet)[colnames(e$comet) == "instrument"] <- "Instrument"
  colnames(e$comet)[colnames(e$comet) == "nConfidentPeptide"] <- "nConfidentPeptides"
  
  cc <- c('Md5', 'File.Name', 'time', 'Instrument',  'scanType', 'psmFdrCutoff',
          'nDecoyPSM', 'nConfidentPSM', 'nDecoyPeptide', 'nConfidentPeptides',
          'nDecoyProteins', 'nConfidentProteins', 'fdrPSM', 'fdrPeptide',
          'fdrProtein',  'size',   'nMS2', 'TIC')
  e$comet$time <- as.POSIXct(e$comet$time )
  
  #idx <- order(e$comet$time)
  e$comet[, cc]  |>
    reshape2::melt(id.vars = c("Md5", "File.Name", "time", "Instrument", "scanType")) 
}

.status <- function(f){
  dt <- difftime(Sys.time(), file.mtime(f),  unit='hour') 
  if (dt > 12){
    "danger"
  }else if (dt < 1){
    "success"
  }else{
    "primary"
  }
}
