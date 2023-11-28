#R


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