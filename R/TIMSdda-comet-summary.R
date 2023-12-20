#R

## Christian Panse <cp@fgcz.ethz.ch> 2023-12-20
## to be used to harmonize the output of the comet summary
## for Orbitraps and TIMS devices

assignInstrument <- function(df){
  stopifnot(is.data.frame(df))
  
  df$instrument <- NA
  for (p in c("TIMSTOF_1", "TIMSTOFFLEX_1")){
    df$instrument[grep(p, df$filename)] <- p
  }
  df
}


#' Harmonize comet DDA searches from Bruker TIMS devices
#'
#' @return a data.frame with the following to be used for the qc dashboard
#'
#' @examples
harmonizeTimsDDA <- function(input){
  input |>
    lapply(FUN = function(f){
      stopifnot(file.exists(f))
      x <- scan(f, what = character())
      data.frame(
        md5 = rep(NA, length(x)),
        time = file.mtime(x),
        size = file.size(x),
        filename = as.character(gsub("/srv/www/htdocs/" ,"", x))
      )
    }) |>
    Reduce(f = rbind) |>
    assignInstrument() -> x
  
  
  
  input |>
    dirname() |>
    lapply(FUN = function(d){file.path(d, list.files(d, pattern = ".*comet.txt", recursive = TRUE))}) |>
    unlist() |>
    parallel::mclapply(FUN = function(f){
      x <- read.table(file.path(f), header = TRUE, fill = TRUE, skip=1)  
      protViz::summary.cometdecoy(x) -> S
      S$scanType <- "TIMS"
      S$nMS2 <- nrow(x)
      S$filename <- f
      S$TIC <- NA
      S$m <- f
      S
    }) |> 
    Reduce(f = rbind) -> y
  
  y$m <- strsplit(y$m, '/') |> vapply(FUN = function(x){x[length(x)-1]}, FUN.VALUE="20210302_015_autoQC4L_88min_S1-A2_1_1424")
  
  x$m <- gsub(".d.zip", "", x$filename) |> strsplit('/')  |> vapply(FUN = function(x){x[length(x)]}, FUN.VALUE="20210302_015_autoQC4L_88min_S1-A2_1_1424")
  
  merge(x, y, by = "m") -> S
  
  S$POSIXct <- S$time |> as.POSIXct()
  S$file <- basename(S$filename.x)
  
  cc <- c('m', 'nPSM', 'psmFdrCutoff', 'nDecoyPSM', 'nConfidentPSM', 'nDecoyPeptide', 'nConfidentPeptide', 'nDecoyProteins', 'nConfidentProteins', 'fdrPSM', 'fdrPeptide', 'fdrProtein', 'filename.x', 'scanType', 'md5', 'time', 'size', 'filename.y', 'POSIXct', 'instrument', 'file', 'nMS2', 'TIC')
  S[, cc]
}

input <- c('/scratch/cpanse/autoQC4L/TIMSTOF_1/pfiles_TIMSTOF_1_autoQC4L.txt',
           '/scratch/cpanse/autoQC03dda/TIMSTOFFLEX_1/pfiles_TIMSTOF_1_autoQC03dda.txt')

harmonizeTimsDDA(input)


