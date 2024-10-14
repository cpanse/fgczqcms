#Rst.files('/scratch/cpanse/comet/', pattern="comet.txt$")
## $ Christian Panse <cp@fgcz.ethz.ch>, 2019-2024$
## 2024-10-14

COMETOUTPUTDIR <- '/scratch/cpanse/comet/'

stopifnot(dir.exists(COMETOUTPUTDIR),
  require(protViz),
  require(parallel))

# INPUT
F <- list.files(COMETOUTPUTDIR, pattern=".*comet.txt$")

assignInstrument <- function(df){
    stopifnot(is.data.frame(df))

    df$instrument <- NA
    for (p in c("FUSION_1", "FUSION_2", "G2HD_1", "LC1100",
                "LTQ_1", "LTQFT_1", "ORBI_1", "ORBI_2", "PROTEONXPR36",
                "QEXACTIVE_1", "QEXACTIVE_2", "QEXACTIVE_3", "QEXACTIVEHF_4",
                "QEXACTIVEHF_2", "QEXACTIVEHFX_1", "QTRAP_1", "T100_1",
                "TOFTOF_2", "TRIPLETOF_1", "TIMSTOFFLEX_1", "TSQ_1", "TSQ_2", "VELOS_1", "LUMOS_1", "LUMOS_2",
                "VELOS_2", "EXPLORIS_1", "EXPLORIS_2"
    )){
        df$instrument[grep(p, df$filename)] <- p
    }
    # df$POSIXct <- as.POSIXct(df$time, origin="1970-01-01")
    df
}


autoQC4L <- read.table("pfiles_autoQC4L.txt", sep=';', col.names=c('md5', 'time', 'size', 'filename'))
autoQC4L$m <- gsub(".raw$", "", sapply(strsplit(as.character(autoQC4L$filename), '/'), function(x){paste(x[3],x[5], sep='.')}))
autoQC4L$POSIXct <- as.POSIXct(autoQC4L$time, origin="1970-01-01")
autoQC4L <- assignInstrument(autoQC4L)

comet <- parallel::mclapply(F, function(f){
    df <- read.table(file.path(COMETOUTPUTDIR, f), header = TRUE, fill = TRUE, skip=1)
    if(nrow(df) > 0){
      # message(paste0("Reading file '", f, "'..."))
      S <- protViz::summary.cometdecoy(df)
      S$filename <- f;
      return(S)
    }
    NULL
  },
  mc.cores=8) |>
  Reduce(f=rbind)

#comet <- do.call('rbind', comet)
comet$m <- sapply(strsplit(comet$filename, split='\\.', perl=TRUE), function(x){paste(x[1], x[2], sep='.')})
comet$scanType <- sapply(strsplit(comet$filename, split='\\.', perl=TRUE), function(x){paste(x[3], x[4], sep='.')})


comet <- merge(comet, autoQC4L, by.x='m', by.y='m')
comet$file <- basename(as.character(comet$filename.y))

## TODO: generalize
(system("cd /scratch/cpanse/comet && cat *.info   > info.txt"))

info <- read.table(file.path(COMETOUTPUTDIR, "info.txt"), sep=",", col.names=c('filename', 'scanType', 'nMS2', 'TIC', 'Ttmp'))
info$filename <- basename(as.character(info$filename))
info$m <- paste(gsub(".raw", "", info$filename), info$scanType, sep='.')
comet$m <- gsub(".comet.txt", "", comet$filename.x)
comet <- merge(comet, info[,c('m', 'nMS2','TIC')], by='m')   

#comet$time.x <- file.mtime(file.path(COMETOUTPUTDIR, comet$filename.x))

idx <- order(comet$instrument, comet$POSIXct)
comet <- comet[rev(idx), ]

source("/scratch/FGCZ-MS-QC-DASHBOARD_A331/qc/R/TIMSdda-comet-summary.R")

input <- c('/scratch/cpanse/autoQC4L/TIMSTOF_1/pfiles_TIMSTOF_1_autoQC4L.txt',
           '/scratch/cpanse/autoQC03dda/TIMSTOFFLEX_1/pfiles_TIMSTOF_1_autoQC03dda.txt')

harmonizeTimsDDA(input) |>
    rbind(comet) -> comet

#comet$nConfidentPSM <- as.integer(comet$nConfidentPSM)
save(comet, file='/scratch/cpanse/comet/comet.RData')

#rmarkdown::render('~/__checkouts/F1000_rawDiag/vignettes/autoQC4L_fgcz-r-033.Rmd', output_file='~/WWW/autoQC4L.html')
