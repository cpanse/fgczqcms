#R
stopifnot(require(rawrr), require(readr))

.extractMaximumPeak <- function (y) 
{
	  idx <- which(y == max(y))[1]
	    
	      seq(idx - 15 , idx + 15)
}

.pickPeak.rawrrChromatogram <- function (x) 
{
	  lapply(x, function(y) {
		      idx <- .extractMaximumPeak(y$intensities)
		          rv <- y
			      rv$times <- y$times[idx]
			          rv$intensities <- y$intensities[idx]
				      rv
				        })
}

.isPeak <- function(FWHM, x) {
	for (i.x in x){
        	if (FWHM$x1 < i.x & i.x < (FWHM$x1 + FWHM$fwhm))
		{
			return(TRUE)
		}
	}
	return(FALSE)
}

.fwhm <- function(x, y){
    ymax <- max(y)
    halfmax <- ymax / 2
    
    idxMax <- which(y == ymax)[1]
    
    for (i in 1:length(y)){
      if (y[i] >= halfmax){
        break
      }
    }
    
    fwhm <- 2 * (x[idxMax] - x[i])
    
    return(list(x1= x[i], y1 = y[i], 
                idxMax = idxMax,
                fwhm = fwhm))
  }

iRTmz <- function(){
  iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384,
             683.8282, 683.8541, 699.3388, 726.8361, 776.9301)
  
  names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                    "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                    "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                    "LFLQFGAQGSPFLK")
  
  iRTmz
}

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

.extractValues <- function(file, time, rootdir = "/srv/www/htdocs/", Ms1ppmError =10, outputfile){
  if (isFALSE(file.exists(file.path(rootdir, file)))){ return(NULL)}

  file.path(rootdir, file) |>
    rawrr::readChromatogram(mass = iRTmz(),
                            tol = as.integer(Ms1ppmError),
                            type = "xic",
                            filter = "ms") |>
    .pickPeak.rawrrChromatogram() |>
    .fitPeak.rawrrChromatogram(delta = 0.5, n = 400) |>
    lapply(function(x){
      AUC <- sum(diff(x$xx) * (head(x$yp, -1) + tail(x$yp,  -1))) / 2
      APEX <- x$xx[which.max(x$yp)[1]]
      peptide <- names(iRTmz())[which(x$mass == iRTmz())]
      FWHM <- .fwhm(x$xx, x$yp)
      r.squared <- x$r.squared[1]
      df <- data.frame(filename = file, time = time, peptide = peptide, auc = NA, apex = NA, FWHM = NA)
      
      if (.isPeak(FWHM, x$times)){
        df <- data.frame(filename = file, time = time, peptide=peptide, auc=AUC, apex=APEX, FWHM=FWHM$fwhm)
      }
      df 
    }) |>
    Reduce(f = rbind) -> rv


    write.table(rv, file = outputfile, sep = ";", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
}


input <- "/scratch/FGCZ-MS-QC-DASHBOARD_A331/qc/dump/autoQC01.input.txt"
output <- "/scratch/FGCZ-MS-QC-DASHBOARD_A331/qc/dump/autoQC01-fit-apex-auc-fwhm.txt"
#output <- "/tmp/autoQC01-fit-apex-auc-fwhm.txt"

Sin <- readr::read_delim(input,
        delim = ";",
        escape_double = FALSE,
        col_names =c('md5', 'time', 'size', 'filename'),
        col_types = cols(time = col_datetime(format = "%s"),
                         size = col_integer(),
                         n = col_integer()),
        trim_ws = TRUE)

Sout <- readr::read_delim(output, delim = ";",
   escape_double = FALSE, 
   col_names =c('filename', 'time', 'peptide', 'AUC', 'APEX', 'FWHM'),
   col_types = readr::cols(time = col_datetime(format = "%s")),
   trim_ws = TRUE) 



missing.idx <- which(!(Sin$filename %in% Sout$filename))

S <- Sin[missing.idx, ] |> tail(1000)

mapply(time = S$time, filename = S$filename,
	FUN = function(time, filename){
		dt <- difftime(Sys.time(), file.mtime(file.path("/srv/www/htdocs/", filename)),  unit = 'hour') |> as.numeric()

		if (dt < 72){
		   message(paste0("processing ", filename, " | mtime=", dt, "hours ..."))
		   try(.extractValues(file = filename, time = format(time, "%s"), outputfile = output))
		}else{
		   message(paste0("skipping ", filename, " because mtime is older than expected ", " | mtime=", dt, "hours ..."))
		}
	}) -> rv
quit('yes')
