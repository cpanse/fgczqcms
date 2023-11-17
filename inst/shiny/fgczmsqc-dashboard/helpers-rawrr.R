#R

.rawrrChromatogramSet <- function (x, diagnostic = TRUE,
                                   xlim = range(unlist(lapply(x,
                                       function(o) { o$times }))), ...) 
{
  stopifnot(attr(x, "class") == "rawrrChromatogramSet")
  if (attr(x, "type") == "xic") {
    plot(0, 0, type = "n", xlim = xlim,
         ylim = range(unlist(lapply(x, function(o) {
                                                        o$intensities
                                                      }))), frame.plot = FALSE, xlab = "Retention Time [min]", 
         ylab = "Intensities", ...)
    cm <- hcl.colors(length(x), "Set 2")
    mapply(function(o, co) {
      lines(o$times, o$intensities, col = co)
      points(o$times, o$intensities, col = co, cex=0.5)
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

.isFWHM<-function(FWHM, x) {
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

.extractMaximumPeak <- function (y) 
{
  idx <- which(y == max(y))[1]
  
  seq(idx - 10 , idx + 10)
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
