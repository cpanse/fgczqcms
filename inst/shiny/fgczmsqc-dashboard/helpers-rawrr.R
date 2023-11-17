#R

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
