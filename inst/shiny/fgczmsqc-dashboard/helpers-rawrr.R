#R

.rawrr_logo <- function(){
  ## compose colormap
  'graphics/rawrr_logo-colormap.rgb' |>
    read.csv(header = TRUE, sep = ' ') -> cm
  
  rgb(cm$red, cm$green, cm$blue, alpha = 1) -> cm.rgb
  
  ## cast to matrix
  scan('graphics/rawrr_logo144x153.txt') |>
    matrix(144, 153) -> m
  
  ## flip x axis
  m %*% diag(153)[,153:1] |> 
    image(asp = 1,axes = FALSE, useRaster = TRUE, col = cm.rgb)
  
  ## print colormap 
  #  points(x <- seq(0, 1, length = length(cm.rgb)),
  #         y <- rep(0, length(cm.rgb)),
  #         col = cm.rgb,
  #         pch=22)
}

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

#' When is a peak a peak?
#'
#' @description
#' A super basic heuristic to descide if x describes a Gaussian peak.
#' 
#' @param FWHM proposal
#' @param x curve
#'
#' @return TRUE if at least one point is in between the FWHM x-range.
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

.plotGaussianPeakProfile <- function(x){
  AUC <- sum(diff(x$xx) * (head(x$yp, -1) + tail(x$yp,  -1))) / 2
  APEX <- x$xx[which.max(x$yp)[1]]
  # peptide <- names(vals$mZ)[which(x$mass == vals$mZ)]
  # progress$set(detail = paste0("Render ", peptide), value = 4/5)
  
  FWHM <- .fwhm(x$xx, x$yp)
  r.squared <- x$r.squared[1]
  
  if (.isPeak(FWHM, x$times)){
    plot(x$times, x$intensities,
         type='p',
         sub = sprintf("AUC: %.1e | APEX: %.1f | FWHM: %.1e", AUC, APEX, FWHM$fwhm),
         ylim = range(c(x$intensities, x$yp)),
         xlim = range(APEX - 2 * FWHM$fwhm, APEX + 2 * FWHM$fwhm),
         main = paste0(x$mass, collapse = " | "));
    # legend("topleft", legend = c(sprintf("R^2: %.1e", r.squared)), cex = 0.5)
    lines(x$xx, x$yp, col='red');
    segments(FWHM$x1, FWHM$y1, FWHM$x1 + FWHM$fwhm, FWHM$y1, col = 'green')
    abline(v = APEX, col = 'blue')
  }else{
    plot(x$times, x$intensities, main = paste0(x$mass), sub = 'fitting failed!')
  }
}
