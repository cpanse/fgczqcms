#R
## some ggplot2 templates for autoQC03 plots


.ggh4x <- function(gp, data){
  message("using 'ggh4x' package ...")
  if ('peptide' %in% colnames(data)){
    gp +
      ggh4x::facet_grid2(peptide ~ variable * Instrument,
                         scales = "free_y", independent = "y") +
      ggplot2::theme(legend.position = "none") -> gp
  }else if ('scanType' %in% colnames(data)){
    gp +
      ggh4x::facet_grid2(scanType ~ variable * Instrument,
                         scales = "free_y", independent = "y") +
      ggplot2::theme(legend.position = "none") -> gp
  }else{
    gp + 
      ggh4x::facet_grid2(. ~ variable * Instrument,
                         scales = "free_y", independent = "y") -> gp
  }
  gp
}

.defaultFacets <- function(gp, data){
  warning("install package 'ggh4x' to free y-axis scales.")
  if ('peptide' %in% colnames(data)){
    gp + ggplot2::facet_grid(peptide ~  variable * Instrument, scales = "free_y") +
      ggplot2::theme(legend.position = "none") -> gp
  }else if ('scanType' %in% colnames(data)){
    gp + ggplot2::facet_grid(scanType ~  variable * Instrument, scales = "free_y") +
      ggplot2::theme(legend.position = "none") -> gp
  }else{
    gp + ggplot2::facet_grid(. ~   variable * Instrument, scales = "free_y")  -> gp
  }
}

.ggplot <- function(data = NULL, variables = NULL){
  stopifnot(!is.null(data),
            !is.null(variables))
  data |> 
    subset(variable %in% variables) -> data 

  if ('peptide' %in% colnames(data)){
    formula <- value ~ peptide * Instrument * variable
  }else if ('scanType' %in% colnames(data)){
    formula <- value ~ scanType * Instrument * variable
  }else{
    formula <- value ~ Instrument * variable
  }

  hlineMedian <- aggregate(formula, FUN = median, data = data)

  irqL <- aggregate(formula,
                    data = data,
                    FUN = function(x){
                      q <- quantile(x, c(0.25, 0.75));
                      irq.low <- q[1] - 1.5 * diff(q)
                    })
  
  irqU <- aggregate(formula,
                    data = data,
                    FUN = function(x){
                      q <- quantile(x, c(0.25, 0.75));
                      irq.up <- q[2]  + 1.5 * diff(q)
                    })
  
  
  q25 <- aggregate(formula,
                   data = data,
                   FUN = function(x){
                     q <- quantile(x, c(0.25));
                     q
                   })
  
  q75 <- aggregate(formula,
                   data = data,
                   FUN = function(x){
                     q <- quantile(x, c(0.75));
                     q
                   })

  ## remove values less equal 0
  irqL |> 
    subset(irqL$value > 0) -> irqL

  data |> 
    ggplot2::ggplot(ggplot2::aes(time, value)) +
    ggplot2::geom_hline(data = hlineMedian, ggplot2::aes(yintercept = value), col = 'darkgreen') +
    ggplot2::geom_hline(data = irqL, ggplot2::aes(yintercept = value), colour = 'grey', lwd = 1) +
    ggplot2::geom_hline(data = irqU, ggplot2::aes(yintercept = value), colour = 'grey', lwd = 1) +
    ggplot2::geom_hline(data = q25, ggplot2::aes(yintercept = value), colour = 'cornflowerblue') +
    ggplot2::geom_hline(data = q75, ggplot2::aes(yintercept = value), colour = 'cornflowerblue') -> gp
#    ggplot2::labs(subtitle = "dargreen shows the median; the lower and upper grey lines indicate the inter-quantile range; the cornflowerblue represents 0.25 and 0.75 quantiles.") -> gp

  if (require("ggh4x")){
    .ggh4x(gp, data) -> gp
  }else{
    .defaultFacets(gp, data) -> gp
  }

  gp
}


.ggplotAutoQC03 <- function(data = NULL, variables = NULL, alpha = 0.7){
  .ggplot(data, variables) -> gp
  
  if ("scanType" %in% names(data)){
    gp +
      ggplot2::geom_point(ggplot2::aes(time, value, colour = scanType), alpha = alpha) +
      ggplot2::geom_line(ggplot2::aes(time, value, colour = scanType), alpha = 0.5 * alpha) -> gp
    
  }else{
    gp +
      ggplot2::geom_point(ggplot2::aes(time, value), alpha = 1.0) +
      ggplot2::geom_line(ggplot2::aes(time, value), alpha = 1.0) -> gp
  }
  gp + ggplot2::theme_light()
}

.ggplotAutoQC01 <- function(data = NULL, variables = NULL, alpha = 0.7){
  .ggplot(data, variables) +
    ggplot2::geom_point(ggplot2::aes(time, value, colour = peptide), alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(time, value, colour = peptide), alpha = alpha) -> gp
  gp + ggplot2::theme_light()
}
