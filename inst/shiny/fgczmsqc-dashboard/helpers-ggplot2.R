#R
## some ggplot2 templates for autoQC03 plots
.ggplot <- function(data = NULL, variables = NULL){
  stopifnot(!is.null(data),
            !is.null(variables))
  
  #q <- quantile(y, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  
    
  data |> 
    subset(variable %in% variables) -> data 
  
  if ('peptide' %in% colnames(data)){
    hlineMedian <- aggregate(value ~ peptide * Instrument * variable, FUN = median, data = data)
  }else{
    hlineMedian <- aggregate(value ~ Instrument * variable, FUN = median, data = data)
  }
  
  irqL <- aggregate(value ~ Instrument * variable, data = data,
                    FUN = function(x){
                      q <- quantile(x, c(0.25, 0.75));
                      irq.low=q[1] - 1.5 * diff(q)
                    })
  irqU <- aggregate(value ~ Instrument * variable, data = data,
                    FUN = function(x){
                      q <- quantile(x, c(0.25, 0.75));
                      q[2]  + 1.5 * diff(q)
                    })
  
  print(irqL)
  data |> 
    ggplot2::ggplot(ggplot2::aes(time, value)) +
    ggplot2::geom_hline(data = hlineMedian, ggplot2::aes(yintercept = value), col = 'darkgreen') +
    ggplot2::geom_hline(data = irqL, ggplot2::aes(yintercept = value), colour = 'grey') +
    ggplot2::geom_hline(data = irqU, ggplot2::aes(yintercept = value), colour = 'grey')  -> gp
    
    
  if ('peptide' %in% colnames(data)){
    gp + ggplot2::facet_grid(peptide ~  variable * Instrument, scales = "free_y") +
      ggplot2::theme(legend.position = "none") -> gp
  }else{
    gp + ggplot2::facet_grid(. ~   variable * Instrument, scales = "free_y")  -> gp
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
  gp
}

.ggplotAutoQC01 <- function(data = NULL, variables = NULL, alpha = 0.7){
  .ggplot(data, variables) +
    ggplot2::geom_point(ggplot2::aes(time, value, colour = peptide), alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(time, value, colour = peptide), alpha = alpha) -> gp
  gp
}