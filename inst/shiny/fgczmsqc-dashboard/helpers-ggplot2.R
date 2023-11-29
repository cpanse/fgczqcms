#R
## some ggplot2 templates for autoQC03 plots
.ggplot <- function(data = NULL, variables = NULL){
  stopifnot(!is.null(data),
            !is.null(variables))
  data |> 
    subset(variable %in% variables) |> 
    ggplot2::ggplot(ggplot2::aes(time, value)) +
    ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1)
}

.ggplotAutoQC03 <- function(data = NULL, variables = NULL){
  .ggplot(data, variables) +
    ggplot2::geom_point(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) -> gp
  gp
}

.ggplotAutoQC01 <- function(data = NULL, variables = NULL){
  .ggplot(data, variables) +
    ggplot2::geom_point(ggplot2::aes(time, value, colour = peptide), alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(time, value, colour = peptide), alpha = 0.4) -> gp
  gp
}