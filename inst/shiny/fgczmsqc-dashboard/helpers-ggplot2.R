#R


.ggplotAutoQC03 <- function(data = NULL, variables = NULL, useBFabric = FALSE, BFabricTime = NULL){
  stopifnot(!is.null(data),
            !is.null(variables))
  
  data |> 
    subset(variable %in% variables) |> 
    ggplot2::ggplot(ggplot2::aes(time, value)) +
    ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1)  -> gp
  
  if ("scanType" %in% names(data)){
    gp +
      ggplot2::geom_point(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) +
      ggplot2::geom_line(ggplot2::aes(time, value, colour = scanType), alpha = 0.4) -> gp
    
  }else{
    gp +
      ggplot2::geom_point(ggplot2::aes(time, value), alpha = 0.4) +
      ggplot2::geom_line(ggplot2::aes(time, value), alpha = 0.4) -> gp
  }
  if (useBFabric & !is.null(BFabricTime)){
    gp + ggplot2::geom_vline(xintercept = BFabricTime,
                             linetype="dashed", 
                             color = "red", size = 1) -> gp
  }
  gp
}


#R


.ggplotAutoQC01 <- function(data = NULL, variables = NULL, useBFabric = FALSE, BFabricTime = NULL){
  stopifnot(!is.null(data),
            !is.null(variables))

  data |> 
    subset(variable %in% variables) |> 
    ggplot2::ggplot(ggplot2::aes(time, value)) +
    ggplot2::facet_wrap(. ~  Instrument * variable, scales="free_y", ncol = 1) +
    ggplot2::geom_point(ggplot2::aes(time, value, colour = peptide), alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(time, value, colour = peptide), alpha = 0.4) -> gp

  if (useBFabric & !is.null(BFabricTime)){
    gp + ggplot2::geom_vline(xintercept = BFabricTime,
                             linetype="dashed", 
                             color = "red", size = 1) -> gp
  }
  gp
}