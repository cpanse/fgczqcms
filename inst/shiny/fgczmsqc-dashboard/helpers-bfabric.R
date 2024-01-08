#R

.getChildInstruments <- function(instrumentIds = (.getInstruments() |>
                                   as.integer() |> as.list())){
  
  rv <- bfabricShiny::readPages(login,
                                webservicepassword,
                                endpoint = 'instrument',
                                query = list(id = instrumentIds),
                                posturl = bfabricposturl) -> rv
  
  rv |>
    lapply(function(x){
      x$child |>
        vapply(function(y)y$`_id`, 1) -> child;
      list(parent=x$`_id`, child=child)
    }) |>
    lapply(function(x){
      ## message("length = ", length(x$child))
      if(length(x$child) >0){
        data.frame(parent = rep(x$parent, length(x$child)), child = x$child)
      }else{NULL}
    }) |> Reduce(f = rbind) -> rv
  
  rv
}

#' @example 
#' .getInstruments() |> unlist()|>as.integer()|>Filter(f = function(x)x>0) -> instrumentIds
#' instrumentIds |> .getChildInstruments() -> rv.childs
#' c(rv.childs$child, instrumentIds)|> unique() |> .getInstrumentEvent() -> instrumentEvents
#' 
.getInstrumentEvent <- function(instrumentIds = NULL){
  
  rv <- bfabricShiny::readPages(login,
                                webservicepassword,
                                endpoint = 'instrumentevent',
                                query = list(instrumentid = instrumentIds),
                                posturl = bfabricposturl)  |>
    lapply(FUN = function(x){
      if (all(c('description', 'datetime') %in% names(x))){
        df <- data.frame(time = (x$datetime |> as.POSIXlt()),
                         instrumentid = as.integer(x$instrument$`_id`),
                         description  = x$description, # (x$description |> gsub(pattern = '\r\n', replacement = '')),
                         instrumenteventtypeid = x$instrumenteventtype$`_id`)
        return(df)
      }
      NULL
    }) |>
    Reduce(f = rbind)  
  
  return(rv[order(rv$time), ])
}


#' @examples
#' 
#' 
.composeInstrumentEvents <- function(){
  .getInstruments() |>
    unlist()|> as.integer() |> Filter(f = function(x)x>0) -> instrumentIds
  
  instrumentIds |>
    .getChildInstruments() -> rv.childs
  
  c(rv.childs$child, instrumentIds) |> unique() |>
    .getInstrumentEvent() -> rv
  
  
  rv$type[rv$instrumentid %in% rv.childs$child] <- 'child'
  rv$type[rv$instrumentid %in% rv.childs$parent] <- 'parent'
  
  t1 <- data.frame(instrumentName = names(.getInstruments()), instrumentId = .getInstruments()|>unlist()) 
  t2 <- merge(rv.childs, t1, by.x = 'parent', by.y = 'instrumentId')
  
  instrumentId <- c(t2$child, t1$instrumentId)
  instrumentName <- c(t2$instrumentName, t1$instrumentName)
  t3 <- data.frame(instrumentName = instrumentName , instrumentid = instrumentId)
  
  
  merge(rv, t3, by = 'instrumentid') -> rv
  
  rv
}


#' @examples
#' .composeInstrumentEvents() |> .plotBfabricEvents()
.plotBfabricEvents <- function(x){
  
  n <- x$instrumentName |> unique() |> length()
  
  lattice::dotplot(~ time | instrumentName,
                   group = x$type,
                   layout = c(1, n),
                   data = x,
                   cex = 1,
                   pch = 22,
                   auto.key = list(space = "right", pch=22),
                   sub = 'B-Fabric instrument events grouped by child|parent instrument type')
}
