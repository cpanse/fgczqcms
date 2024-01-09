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
                                posturl = bfabricposturl)  
  
  instrumentId <- vapply(rv, function(x)x$instrument$`_id`, FUN.VALUE = 1)
  instrumentEventTypeId <- vapply(rv, function(x)x$instrumenteventtype$`_id`, FUN.VALUE = 1)
  description <- vapply(rv, function(x){if('description' %in% names(x)){x$description}else{""}}, FUN.VALUE = "OTC")
  time <- lapply(rv, function(x){if('datetime' %in% names(x)){x$datetime}else{NA}}) |> unlist() |> as.POSIXct()
  
  
  df <- data.frame(
    instrumentid = instrumentId,
    instrumenteventtypeid = instrumentEventTypeId,
    description = description,
    time = time
  )
  na.omit(df) -> df
  
  df[order(df$time), ]
}


#' @examples
#' 
#' 
.composeInstrumentEvents <- function(){
  .getInstruments() |>
    unlist()|> as.integer() |> Filter(f = function(x)x>0) -> instrumentIds
  
  instrumentIds |>
    .getChildInstruments() -> rv.childs
  
  c(rv.childs$child, instrumentIds) |>
    unique() |>
    .getInstrumentEvent() -> rv
  
  
  rv$type[rv$instrumentid %in% rv.childs$child] <- 'child'
  rv$type[rv$instrumentid %in% rv.childs$parent] <- 'parent'
  
  t1 <- data.frame(instrumentName = names(.getInstruments()), instrumentId = .getInstruments()|>unlist()) 
  t2 <- merge(rv.childs, t1, by.x = 'parent', by.y = 'instrumentId')
  
  instrumentId <- c(t2$child, t1$instrumentId)
  instrumentName <- c(t2$instrumentName, t1$instrumentName)
  t3 <- data.frame(instrumentName = instrumentName , instrumentid = instrumentId)
  
  
  merge(rv, t3, by = 'instrumentid') -> rv
  colnames(rv)[colnames(rv) == "instrumentName"] <- "Instrument"
  rv
}


#' @examples
#' .composeInstrumentEvents() |> .plotBfabricEvents()
.plotBfabricEvents <- function(x){
  
  n <- x$Instrument |> unique() |> length()
  
  lattice::dotplot(~ time | Instrument,
                   group = x$type,
                   layout = c(1, n),
                   data = x,
                   cex = 1,
                   pch = 22,
                   auto.key = list(space = "right", pch=22),
                   sub = 'B-Fabric instrument events grouped by child|parent instrument type')
}
