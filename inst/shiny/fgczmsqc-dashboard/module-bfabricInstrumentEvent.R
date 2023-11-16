#R

bfabricInstrumentEventUI <- function(id){
  fluidRow(box(plotOutput(NS(id, "plotSummaryBfabricEvents"),
                          height = 600), width = "95%"))
}

bfabricInstrumentEventServer <- function(id, filterValues){
  moduleServer(id,
               function(input, output, session) {
                 ## Fetch bfabricInstrumentEventType
                 bfabricInstrumentEventTypeFetch <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Fetching instrument event types ...")
                   on.exit(progress$close())
                   
                   bfabricShiny::readPages(login,
                                           webservicepassword,
                                           endpoint='instrumenteventtype',
                                           query=list(),
                                           posturl=bfabricposturl) |>
                     lapply(function(x){data.frame(id=x$`_id`, name=x$name)}) |>
                     Reduce(f = rbind)
                 })
                 
                 ## Fetch bfabricInstrumentEvents  -------------
                 bfabricInstrumentFetch <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Fetching instrument events ...")
                   on.exit(progress$close())
                   
                   rv <- bfabricShiny::readPages(login,
                                                 webservicepassword,
                                                 endpoint = 'instrumentevent',
                                                 query = list(instrumentid = .getInstruments() |>
                                                                as.integer() |> as.list()),
                                                 posturl = bfabricposturl)  |>
                     lapply( FUN=function(x){
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
                 })
                 
                 bfabricInstrumentEvents <- reactive({
                   
                   
                   S <- bfabricInstrumentFetch()
                   II <- .getInstruments()
                   IET <- bfabricInstrumentEventTypeFetch()
                   
                   S$InstrumentName <- names(sapply(S$instrumentid, function(x){which(x == II)}))
                   S$InstrumentEventTypeName <- IET[match(S$instrumenteventtypeid, IET$id), 'name']
                   return(S)
                   
                 })
                 
                 bfabricInstrumentEventsFiltered <- reactive({
                   shiny::req(bfabricInstrumentEvents())
                   # shiny::req(input$instrument)
                   
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Filter B-Fabric instrument events ...")
                   on.exit(progress$close())
                   
                   now <- Sys.time()
                   
                   Filter <- as.integer(bfabricInstrumentEvents()$instrumentid) %in% (.getInstruments()[filterValues$instrument] |> unlist() |> as.integer()) &
                     filterValues$timeMin <= bfabricInstrumentEvents()$time & bfabricInstrumentEvents()$time < filterValues$timeMax
                   
                   bfabricInstrumentEvents()[Filter, ]
                 })
                 
                 ## Plots ==================
                 output$plotSummaryBfabricEvents <- renderPlot({
                   #shiny::req(input$useBfabric)
                   shiny::req(bfabricInstrumentEvents())
                   
                   n <- length(unique(bfabricInstrumentEvents()$instrumentid))
                   
                   lattice::dotplot(~ time | InstrumentName,
                                    group = InstrumentEventTypeName,
                                    layout = c(1, n),
                                    data = bfabricInstrumentEvents(),
                                    cex = 1,
                                    pch = 22,
                                    auto.key = list(space = "right", pch=22),
                                    main = 'B-Fabric instrument events grouped by event type',
                   )
                 })
                 
                 
                 
                 return(list(bfabricInstrumentEvents = bfabricInstrumentEvents,
                             bfabricInstrumentEventsFiltered = bfabricInstrumentEventsFiltered))
               }
  )
}

