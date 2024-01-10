#R


bfabricInstrumentEventUI <- function(id){
  ns <- NS(id)
  
  htmlOutput(ns("ui"))
  
}

bfabricInstrumentEventServer <- function(id, filterValues){
  moduleServer(id,
               function(input, output, session) {
                 
                 output$ui <- renderUI({
                   
                   status <- function(){
                     if (filterValues$useBFabric){
                       "info"
                     }else{
                       "warning"
                     }
                   }
                   
                   tagList(
                     shinydashboard::box(title = "B-Fabric instrument events frequency",
                                         tableOutput(NS(id, "bfabricInstrumentEventsFrequency")),
                                         status = status(),
                                         solidHeader = TRUE,
                                         collapsible = TRUE,
                                         collapsed = !filterValues$useBFabric,
                                         width = 6,
                     ),
                     shinydashboard::box(title = "B-Fabric last instrument event",
                                         tableOutput(NS(id, "lastEntryBfabricInstrumentEvents")),
                                         status = status(),
                                         solidHeader = TRUE,
                                         collapsible = TRUE,
                                         collapsed = !filterValues$useBFabric,
                                         width = 6,
                                         footer = "display the last entry for each instrument."),
                     shinydashboard::box(title = "B-Fabric instrument events",
                                           plotOutput(NS(id, "plotBfabricInstrumentEvents"), height = 600),
                                         status = status(),
                                         solidHeader = TRUE,
                                         collapsible = TRUE,
                                         collapsed = !filterValues$useBFabric,
                                         width = 12,
                                         footer = "Once enabled, the graph displays all (including
                       child instruments) instrument events fetched by the
                       bfabric system.")
                     
                     
                   )
                 })
                 
                 
                 ## TODO(cpanse): move to helper function
                 ## Fetch bfabricInstrumentEventType ========
                 bfabricInstrumentEventTypeFetch <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Fetching instrument event types ...")
                   on.exit(progress$close())
                   ## TODO(cp): considering adding it as helper function
                   ## to make it possible to test upfront
                   bfabricShiny::readPages(login,
                                           webservicepassword,
                                           endpoint = 'instrumenteventtype',
                                           query = list(),
                                           posturl = bfabricposturl) |>
                     lapply(function(x){data.frame(id = x$`_id`, name=x$name)}) |>
                     Reduce(f = rbind) -> rv
                   
                   rv
                 })
                 
               
                 ## TODO(cp): give meaningful names
                 bfabricInstrumentEvents <- reactive({
                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Fetching instrument event ...")
                   on.exit(progress$close())
                   
                   S <- .composeInstrumentEvents()
                   IET <- bfabricInstrumentEventTypeFetch()
                   S$InstrumentEventTypeName <- IET[match(S$instrumenteventtypeid, IET$id), 'name']
                   return(S)
                 })
                 
                 bfabricInstrumentEventsFiltered <- reactive({
                   shiny::req(bfabricInstrumentEvents())

                   progress <- shiny::Progress$new(session = session)
                   progress$set(message = "Filter B-Fabric instrument events ...")
                   on.exit(progress$close())


                   Filter <- bfabricInstrumentEvents()$Instrument %in% filterValues$instrument &
                     filterValues$timeMin <= bfabricInstrumentEvents()$time & bfabricInstrumentEvents()$time < filterValues$timeMax
                   
                   bfabricInstrumentEvents()[Filter, ] -> rv
                   idx <- order(rv$time, decreasing = TRUE)
                   rv$time <- format(rv$time, format = "%Y-%m-%d %H:%M")

                   rv[idx, ]
                 })
                 
                 
                 output$bfabricInstrumentEventsFrequency <- renderTable({
                   if (filterValues$useBFabric){
                     shiny::req(bfabricInstrumentEvents())
                     aggregate(. ~ Instrument + type,
                               data = bfabricInstrumentEvents(),
                               FUN=length)[, 1:3] -> df
                     
                     return(df)
                   }else{
                     return(NULL)
                   }
                 })
                 
                 output$lastEntryBfabricInstrumentEvents <- renderTable({
                   
                  
                   if (filterValues$useBFabric){
                     shiny::req(bfabricInstrumentEvents())
                     
                     progress <- shiny::Progress$new(session = session)
                     progress$set(message = "determine last ", detail = "instrument event ...")
                     on.exit(progress$close())
                     rv <- .determineLastEntry(bfabricInstrumentEvents())
                     print(rv)
                     return(rv)
                   }else{NULL}
                   
                 })
                 
                 ## Plots ==================
                 output$plotBfabricInstrumentEvents <- renderPlot({
                   if (filterValues$useBFabric){
                     shiny::req(bfabricInstrumentEvents())

                     bfabricInstrumentEvents() |> .plotBfabricEvents()
                  
                   }else{NULL}
                 })
                 
                 
                 return(list(bfabricInstrumentEvents = bfabricInstrumentEvents,
                             bfabricInstrumentEventsFiltered = bfabricInstrumentEventsFiltered))
               }
  )
}

