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
                   shinydashboard::box(title = "B-Fabric instrument events",
                                       tagList(
                                         plotOutput(NS(id, "plotBfabricInstrumentEvents"), height = 600)
                                       ),
                                       status = status(),
                                       solidHeader = TRUE,
                                       collapsible = TRUE,
                                       collapsed = !filterValues$useBFabric,
                                       width = 12,
                                       footer = "Once enabled, the graph displays all (including
                       child instruments) instrument events fetched by the
                       bfabric system.")
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


                   Filter <- bfabricInstrumentEvents()$instrumentName %in% filterValues$instrument &
                     filterValues$timeMin <= bfabricInstrumentEvents()$time & bfabricInstrumentEvents()$time < filterValues$timeMax
                   
                   bfabricInstrumentEvents()[Filter, ]
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

