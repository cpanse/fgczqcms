#R
#' @noRd
#' @export
configServer <- function(id){
  moduleServer(id,
               function(input, output, session) {
                 rootdirraw <- reactive({
                   cands <- c("/srv/www/htdocs/", "/scratch/DIAQC/qc/dump", "/Users/cp/Downloads/dump/")
                   for (d in cands){
                     if (dir.exists(d))return(d)
                   }
                   NULL
                 })
                 
                 rootdir <- reactive({
                   cands <- c("/scratch/FGCZMSQC-DASHBOARD/", "/Users/cp/Downloads/dump/", "/scratch/DIAQC/qc/dump")
                   for (d in cands){
                     if (dir.exists(d))return(d)
                   }
                   NULL
                 })
                 return(list(rootdir = rootdir(),
                             rootdirraw = rootdirraw()))
               }
  )
}
