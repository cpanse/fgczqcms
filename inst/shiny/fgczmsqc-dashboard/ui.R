#R
## 
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "FCGZ MS QC"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("DIA-NN stat", tabName = "diannplots", icon = icon("chart-line")),
      menuItem("DIA-NN stat data", tabName = "dianndata", icon = icon("table")),
      menuItem("iRT profiles", tabName = "iRTprofiles", icon = icon("chart-line")),
      menuItem("TIC", tabName = "tic", icon = icon("chart-line")),
      htmlOutput("instrument"),
      htmlOutput("variable"),
      htmlOutput("file")),
    br(),
    a(img(src="https://img.shields.io/badge/JIB-10.1515%2Fjib.2022.0031-brightgreen"),
      href='https://www.degruyter.com/document/doi/10.1515/jib-2022-0031/html'),
    a(img(src="https://img.shields.io/badge/JPR-10.1021%2Facs.jproteome.0c00866-brightgreen"),
      href='http://dx.doi.org/10.1021/acs.jproteome.0c00866')
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabItems(
      # First tab content
      tabItem(tabName = "diannplots",
              fluidRow(
                sliderInput("days", "Observation range in days:", min = 0, max = 365, value = c(0, 28), width = 1000),
                box(plotOutput("plot1", height = 500, width = 1000))
              )
      ),
      # Second tab content
      tabItem(tabName = "dianndata",
              fluidRow(
                h2("DIA-NN stat.tsv"),
                dataTableOutput('table')
              )
      ),
      tabItem(tabName = "iRTprofiles",
              fluidRow(
                box(plotOutput("plotiRTDDAChromatograms", height = 300, width = 1000),
                )
                
              ),
              fluidRow(
                
                box(plotOutput("plotiRTfits", height = 300, width = 300))
              ),
              fluidRow(
               
                  box(plotOutput("plotiRTprofiles", height = 400, width = 1000))
                
              ),
              fluidRow(radioButtons("ppmError", "ppmError",
                                    choices = c(10, 15, 20, 100),
                                    selected = 10,
                                    inline = TRUE,
                                    width = NULL
              ))
      )
    )
  )
)
