#R
## 
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "FCGZ MS QC"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("autoQC01", tabName = "autoQC01Plots", icon = icon("chart-line"), badgeLabel = "new", badgeColor = "green"),
      menuItem("DIA-NN stat", tabName = "diannplots", icon = icon("chart-line")),
      menuItem("DIA-NN stat data", tabName = "dianndata", icon = icon("table")),
      menuItem("DDA-comet stat", tabName = "cometplots", icon = icon("chart-line")),
      menuItem("DDA-comet stat data", tabName = "cometdata", icon = icon("table")),
      hr(),
      menuItem("raw file", tabName = "rawFile", icon = icon("chart-line")),
      menuItem("TIC", tabName = "tic", icon = icon("chart-line")),
      htmlOutput("instrument"),
      
      selectInput('regex', 'file regex',
                  c(".*", ".*raw$", ".*autoQC.*dia.*raw$", ".*autoQC.*dda.*raw$", "*.zip"),
                  multiple = FALSE,
                  selected = ".*raw$")
      ),
    br(),
    a(img(src="https://img.shields.io/badge/JIB-10.1515%2Fjib.2022.0031-brightgreen"),
      href='https://www.degruyter.com/document/doi/10.1515/jib-2022-0031/html'),
    a(img(src="https://img.shields.io/badge/JPR-10.1021%2Facs.jproteome.0c00866-brightgreen"),
      href='http://dx.doi.org/10.1021/acs.jproteome.0c00866'),
    a(img(src="https://img.shields.io/badge/git-proteomics/fgczmsqc-pink"),
      href='https://gitlab.bfabric.org/proteomics/qc/'),
    br(),
    HTML("tested on"),
    img(src='https://upload.wikimedia.org/wikipedia/commons/2/28/Firefox_logo%2C_2017.svg', width = '30px'),
    sidebarMenu(menuItem("sessionInfo", tabName = "sessionInfo", icon = icon("laptop")))
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabItems(
      tabItem(tabName = "autoQC01Plots",
              fluidRow(h2("autoQC01")),
              fluidRow(htmlOutput("autoQC01Variable")),
              fluidRow(htmlOutput("autoQC01TimeSlider")),
              #fluidRow(htmlOutput("autoQC01TimeSlider")),
              fluidRow(box(plotOutput("autoQC01Plot"), height = "55%", width = "100%"))
      ),
      tabItem(tabName = "diannplots",
              htmlOutput("variable"),
              fluidRow(htmlOutput("diannTimeSlider")),
              fluidRow(box(plotOutput("diannPlot"), height = "75%", width = "100%"))
      ),
      tabItem(tabName = "dianndata",
              fluidRow(
                h2("DIA-NN stat.tsv"),
                DT::dataTableOutput('tableDIANN')
                #htmlOutput("bfabricInstrumentEventsOutput")
              )
      ),
      tabItem(tabName = "cometplots",
              fluidRow(htmlOutput("cometVariable")),
              fluidRow(htmlOutput("cometTimeSlider")),
              fluidRow(box(plotOutput("cometPlot"), height = "75%", width = "100%")),
              hr()
              fluidRow(box(htmlOutput("bfabricInstrumentEventsOutput")))
      ),
      tabItem(tabName = "cometdata",
              fluidRow(
                h2("comet RData"),
                fluidRow(DT::dataTableOutput('tableComet'))
              )
      ),
      tabItem(tabName = "rawFile",
              tagList(
            htmlOutput("fileInput", width = "100%"),
            htmlOutput("fileOutput"))
      ),
      tabItem(tabName = "tic",
              fluidRow(
                h2("Total ion count"),
                fluidRow(box(htmlOutput("ticFileInput"), width = "100%")),
                fluidRow(
                  box(plotOutput("plotTIC"), width = "100%")
                )
              )
      ),
      tabItem(tabName = "sessionInfo",
              fluidRow(
                h2("session information"),
                fluidRow(
                  box(verbatimTextOutput("sessionInfo"), width = 800)
                )
              )
      )
    )
  )
)