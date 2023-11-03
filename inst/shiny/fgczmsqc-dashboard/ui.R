#R
## 
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "FCGZ MS QC"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("DIA-NN stat", tabName = "diannplots", icon = icon("chart-line")),
      menuItem("DIA-NN stat data", tabName = "dianndata", icon = icon("table")),
      menuItem("comet stat", tabName = "cometplots", icon = icon("chart-line")),
      menuItem("comet stat data", tabName = "cometdata", icon = icon("table")),
      menuItem("iRT profiles", tabName = "iRTprofiles", icon = icon("chart-line")),
      menuItem("TIC", tabName = "tic", icon = icon("chart-line")),
      htmlOutput("instrument"),
      htmlOutput("file"),
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
      href='https://gitlab.bfabric.org/proteomics/qc/-/tree/main/inst/shiny/fgczmsqc-dashboard?ref_type=heads'),
    br(),
    HTML("tested on"),
    img(src='https://upload.wikimedia.org/wikipedia/commons/2/28/Firefox_logo%2C_2017.svg', width = '30px'),
    sidebarMenu(menuItem("sessionInfo", tabName = "sessionInfo", icon = icon("laptop")))
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabItems(
      # First tab content
      tabItem(tabName = "diannplots",
              htmlOutput("variable"),
              fluidRow(htmlOutput("diannTimeSlider")),
              fluidRow(box(plotOutput("diannPlot", height = 500, width = 1000)))
      ),
      # Second tab content
      tabItem(tabName = "dianndata",
              fluidRow(
                h2("DIA-NN stat.tsv"),
                dataTableOutput('tableDIANN')
              )
      ),
      tabItem(tabName = "cometplots",
              fluidRow(htmlOutput("cometVariable")),
              fluidRow(htmlOutput("cometTimeSlider")),
              fluidRow(box(plotOutput("cometPlot",width = 1000)))
      ),
      tabItem(tabName = "cometdata",
              fluidRow(
                h2("comet RData"),
                fluidRow(DT::dataTableOutput('tableComet'))
              )
      ),
      tabItem(tabName = "iRTprofiles",
              fluidRow(box(plotOutput("plotiRTDDAChromatograms",
                                      height = 300, width = 1000))),
              fluidRow(box(plotOutput("plotDDAiRTfits", height = 300, width = 300))),
              fluidRow(box(plotOutput("plotDDAiRTprofiles", height = 400, width = 1000))),
              fluidRow(radioButtons("ppmError", "ppmError",
                                    choices = c(10, 15, 20, 100),
                                    selected = 10,
                                    inline = TRUE,
                                    width = NULL)),
              checkboxInput("plotDiannMs2", "analyse ms2 profiles", value = FALSE, width = NULL),
              fluidRow(htmlOutput("scanType")),
              fluidRow(plotOutput("plotDIAiRTprofiles", height = 400, width = 1000))
      ),
      tabItem(tabName = "tic",
              fluidRow(
                h2("Total ion count"),
                fluidRow(
                  box(plotOutput("plotTIC", height = 400, width = 1000))
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