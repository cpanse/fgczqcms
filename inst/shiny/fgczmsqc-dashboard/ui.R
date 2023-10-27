#R
## 
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "FCGZ MS QC"),
  dashboardSidebar(
    htmlOutput("instrument"),
    htmlOutput("variable"),
    htmlOutput("file"),
    radioButtons(
      "method",
      "Method",
      choices = c("DIANN", "iRT profile", "TIC"),
      selected = c("DIANN"),
      inline = FALSE,
      width = NULL
  #    choiceNames = c("DIANN"),
  #    choiceValues = c("DIANN")
    ),
  br(),
  a(img(src="https://img.shields.io/badge/JIB-10.1515%2Fjib.2022.0031-brightgreen"),
    href='https://www.degruyter.com/document/doi/10.1515/jib-2022-0031/html'),
  a(img(src="https://img.shields.io/badge/JPR-10.1021%2Facs.jproteome.0c00866-brightgreen"),
    href='http://dx.doi.org/10.1021/acs.jproteome.0c00866')
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      sliderInput("days", "Observation range in days:", min = 0, max = 365, value = c(0, 28), width = 1000),
      box(plotOutput("plot1", height = 500, width = 1000))
    )
  )
)