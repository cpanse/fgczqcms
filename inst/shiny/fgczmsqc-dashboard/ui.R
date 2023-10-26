#R
## 
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "FCGZ MS QC"),
  dashboardSidebar(
    htmlOutput("instrument"),
    htmlOutput("variable")
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      sliderInput("days", "Observation range in days:", min = 0, max = 365, value = c(0, 28), width = 1000),
      box(plotOutput("plot1", height = 500, width = 1000))
    )
  )
)