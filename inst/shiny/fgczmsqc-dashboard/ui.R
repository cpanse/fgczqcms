#R
## Christian Panse <cp@fgcz.ethz.ch> November 2023
stopifnot(require(shinydashboard))

source('module-bfabricInstrumentEvent.R')
source('module-autoQC03.R')

dashboardPage(
  dashboardHeader(title = "FGCZ MS QC"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("autoQC03", tabName = "autoQC03", icon = icon("chart-line"),
               badgeLabel = "alpha", badgeColor = "fuchsia"),
      menuItem("autoQC01", tabName = "autoQC01", icon = icon("chart-line")),
      menuItem("Summary", tabName = "summary", icon = icon("table")),
      hr(),
      htmlOutput("instrument"),
      htmlOutput("useBfabric"),
      selectInput('timeRange', 'time range in days',
                  c(7, 14, 30, 60, 90, 180, 365, 2:10*365),
                  multiple = FALSE,
                  selected = 90),
      selectInput('regex', 'file regex',
                  c(".*", ".*raw$", ".*autoQC.*dia.*raw$", ".*autoQC.*dda.*raw$", "*.zip"),
                  multiple = FALSE,
                  selected = ".*raw$")
    ), # sidebarMenu
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
    sidebarMenu(
      menuItem("sessionInfo", tabName = "sessionInfo", icon = icon("laptop"))
    )
  ), # dashboardSidebar
  dashboardBody(
    tabItems(
      tabItem(tabName = "autoQC01",
              fluidRow(h2("autoQC01 - Biognosys iRT peptides runs")),
              fluidRow(htmlOutput("autoQC01TimeSlider")), 
              # TODO(cp): can't we call the autoQC01 right from here?
              # autoQC01UI("autoQC01")
              fluidRow(htmlOutput("autoQC01"), width = "100%"),
      ),
      tabItem(tabName = "autoQC03",
              tagList(
                fluidRow(h2("autoQC03 - runs")),
                htmlOutput("cometTimeSlider"),
                autoQC03UI("autoQC03-DDA"),
                autoQC03UI("autoQC03-DIA")
              )
      ),
      tabItem(tabName = "summary",
              fluidRow(
                h2("Summary"),
                shinydashboard::box(title = "Frequency of MS QC events",
                                    verbatimTextOutput("summary"),
                                    solidHeader = TRUE,
                                    collapsible = TRUE,
                                    status = "primary",
                                    width = 12,
                                    footer = "MS QC event frequency in numbers and graphics"),
                htmlOutput("plotSummaryLCMSruns"),
                shinydashboard::box(plotOutput("plotSummaryCumsum", height = 250),
                                             width = 12,
                                             status = "primary",
                                             solidHeader = TRUE,
                                             collapsible = TRUE,),
                bfabricInstrumentEventUI("bfabric01"),
              )
      ),
      tabItem(tabName = "sessionInfo",
              fluidRow(
                h2("Session information"),
                fluidRow(
                  shinydashboard::box(tagList(
                    verbatimTextOutput("consolenNodename"),
                    verbatimTextOutput("sessionInfo")
                  ), width = 12)
                )
              )
      ) # tabItem 
    ) # tabItems
  ) # dashboardBody
) # dashboardPage
