#R
## Christian Panse <cp@fgcz.ethz.ch> November 2023
library(shinylogs)
library(shinydashboard)
use_tracking()
library(fgczqcms)
imgBanner <- "graphics/fgcz-header-background.png"
imgBanner <- "http://fgcz-ms.uzh.ch/~cpanse/fgcz-cropped.svg"
# imgBanner <- "/Users/cp/src/gitlab.bfabric.org/proteomics/qc/inst/shiny/fgczmsqc-dashboard/graphics/fgcz-header-background.png"
#stopifnot(file.exists(imgBanner))

tl <- tagList(
  tags$li(
    a(href = 'http://www.fgcz.ch', 
      target = "_blank",
      img(src = imgBanner, title = "FGCZ", height = "30px"),
      style = "padding-top:10px; padding-bottom:5px;"),
    class = "dropdown"),
)


dashboardPage(
  skin = "black",
  dashboardHeader(title = paste0("fgczqcms v",
                                 packageVersion('fgczqcms')), .list = tl),
  dashboardSidebar(
    sidebarMenu(
      menuItem("autoQC01", tabName = "autoQC01beta", icon = icon("chart-line"),
        badgeLabel = "ready", badgeColor = "green"),
      menuItem("autoQC03 - DDA | DIA", tabName = "autoQC03", icon = icon("chart-line"),
               badgeLabel = "now with tims data", badgeColor = "yellow"),
      menuItem("autoQC01", tabName = "autoQC01", icon = icon("chart-line"),
                 badgeLabel = "deprecated module", badgeColor = 'blue'),
      menuItem("summary | status", tabName = "summary", icon = icon("table"),
               badgeLabel = "ready", badgeColor = "green"),
      hr(),
      selectInput('instrument', 'instruments',
                  names(.getInstruments()),
                  multiple = FALSE,
                  selected = names(.getInstruments())[2]),
      htmlOutput("useBfabric"),
      selectInput('timeRange', 'time range in days',
                  c(7, 14, 30, 60, 90, 180, 365, 2:10*365),
                  multiple = FALSE,
                  selected = 90),
      checkboxInput("addSmoothing", "add smoothing", value = FALSE),
      selectInput("heightPx", "height in px",
                  c(250, 300, 350, 400, 450, 500),
                  multiple = FALSE,
                  selected = 400),
      selectInput("peptide", "peptide", 
                  names(.iRTmz()), 
                  selected = c('LGGNEQVTR', 'YILAGVENSK', 'DGLDAASYYAPVR',
                               'GTFIIDPAAVIR'),
                  multiple = TRUE),
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
    a(img(src="https://img.shields.io/badge/JPR-10.1021%2Facs.jproteome.8b00173-brightgreen"),
      href='http://dx.doi.org/10.1021/acs.jproteome.8b00173'),
    a(img(src="https://img.shields.io/badge/git-proteomics/fgczqcms-pink"),
      href='https://gitlab.bfabric.org/proteomics/qc/'),
    br(),
    HTML("tested on"),
    img(src='https://upload.wikimedia.org/wikipedia/commons/2/28/Firefox_logo%2C_2017.svg', width = '30px'),
    sidebarMenu(
      tagList(
      menuItem("session | debug", tabName = "sessionInfo", icon = icon("laptop")),
      menuItem("readme", tabName = "README", icon = icon("laptop")))
    )
  ), # dashboardSidebar
  dashboardBody(
    tabItems(
      tabItem(tabName = "autoQC01beta",
              tagList(
                htmlOutput("autoQC01TimeSlider"),
                autoQC03UI("__autoQC01__"),
              )
      ),
      tabItem(tabName = "autoQC03",
              tagList(
                htmlOutput("cometTimeSlider"),
                autoQC03UI("autoQC03-DDA"),
                autoQC03UI("autoQC03-DIA")
              )
      ),
      tabItem(tabName = "autoQC01",
              fluidRow(h2("autoQC01 - Biognosys iRT peptides runs (old)")),
              # TODO(cp): can't we call the autoQC01 right from here?
              # autoQC01UI("autoQC01")
              fluidRow(htmlOutput("autoQC01"), width = "100%"),
      ),
      tabItem(tabName = "summary",
              fluidRow(
                h2("Summary | Status"),
                shinydashboard::box(title = "Last entries",
                                    status = "primary",
                                    solidHeader = TRUE,
                                    width = 12,
                                    collapsible = TRUE,
                                    footer = "Entries are reverse sorted by date.",
                                    tagList(
                                      shinydashboard::box(title = "autoQC01",
                                                          tableOutput('lastEntryAutoQC01'),
                                                          status = "primary",
                                                          solidHeader = TRUE,
                                                          width = 4,
                                                          collapsible = FALSE),
                                      shinydashboard::box(title = "autoQC03 DDA",
                                                          tableOutput('lastEntryAutoQC03dda'),
                                                          status = "primary",
                                                          solidHeader = TRUE,
                                                          width = 4,
                                                          collapsible = FALSE),
                                      shinydashboard::box(title = "autoQC03 DIA",
                                                          tableOutput('lastEntryAutoQC03dia'),
                                                          status = "primary",
                                                          solidHeader = TRUE,
                                                          width = 4,
                                                          collapsible = FALSE))),
              #  shinydashboard::box(title = "Frequency of MS QC events",
              #                      htmlOutput("summaryFrequency"),
              #                      solidHeader = TRUE,
              #                      collapsible = TRUE,
              #                      status = "primary",
              #                      width = 12,
              #                      footer = "MS QC event frequency in numbers and graphics"),
                shinydashboard::box(title = "Overview: Instrument ~ time | method",
                                    plotOutput("plotSummaryLCMSruns", height = 900),
                                    width = 12,
                                    status = "primary",
                                    solidHeader = TRUE,
                                    collapsed = TRUE,
                                    collapsible = TRUE,
                                    footer = "Each dot represents an LC-MS QC run."),
              shinydashboard::box(title = "Frequency",
                                  plotOutput("plotFrequency", height = 1500),
                                  width = 12,
                                  status = "primary",
                                  solidHeader = TRUE,
                                  collapsed = TRUE,
                                  collapsible = TRUE,
                                  footer = "Graphs monthly number of LC-MS QC runs for each instrument and method."),
                bfabricInstrumentEventUI("bfabric01"),
              )
      ),
      tabItem(tabName = "README",
                  shinydashboard::box(tagList(
                    includeMarkdown("README.md"),
                  ), width = 12)
      ),
      tabItem(tabName = "sessionInfo",
              fluidRow(
                h2("Session | debug - information"),
                fluidRow(
                  shinydashboard::box(tagList(
                    
                    verbatimTextOutput("sessionInfo", placeholder = TRUE),
                    verbatimTextOutput("clientDataText"),
                  ), width = 12)
                )
              )
      ) # tabItem 
    ) # tabItems
  )
) # dashboardPage
