#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))
options(spinner.type=5)
suppressPackageStartupMessages(library(DT))
source("R/resources.R")


appLoadCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"

wellpanel.settings.style = "background: white; margin-top: 15px; margin-bottom: 0px; padding-top: 3px; padding-bottom: 3px;"

# Define UI for application that draws a histogram
shinyUI(fluidPage(title = "SCAN-B Mutation Explorer",
    useShinyjs(),
    inlineCSS(appLoadCSS),

    # Loading message
    div(id = "loading-content",
        h2("Loading...")
    ),

    hidden(
        div(id = "app-content",

            # Application title, rendered by the server
            uiOutput("header_panel"),

            sidebarLayout(
                sidebarPanel(
                    tabsetPanel(type = "pills",
                                tabPanel("Plot and Data Selection",
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Selection", style = "padding-top: 5px;"),

                                                   selectInput("plotType", "Plot Type", plot.type.options),
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.gene.plot'",
                                                       selectizeInput("gene.input", "Gene", choices = c(), multiple=TRUE,  # choices updated from the server side
                                                                      options = list(maxItems = 9,
                                                                                     plugins = list('remove_button')))  # enable deselection of items
                                                   ),
                                                   # pathway plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.pathway.plot'",

                                                       # choice between existing and custom pathway
                                                       selectInput("pathwayType", "Pathway Definition Source", pathway.type.options),
                                                       conditionalPanel(
                                                           condition = "input.pathwayType == 'pathway.reactome'",
                                                           selectInput("pathway.input", "Reactome",
                                                                       choices = pathway.ui.options,
                                                                       multiple=FALSE)
                                                       ),
                                                       conditionalPanel(
                                                           condition = "input.pathwayType == 'pathway.custom'",
                                                           selectInput("custom.pathway.input", "Custom Definition",
                                                                       choices = c(),  # updated from the server side
                                                                       multiple=TRUE)
                                                       )
                                                   ),
                                                   # burden plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.burden.plot'",
                                                       numericInput("tmb.cutoff", "", 75)  # updated from the server side
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Sample Selection"),

                                                   selectInput("treatment.input", "Treatment Group",
                                                               choices = treatment.options),

                                                   radioGroupButtons("er.status",
                                                                     label = "ER Status", choices = ui.options$ER, selected = ui.options$ER[["Any"]],
                                                                     status = "primary", individual = TRUE),
                                                   radioGroupButtons("pgr.status",
                                                                     label = "PgR Status", choices = ui.options$PgR, selected = ui.options$PgR[["Any"]],
                                                                     status = "primary", individual = TRUE),
                                                   radioGroupButtons("her2.status",
                                                                     label = "HER2 Status", choices = ui.options$HER2, selected = ui.options$HER2[["Any"]],
                                                                     status = "primary", individual = TRUE),
                                                   radioGroupButtons("ki67.status",
                                                                     label = "Ki67 Status", choices = ui.options$Ki67, selected = ui.options$Ki67[["Any"]],
                                                                     status = "primary", individual = TRUE),
                                                   radioGroupButtons("nhg",
                                                                     label = "NHG", choices = ui.options$NHG, selected = ui.options$NHG[["Any"]],
                                                                     status = "primary", individual = TRUE),
                                                   pickerInput(inputId = "pam50",
                                                               label = "PAM50 Subtype",
                                                               choices = ui.options$PAM50, choicesOpt = list(content = pam50.html.labels))
                                         )
                                ),
                                tabPanel("Plot Settings",
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Dimensions"),
                                                   splitLayout(
                                                       numericInput("height",
                                                                    label = "Height",
                                                                    value = 500
                                                       ),
                                                       numericInput("width",
                                                                    label = "Width",
                                                                    value = 500
                                                       )
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Color Settings"),
                                                   selectInput("color.palette", "Color Palette",
                                                               choices = color.palette.options, selected = color.palette.options[["JCO"]],
                                                               multiple = FALSE)
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Legend Settings"),
                                                   awesomeCheckbox("showLegend",
                                                                   label = "Show legend",
                                                                   value = TRUE),
                                                   conditionalPanel(
                                                       condition = "input.showLegend == true",
                                                       selectInput("legendLoc", "Legend Location",
                                                                   choices = plot.legend.loc.options),
                                                       conditionalPanel(
                                                           condition = "input.legendLoc == 'custom'",
                                                           splitLayout(
                                                               numericInput("legend.coord.x", "X Coord (0-1)", 0.25),
                                                               numericInput("legend.coord.y", "Y Coord (0-1)", 0.25)
                                                           )
                                                       )
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Features"),
                                                   awesomeCheckbox("show.risk.table",
                                                                   label = "Show risk table",
                                                                   value = TRUE),
                                                   awesomeCheckbox("show.pval",
                                                                   label = "Show logrank p-value",
                                                                   value = TRUE),
                                                   awesomeCheckbox("show.censors",
                                                                   label = "Show censors",
                                                                   value = FALSE),
                                                   awesomeCheckbox("show.conf.int",
                                                                   label = "Show confidence intervals",
                                                                   value = FALSE)
                                         )
                                )
                    ),
                    wellPanel(style = wellpanel.settings.style,
                           h3("Download"),
                           splitLayout(
                               downloadButton("downloadPlot", label = "Plot as PDF")
                           )
                    )
                ),

                # Survival plot
                mainPanel(
                    tabsetPanel(type = "pills",
                                tabPanel("Plots",
                                         withSpinner(plotOutput("survplot"))
                                ),
                                tabPanel("Sample Table",
                                         DT::dataTableOutput("sample.table")
                                ),
                                tabPanel("Mutation Table",
                                         DT::dataTableOutput("mut.table")
                                ),
                                tabPanel("Citation and About",
                                         htmlOutput("appCiteAbout"))
                    )
                )
            )
        ))
))
