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

# Define UI for application that draws a histogram
shinyUI(fluidPage(
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

            # Sidebar with plot settings
            sidebarPanel(
                tabsetPanel(type = "pills",
                            tabPanel("Data Selection",
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
                                     ),
                                     selectInput("treatment.input", "Treatment Group",
                                                 choices = treatment.options),

                                     radioGroupButtons("radio.er.status",
                                                       label = "ER Status", choices = ui.options$ER, selected = ui.options$ER[["Any"]],
                                                       status = "primary", individual = TRUE),
                                     radioGroupButtons("radio.pgr.status",
                                                       label = "PgR Status", choices = ui.options$PgR, selected = ui.options$PgR[["Any"]],
                                                       status = "primary", individual = TRUE),
                                     radioGroupButtons("radio.her2.status",
                                                       label = "HER2 Status", choices = ui.options$HER2, selected = ui.options$HER2[["Any"]],
                                                       status = "primary", individual = TRUE),
                                     radioGroupButtons("radio.ki67.status",
                                                       label = "Ki67 Status", choices = ui.options$Ki67, selected = ui.options$Ki67[["Any"]],
                                                       status = "primary", individual = TRUE),
                                     radioGroupButtons("radio.nhg",
                                                       label = "NHG", choices = ui.options$NHG, selected = ui.options$NHG[["Any"]],
                                                       status = "primary", individual = TRUE),
                                     pickerInput(inputId = "radio.pam50",
                                                 label = "PAM50 Subtype",
                                                 choices = ui.options$PAM50, choicesOpt = list(content = pam50.html.labels))
                            ),
                            tabPanel("Plot Settings",
                                     numericInput("height",
                                                  label = "Plot Height",
                                                  value = 500),
                                     numericInput("width",
                                                  label = "Plot width",
                                                  value = 500),
                                     selectInput("color.palette", "Color Palette",
                                                 choices = color.palette.options, selected = color.palette.options[["JCO"]],
                                                 multiple = FALSE),
                                     checkboxInput("showLegend",
                                                   label = "Show legend",
                                                   value = TRUE),
                                     conditionalPanel(
                                         condition = "input.showLegend == true",
                                         selectInput("legendLoc", "Legend Location",
                                                     choices = plot.legend.loc.options),
                                         conditionalPanel(
                                             condition = "input.legendLoc == 'custom'",
                                             fluidRow(
                                                 column(6,
                                                        numericInput("legend.coord.x", "X Coordinate (0-1)", 0.25)),
                                                 column(6,
                                                        numericInput("legend.coord.y", "Y Coordinate (0-1)", 0.25))
                                             )
                                         )
                                     ),
                                     checkboxInput("show.risk.table",
                                                   label = "Show risk table",
                                                   value = TRUE),
                                     checkboxInput("show.pval",
                                                   label = "Show logrank p-value",
                                                   value = TRUE),
                                     checkboxInput("show.censors",
                                                   label = "Show censors",
                                                   value = FALSE),
                                     checkboxInput("show.conf.int",
                                                   label = "Show confidence intervals",
                                                   value = FALSE)
                            )
                ),
                column(12,
                       h3("Save the plot")),
                ## Code is prepared for a selectInput with the option formats, but users
                ## requested these 3 buttons as they find it easier.
                ## Harcoding fixes problem with plotting outdated data,
                ## but needs to be studied.
                column(12,
                       p("Select the format in which you wish to save the generated plot."),
                       column(4,
                              downloadButton("png",
                                             label = "png")),
                       column(4,
                              downloadButton("tiff",
                                             label = "tiff")),
                       column(4,
                              downloadButton("pdf",
                                             label = "pdf"))),
                ## Clearfix
                tags$div(class = 'clearfix')
            ),

            # Survival plot
            mainPanel(
                tabsetPanel(type = "pills",
                            tabPanel("Plots",
                                     withSpinner(plotOutput("survplot"))
                            ),
                            tabPanel("Sample Table",
                                     withSpinner(DT::dataTableOutput("sample.table"))
                            ),
                            tabPanel("Mutation Table",
                                     withSpinner(DT::dataTableOutput("mut.table"))
                            ),
                            tabPanel("Citation and About",
                                     htmlOutput("appCiteAbout"))
                )
            )
        )
    )
))
