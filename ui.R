#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyhelper))
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
                                                   h4("Mutation Selection") %>%
                                                       helper(content = "mutation_selection"),

                                                   radioGroupButtons("mutationSelection",
                                                                     label = "", choices = mutation.selection.options, selected = "mutations.all",
                                                                     status = "primary", individual = TRUE),
                                                   pickerInput(inputId = "mutationEffect",
                                                               label = "Mutation Effect",
                                                               choices = mut_effects,
                                                               selected = mut_effects,
                                                               options = list(
                                                                   `actions-box` = TRUE
                                                               ),
                                                               multiple = TRUE
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Selection"),

                                                   selectInput("plotType", "Plot Type", plot.type.options) %>%
                                                       helper(content = "plottype"),
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.gene.plot'",
                                                       selectizeInput("gene.input", "Genes", choices = c(), multiple=TRUE,  # choices updated from the server side
                                                                      options = list(maxItems = 9,
                                                                                     plugins = list('remove_button'))) %>%  # enable deselection of items
                                                           helper(content = "gene_selection")
                                                   ),
                                                   # pathway plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.pathway.plot'",

                                                       # choice between existing and custom pathway
                                                       selectInput("pathwayType", "Pathway Definition Source", pathway.type.options) %>%
                                                           helper(content = "pathway_definition_source"),
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
                                                       sliderInput("tmb.cutoff", "Burden Cutoff", min = 1, max = 10, value = 5, step = 1) %>%
                                                           helper(content = "tmb_cutoff")
                                                   ),
                                                   # burden plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.waterfall'",
                                                       numericInput("waterfall.cutoff", "Number of most mutated genes", min = 1, max = 50, value = 20) %>%
                                                           helper(content = "waterfall_cutoff")
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Sample Selection") %>%
                                                       helper(content = "sample_selection"),

                                                   selectInput("treatment.input", "Treatment Group",
                                                               choices = treatment.options),

                                                   radioGroupButtons("hist.subtype",
                                                                     label = "Histological Subtype", choices = ui.options$HistSubtype, selected = ui.options$HistSubtype[["Any"]],
                                                                     status = "primary", individual = TRUE),
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
                                                   h4("Plot Dimensions") %>%
                                                       helper(content = "plot_dimensions"),
                                                   splitLayout(
                                                       numericInput("height",
                                                                    label = "Height",
                                                                    value = 700
                                                       ),
                                                       numericInput("width",
                                                                    label = "Width",
                                                                    value = 700
                                                       )
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Color Settings") %>%
                                                       helper(content = "color_settings"),
                                                   selectInput("color.palette", "Color Palette",
                                                               choices = color.palette.options, selected = color.palette.options[["JCO"]],
                                                               multiple = FALSE)
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Legend Settings") %>%
                                                       helper(content = "legend_settings"),
                                                   awesomeCheckbox("showLegend",
                                                                   label = "Show Legend",
                                                                   value = TRUE),
                                                   conditionalPanel(
                                                       condition = "input.showLegend == true",
                                                       textInput("legendTitle", "Legend Title", value = "strata"),
                                                       selectInput("legendLoc", "Legend Location",
                                                                   choices = plot.legend.loc.options),
                                                       conditionalPanel(
                                                           condition = "input.legendLoc == 'custom'",
                                                           splitLayout(
                                                               sliderInput("legend.coord.x", "X Coordinate", min = 0, max = 1,
                                                                           value = 0.25, step = 0.01),
                                                               sliderInput("legend.coord.y", "Y Coordinate", min = 0, max = 1,
                                                                           value = 0.25, step = 0.01)
                                                           )
                                                       )
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Features") %>%
                                                       helper(content = "plot_features"),
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
                           h3("Downloads") %>%
                               helper(content = "downloads"),
                           splitLayout(
                               downloadButton("downloadPlot", label = "Plot PDF"),
                               downloadButton("downloadSamples", label = "Sample TSV"),
                               downloadButton("downloadMutations", label = "Mutation TSV")
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
                                tabPanel("Dataset Statistics",
                                         htmlOutput("datasetStats")
                                ),
                                tabPanel("Citation and About",
                                         htmlOutput("appCiteAbout")
                                )
                    )
                )
            )
        ))
))
