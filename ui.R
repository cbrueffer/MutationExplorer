#
# This is the user-interface definition of the web application.
#

suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))
options(spinner.type=5)


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

tabControlJS <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

tabControlCSS <- "
  .nav li a.disabled {
  background-color: unset !important;
  color: grey !important;
  cursor: not-allowed !important;
}"


wellpanel.settings.style = "background: white; margin-top: 15px; margin-bottom: 0px; padding-top: 3px; padding-bottom: 3px;"

# Define UI for application that draws a histogram
shinyUI(fluidPage(title = "SCAN-B MutationExplorer",
    useShinyjs(),
    extendShinyjs(text = tabControlJS, functions = c("enableTab", "disableTab")),
    inlineCSS(c(tabControlCSS, appLoadCSS)),

    # Enable Google Analytics
    tags$head(includeScript("google-analytics.js")),

    # Loading message
    div(id = "loading-content",
        h2("Loading...")
    ),

    hidden(
        div(id = "app-content",

            titlePanel(h1(config$app_title,
                          h3(paste(prettyNum(n.mut, big.mark=","),
                                   "mutations in",
                                   prettyNum(n.samples, big.mark=","),
                                   "Primary Breast Cancer Transcriptomes")))
            ),
            sidebarLayout(
                sidebarPanel(
                    tabsetPanel(type = "pills",
                                tabPanel("Plot and Data Selection",
                                         wellPanel(style = wellpanel.settings.style,
                                                   h4("Plot Selection"),

                                                   helper(selectInput("plotType", "Plot Type", plot.type.options),
                                                          content = "plot_type"),
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.gene.plot'",
                                                       helper(selectizeInput("gene.input", "Genes", choices = names(mutated.gene.columns), multiple=TRUE,
                                                                             options = list(maxItems = 9,
                                                                                            plugins = list('remove_button'))),  # enable deselection of items
                                                              content = "gene_selection")
                                                   ),
                                                   # pathway plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.pathway.plot'",

                                                       # choice between existing and custom pathway
                                                       helper(selectInput("pathwayType", "Pathway Definition Source", pathway.type.options),
                                                              content = "pathway_definition_source"),
                                                       conditionalPanel(
                                                           condition = "input.pathwayType == 'pathway.reactome'",
                                                           selectizeInput("pathway.input", "Reactome Pathways",
                                                                          choices = pathway.ui.options,
                                                                          multiple = TRUE,
                                                                          options = list(maxItems = 9,
                                                                                         plugins = list('remove_button')))  # enable deselection of items
                                                       ),
                                                       conditionalPanel(
                                                           condition = "input.pathwayType == 'pathway.custom'",
                                                           selectizeInput("custom.pathway.input", "Genes in Custom Pathway (1-20)",
                                                                          choices = names(mutated.gene.columns),
                                                                          multiple = TRUE,
                                                                          options = list(maxItems = 20,
                                                                                         plugins = list('remove_button')))  # enable deselection of items
                                                       )
                                                   ),
                                                   # burden plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.burden.plot'",
                                                       helper(radioGroupButtons("tmb.type",
                                                                                label = "Burden Type", choices = tmb.type.options, selected = tmb.type.options[1],
                                                                                status = "primary", individual = TRUE),
                                                              content = "tmb_type"),
                                                       helper(sliderInput("tmb.cutoff", "Burden Cutoff", min = 1, max = 10, value = 5, step = 1),
                                                              content = "tmb_cutoff")
                                                   ),
                                                   # waterfall plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.waterfall.plot'",
                                                       helper(numericInput("waterfall.cutoff", "Number of most mutated genes (1-30)", min = 1, max = 30, value = 20),
                                                              content = "waterfall_cutoff")
                                                   ),
                                                   # protein plot specific settings
                                                   conditionalPanel(
                                                       condition = "input.plotType == 'mut.protein.plot'",
                                                       selectizeInput("protein.plot.gene", "Gene", choices = names(mutated.gene.columns), multiple=FALSE)
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   helper(h4("Mutation Selection"),
                                                          content = "mutation_selection"),

                                                   radioGroupButtons("mutationSelection",
                                                                     label = "", choices = mutation.selection.options, selected = "mutations.all",
                                                                     status = "primary", individual = TRUE),
                                                   pickerInput(inputId = "mutationEffect",
                                                               label = "Mutation Effect",
                                                               choices = mut.effect.tbl$effect,
                                                               choicesOpt = list(content = mut.effect.html.labels),
                                                               selected = mut.effect.tbl$effect,
                                                               options = list(
                                                                 `actions-box` = TRUE
                                                               ),
                                                               multiple = TRUE
                                                   )
                                         ),
                                         wellPanel(style = wellpanel.settings.style,
                                                   helper(h4("Sample Selection"),
                                                          content = "sample_selection"),

                                                   selectInput("treatment.input", "Treatment Group",
                                                               choices = treatment.options),

                                                   radioGroupButtons("hist.type",
                                                                     label = "Histological Type", choices = ui.options$HistType, selected = ui.options$HistType[["Any"]],
                                                                     status = "primary", individual = TRUE),

                                                   radioGroupButtons("hr.cutoff",
                                                                     label = "ER/PgR Cutoff", choices = hr.cutoff.options, selected = hr.cutoff.options[1],
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
                                         conditionalPanel(
                                             condition = "input.plotType == 'mut.gene.plot' ||
                                                          input.plotType == 'mut.pathway.plot' ||
                                                          input.plotType == 'mut.burden.plot'",
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Plot Dimensions"),
                                                              content = "plot_dimensions"),
                                                       splitLayout(
                                                           numericInput("height.survival",
                                                                        label = "Height",
                                                                        value = 700
                                                           ),
                                                           numericInput("width.survival",
                                                                        label = "Width",
                                                                        value = 700
                                                           )
                                                       )
                                             ),
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Color Settings"),
                                                              content = "color_settings"),
                                                       selectInput("color.palette", "Color Palette",
                                                                   choices = color.palette.options, selected = color.palette.options[["JCO"]],
                                                                   multiple = FALSE)
                                             ),
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Legend Settings"),
                                                              content = "legend_settings"),
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
                                                       helper(h4("Plot Features"),
                                                              content = "plot_features"),
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
                                         ),
                                         conditionalPanel(
                                             condition = "input.plotType == 'mut.waterfall.plot'",
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Plot Dimensions"),
                                                              content = "plot_dimensions"),
                                                       splitLayout(
                                                           numericInput("height.waterfall",
                                                                        label = "Height",
                                                                        value = 700
                                                           ),
                                                           numericInput("width.waterfall",
                                                                        label = "Width",
                                                                        value = 1200
                                                           )
                                                       )
                                             ),
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("General Settings"),
                                                              content = "plot_settings_waterfall"),
                                                       numericInput("waterfall.round.digits",
                                                                    label = "Frequency Digits (0-3)",
                                                                    value = 0,
                                                                    min = 0,
                                                                    max = 3
                                                       ),
                                                       awesomeCheckbox("includeWildtypeSamples",
                                                                       label = "Include Wildtype Samples",
                                                                       value = TRUE)
                                             )
                                         ),
                                         conditionalPanel(
                                             condition = "input.plotType == 'mut.protein.plot'",
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Plot Dimensions"),
                                                              content = "plot_dimensions"),
                                                       splitLayout(
                                                           numericInput("height.protein",
                                                                        label = "Height",
                                                                        value = 300
                                                           ),
                                                           numericInput("width.protein",
                                                                        label = "Width",
                                                                        value = 1000
                                                           )
                                                       )
                                             ),
                                             wellPanel(style = wellpanel.settings.style,
                                                       helper(h4("Plot Settings"),
                                                              content = "protein_plot_settings"),
                                                       numericInput("protein.plot.mutation.cutoff",
                                                                    label = "Mutation Cutoff",
                                                                    value = 0
                                                       ),
                                                       numericInput("protein.plot.anno.cutoff",
                                                                    label = "Annotation Cutoff",
                                                                    value = 10
                                                       )
                                             )
                                         )
                                )
                    ),
                    wellPanel(style = wellpanel.settings.style,
                              helper(h3("Downloads"),
                                     content = "downloads"),
                           splitLayout(
                               downloadButton("downloadPlot", label = "Plot PDF"),
                               downloadButton("downloadSamples", label = "Sample TSV"),
                               downloadButton("downloadMutations", label = "Mutation TSV")
                           ),
                           HTML("<br>")  # for vertical space below the buttons
                    )
                ),

                # Plots, data tables, and statistics
                mainPanel(
                    tabsetPanel(type = "pills",
                                tabPanel("Plots", value = "plotTab", style = "padding-top:30px",
                                         withSpinner(plotOutput("plot"))
                                ),
                                tabPanel("Sample Table", value = "sampleTab",
                                         DT::dataTableOutput("sample.table")
                                ),
                                tabPanel("Mutation Table", value = "mutationTab",
                                         DT::dataTableOutput("mut.table")
                                ),
                                tabPanel("Dataset Statistics", value = "statsTab",
                                         htmlOutput("datasetStats")
                                ),
                                tabPanel("Citation and About", value = "aboutTab",
                                         htmlOutput("appCiteAbout"),
                                         h2("Session Information"),
                                         verbatimTextOutput("sessionInfo")
                                )
                    )
                )
            )
        ))
))
