#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(reactome.db))
source("R/resources.R")
source("R/plot_survival_ggplot.R")
source("R/plot_waterfall.R")
source("R/plot_protein.R")


# Return a list of samples adhering to the specified filters.
filter.sample.tbl <- function(input, sample.tbl) {
    filters = list()

    # Filter by treatment
    if (input$treatment.input != "any") {
        filters[["Treatment"]] = "get(input$treatment.input) == 1"
    }
    # Filter down by biomarker selection
    if (input$hist.type != ui.options$HistType[["Any"]]) {
        filters[["HistType"]] = "Histological_Type == selection.to.label.list$HistType[[as.character(input$hist.type)]]"
    }
    if (input$er.status != ui.options$ER[["Any"]]) {
        filters[["ER"]] = "ER_1perc == selection.to.label.list$ER[[as.character(input$er.status)]]"
    }
    if (input$pgr.status != ui.options$PgR[["Any"]]) {
        filters[["PgR"]] = "PgR_1perc == selection.to.label.list$PgR[[as.character(input$pgr.status)]]"
    }
    if (input$her2.status != ui.options$HER2[["Any"]]) {
        filters[["HER2"]] = "HER2 == selection.to.label.list$HER2[[as.character(input$her2.status)]]"
    }
    if (input$ki67.status != ui.options$Ki67[["Any"]]) {
        filters[["Ki67"]] = "Ki67 == selection.to.label.list$Ki67[[as.character(input$ki67.status)]]"
    }
    if (input$nhg != ui.options$NHG[["Any"]]) {
        filters[["NHG"]] = "NHG == selection.to.label.list$NHG[[as.character(input$nhg)]]"
    }
    if (input$pam50 != ui.options$PAM50[["Any"]]) {
        filters[["PAM50"]] = "PAM50 == selection.to.label.list$PAM50[[as.character(input$pam50)]]"
    }

    # apply filters
    if (length(filters) > 0) {
        filter_str = paste(filters, collapse = " & ")
        sample.tbl <- filter_(sample.tbl, filter_str)
    }

    samples <- as.character(sample.tbl$SAMPLE)

    return(samples)
}

# Determine mutation status of samples for selected genes
add.gene.mut.status <- function(input, sample.tbl, mut.tbl, gene.column.map) {
    for (gene in input$gene.input) {
        mut.var = gene.column.map[[gene]]
        if (!(mut.var %in% colnames(sample.tbl))) {
            sample.tbl <- mutate(sample.tbl, !!mut.var := as.factor(ifelse(SAMPLE %in% mut.tbl$SAMPLE[mut.tbl$gene.symbol == gene], "mut", "wt")))
        }
    }
    return(sample.tbl)
}

# Set mutation count depending on mutation set selection.
set.mutation.counts <- function(input, sample.tbl) {
    if (input$mutationSelection == "mutations.cosmic") {
        sample.tbl$current_mutation_count = sample.tbl$COSMIC_Mutation_Count
        sample.tbl$current_mutation_nonsynon_count = sample.tbl$COSMIC_Mutation_Nonsynon_Count
    } else {
        sample.tbl$current_mutation_count = sample.tbl$Mutation_Count
        sample.tbl$current_mutation_nonsynon_count = sample.tbl$Mutation_Nonsynon_Count
    }
    return(sample.tbl)
}

# Add mutational burden given a cutoff.
add.mutation.burden <- function(input, sample.tbl) {
    sample.tbl <- mutate(sample.tbl, tumor_mutational_burden = as.factor(ifelse(current_mutation_count > input$tmb.cutoff, "High", "Low")))
    return(sample.tbl)
}

# Determine pathway mutation status for each sample.
add.pathway.mut.status <- function(input, sample.tbl, mut.tbl) {
    if (input$pathwayType == "pathway.reactome") {
        sample.tbl <- mutate(sample.tbl, mut.pathway.status = as.factor(ifelse(sample.tbl$SAMPLE %in% mut.tbl$SAMPLE, "mut", "wt")))
    } else {  # pathway.custom
        pathway.genes <- input$custom.pathway.input

        if (!is.null(pathway.genes)) {
            mut.status = apply(sample.tbl, 1, function(sample) {
                pathway.genes %in% mut.tbl$gene.symbol[mut.tbl$SAMPLE %in% sample[["SAMPLE"]]]
                })
            if (is.vector(mut.status)) {  # one gene only, logical vector
                sample.tbl$mut.pathway.status = as.factor(ifelse(mut.status, "mut", "wt"))
            } else {
                sample.tbl$mut.pathway.status = as.factor(ifelse(apply(mut.status, 2, any), "mut", "wt"))
            }
        } else {
            sample.tbl$mut.pathway.status = as.factor(logical(nrow(sample.tbl)))  # FALSE vector
        }
    }
    return(sample.tbl)
}

filter.mut.tbl <- function(input, sample.list, mut.tbl, gene.column.map) {
    mut.tbl <- filter(mut.tbl, SAMPLE %in% sample.list)

    #
    # Filter mutations based on input selections.
    #
    if (input$mutationSelection == "mutations.cosmic") {
        mut.tbl <- filter(mut.tbl, COSMIC_ID != ".")
    }
    if (!is.null(input$mutationEffect)) {
        mut.tbl <- filter(mut.tbl, ANN.effect.class %in% input$mutationEffect)
    }

    #
    # Filter mutations based on plot type.
    #
    if (input$plotType == "mut.gene.plot") {
        mut.tbl <- filter(mut.tbl, gene.symbol %in% input$gene.input)
    } else if (input$plotType == "mut.pathway.plot") {
        if (input$pathwayType == "pathway.reactome") {
            mut.tbl <- filter(mut.tbl, grepl(input$pathway.input, Pathways.Reactome))
        } else {  # pathway.custom
            pathway.genes <- input$custom.pathway.input

            if (!is.null(pathway.genes)) {
                # keep all genes in the selected custom pathway
                mut.tbl = filter(mut.tbl, gene.symbol %in% pathway.genes)
            } else {
                mut.tbl = filter(mut.tbl, FALSE)  # no genes defined -> empty table
            }
        }
    } else if (input$plotType == "mut.burden.plot") {
        mut.tbl <- filter(mut.tbl, ANN.effect.class != "Synonymous")
    }

    return(mut.tbl)
}

# Returns a table with various sample/mutation descriptive statistics.
get.dataset.stats <- function(sample.tbl, mut.tbl) {
    stats = data.frame(Value=character(), Stat=numeric()) %>%
        add_row(Value="Total Samples", Stat=nrow(sample.tbl)) %>%
        add_row(Value="Total Mutations", Stat=nrow(mut.tbl)) %>%
        add_row(Value="Total COSMIC Mutations", Stat=nrow(filter(mut.tbl, COSMIC_ID != "."))) %>%
        add_row(Value="Mean Overall Mutations per Sample", Stat=mean(sample.tbl$current_mutation_count)) %>%
        add_row(Value="Median Overall Mutations per Sample", Stat=median(sample.tbl$current_mutation_count)) %>%
        add_row(Value="Mean Coding Mutations per Sample", Stat=mean(sample.tbl$current_mutation_nonsynon_count)) %>%
        add_row(Value="Median Coding Mutations per Sample", Stat=median(sample.tbl$current_mutation_nonsynon_count)) %>%
        add_row(Value="Median Overall Survival (in Months)", Stat=median(sample.tbl$OS_months))
    return(stats)
}


############################################################################
#
# Server logic
#
############################################################################S

shinyServer(function(input, output, session) {

    config = read_yaml("config.yaml")

    con <- DBI::dbConnect(RSQLite::SQLite(), config$db_file)
    samples <- tbl(con, "samples") %>%
        collect()
    mutations <- tbl(con, "mutations") %>%
        collect()
    dbDisconnect(con)

    # Gene<->Protein map, only the first mapping for each gene is retained
    gene_protein_mapping <- read.csv(config$gene_protein_map_file, sep='\t', header = F, stringsAsFactors = F) %>%
        select(Protein = 1, Gene = 2) %>%
        group_by(Gene) %>%
        filter(row_number() == 1)

    # Generate a mapping from gene to mutation status column names and add the custom columns
    mutated.genes <- sort(unique(mutations$gene.symbol))
    mutated.gene.columns <- paste0("mut.status.", mutated.genes)
    names(mutated.gene.columns) = mutated.genes
    for (val in config$custom_genes) {
        mutated.gene.columns[val$label] = val$column
    }

    n.mut = nrow(mutations)
    n.samples = nrow(samples)

    # set default directory for help files
    observe_helpers(session, "helpfiles")

    sample.list <- reactive({
        samples <- filter.sample.tbl(input, samples)
        return(samples)
    })

    # Add additional information to the current sample set as specified in the input controls.
    sample.tbl <- reactive({
        filtered.samples <- filter(samples, SAMPLE %in% sample.list())
        filtered.samples <- add.gene.mut.status(input, filtered.samples, mut.tbl(), mutated.gene.columns)
        filtered.samples <- set.mutation.counts(input, filtered.samples)
        filtered.samples <- add.mutation.burden(input, filtered.samples)
        filtered.samples <- add.pathway.mut.status(input, filtered.samples, mut.tbl())
        return(filtered.samples)
    })

    mut.tbl <- reactive({
        filtered.muts <- filter.mut.tbl(input, sample.list(), mutations, mutated.gene.columns)
        return(filtered.muts)
    })

    plot.height = reactive({
        if (input$plotType == "mut.survival.plot") {
            height = input$height.survival
        } else if (input$plotType == "mut.waterfall.plot") {
            height = input$height.waterfall
        } else if (input$plotType == "mut.protein.plot") {
            height = input$height.protein
        } else {
            height = 700
        }
        return(height)
    })

    plot.width = reactive({
        if (input$plotType == "mut.survival.plot") {
            width = input$width.survival
        } else if (input$plotType == "mut.waterfall.plot") {
            width = input$width.waterfall
        } else if (input$plotType == "mut.protein.plot") {
            width = input$width.protein
        } else {
            width = 700
        }
        return(width)
    })

    ######################################################
    #
    # Update UI elements
    #
    ######################################################
    output$header_panel <- renderUI({
        titlePanel(
            h1(config$app_title,
               h3(paste(prettyNum(n.mut, big.mark=","),
                        "mutations in",
                        prettyNum(n.samples, big.mark=","),
                        "Primary Breast Cancer Samples")))
        )
    })
    outputOptions(output, "header_panel", suspendWhenHidden=FALSE)  # make sure the header is shown before the loading screen is gone
    updateSelectInput(session, "gene.input",
                      choices = names(mutated.gene.columns)
    )
    updateSelectInput(session, "protein.plot.gene",
                      choices = names(mutated.gene.columns)
    )
    updateSelectInput(session, "custom.pathway.input",
                      choices = names(mutated.gene.columns)
    )
    observeEvent(input$mutationSelection, {
        updateSliderInput(session, "tmb.cutoff",
                          min = min(sample.tbl()[["current_mutation_count"]]),
                          max = max(sample.tbl()[["current_mutation_count"]]),
                          value = median(sample.tbl()[["current_mutation_count"]]))
        updateNumericInput(session, "waterfall.cutoff",
                           max = length(mutated.genes))
    })

    # Hide the loading message when the rest of the server function has executed
    hideElement(id = "loading-content", anim = TRUE, animType = "fade")
    showElement(id = "app-content")

    # Survival plot using ggplot2
    output$survplot <- renderPlot(
        height = function(x) plot.height(),
        width = function(x) plot.width(),
        {
            current.plot <<- create.plot()
            return(current.plot)
        })

    create.plot <- function() {
        sample.data <- sample.tbl()
        treatment.label = treatment.tbl$plot.label[which(treatment.tbl$var == input$treatment.input)]

        if (input$plotType == "mut.burden.plot") {
            fit = survfit(Surv(OS_years, OS_event) ~ tumor_mutational_burden, data = sample.data)
            title = paste0("Mutation Burden (Cutoff ", input$tmb.cutoff, ") in ", treatment.label, " Treated Patients")
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else if (input$plotType == "mut.pathway.plot") {
            fit = survfit(Surv(OS_years, OS_event) ~ mut.pathway.status, data = sample.data)
            title = paste("Pathway in", treatment.label, "Treated Patients")
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else if (input$plotType == "mut.gene.plot") {
            plot.list = list()

            # brute-force determine the row/col counts
            n.cols = n.rows = 1
            while (n.cols * n.rows < length(input$gene.input)) {
                if (n.cols == n.rows) n.cols = n.cols + 1
                else n.rows = n.rows + 1
            }
            # scale plot dimensions to new settings
            # XXX currently resets user-specified dimensions
            updateNumericInput(session, "height.survival", value = round(500 + ((n.rows + log(n.rows)) * 100)))
            updateNumericInput(session, "width.survival", value = round(500 + ((n.cols + log(n.cols)) * 100)))

            for (gene in input$gene.input) {
                title.gene = paste(gene, "Mutation Status")

                mut.var = mutated.gene.columns[[gene]]

                # Call survfit with do.call to avoid a problem with ggsurvplot later on.
                # See: https://github.com/kassambara/survminer/issues/125
                fit <- do.call(survfit,
                               list(formula = Surv(OS_years, OS_event) ~ get(mut.var), data = sample.data))

                plot.list[[gene]] = surv.plot(input, fit, data=sample.data, gene=gene, title=title.gene)
            }
            title.main = paste("Treatment Group: ", treatment.label)
            title.grob = text_grob(title.main, size = 23, face = "bold")
            plot = arrange_ggsurvplots(plot.list, nrow=n.rows, ncol=n.cols, byrow=TRUE, title=title.grob)
        } else if (input$plotType == "mut.waterfall.plot") {
            plot = plot.waterfall(input, sample.data, mut.tbl(), mutated.gene.columns)
        } else if (input$plotType == "mut.protein.plot") {
            plot = plot.protein(input, mut.tbl(), gene_protein_mapping)
        } else {
            # should not happen
        }
        return(plot)
    }

    output$sample.table = DT::renderDataTable({
        DT::datatable(sample.tbl(),
                      caption = "Samples",
                      extensions = c('FixedColumns'),
                      selection = "none",
                      options = list(
                          orderClasses = TRUE,
                          searching = FALSE,
                          scrollX = TRUE,
                          fixedColumns = TRUE
                      ))
    })

    output$mut.table = DT::renderDataTable({
        display.mut.tbl = select(mut.tbl(), -c(Pathways.Reactome))
        DT::datatable(display.mut.tbl,
                      caption = "Mutations",
                      extensions = c('FixedColumns'),
                      selection = "none",
                      options = list(
                          orderClasses = TRUE,
                          searching = FALSE,
                          scrollX = TRUE,
                          fixedColumns = TRUE
                      ))
    })

    output$downloadPlot <- downloadHandler(
        filename = function() { paste("mutation_explorer_survival_plot", "pdf", sep='.') },
        content = function(file) {
            # if no plot, length(current.plot) == 0
            pdf(file, useDingbats = FALSE, width = plot.width() / 72, height = plot.height() / 72)
            print(current.plot, newpage = FALSE)
            dev.off()
        },
        contentType = "application/pdf"
    )

    output$downloadSamples <- downloadHandler(
        filename = function() { paste("samples", "zip", sep='.') },
        content = function(file) {
            tmpfile = gsub("(.+)\\..+", "\\1\\.tsv", file)
            write.table(sample.tbl(), tmpfile, sep = "\t", quote = FALSE, na = "", row.names = FALSE)
            zip(zipfile = file, files = tmpfile, flags = "-r9Xj")  # r9X is default; -j to trim input file names
        },
        contentType = "application/zip"
    )

    output$downloadMutations <- downloadHandler(
        filename = function() { paste("mutations", "zip", sep='.') },
        content = function(file) {
            tmpfile = gsub("(.+)\\..+", "\\1\\.tsv", file)
            write.table(mut.tbl(), tmpfile, sep = "\t", quote = FALSE, na = "", row.names = FALSE)
            zip(zipfile = file, files = tmpfile, flags = "-r9Xj")  # r9X is default; -j to trim input file names
        },
        contentType = "application/zip"
    )

    output$datasetStats <- renderUI ({
        div(
            h2("Statistics for Total Sample/Mutation Set"),
            renderTable(get.dataset.stats(set.mutation.counts(input, samples), mutations), colnames = FALSE),
            h2("Statistics for Selected Sample/Mutation Set"),
            renderTable(get.dataset.stats(sample.tbl(), mut.tbl()), colnames = FALSE)
        )
    })

    output$appCiteAbout <- renderUI ({
        includeMarkdown("about.md")
    })

})
