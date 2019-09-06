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


# Generate header to be written to downloaded flat files.
get.download.table.head <- function(input, header) {
    table_header = paste0(sprintf("## Download timestamp: %s\n", date()), header)

    # Write sample filter settings
    sample_filters = get.sample.filter.descriptions(input)
    table_header = paste0(table_header, "## Sample filters:\n")
    for (idx in seq_along(sample_filters)) {
        table_header = paste0(table_header, "## ", unlist(sample_filters[idx]), "\n")
    }
    table_header = paste0(table_header, "##\n")

    mutation_filters = get.mutation.filter.descriptions(input)
    table_header = paste0(table_header, "## Mutation filters:\n")
    for (idx in seq_along(mutation_filters)) {
        table_header = paste0(table_header, "## ", unlist(mutation_filters[idx]), "\n")
    }
    table_header = paste0(table_header, "##\n")

    return(table_header)
}


# Returns sample filter settings in text form.
get.sample.filter.descriptions <- function(input) {
    hr.cutoff.string = paste0("(", ifelse(input$hr.cutoff == "hr.1perc", "1", "10"), "% cutoff)")

    filter_descriptions = list(
        paste("Treatment:", input$treatment.input),
        paste("Histological Subtype:", selection.to.label.list$HistType[[as.character(input$hist.type)]]),
        paste("ER Status:", selection.to.label.list$ER[[as.character(input$er.status)]], hr.cutoff.string),
        paste("PgR Status:", selection.to.label.list$PgR[[as.character(input$pgr.status)]], hr.cutoff.string),
        paste("HER2 Status:", selection.to.label.list$HER2[[as.character(input$her2.status)]]),
        paste("Ki67 Status:", selection.to.label.list$Ki67[[as.character(input$ki67.status)]]),
        paste("Nottingham Histologial Grade (NHG):", selection.to.label.list$NHG[[as.character(input$nhg)]]),
        paste("PAM50 Subtype:", selection.to.label.list$PAM50[[as.character(input$pam50)]])
    )

    return(filter_descriptions)
}

# Returns mutation filter settings in text form.
get.mutation.filter.descriptions <- function(input) {
    filter_descriptions = list(
        paste("Dataset:", names(mutation.selection.options)[mutation.selection.options == input$mutationSelection]),
        paste("Types:", paste(input$mutationEffect, collapse = ", "))
    )
    return(filter_descriptions)
}

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
        if (input$hr.cutoff == "hr.1perc") {
            filters[["ER"]] = "ER_1perc == selection.to.label.list$ER[[as.character(input$er.status)]]"
        } else {
            filters[["ER"]] = "ER_10perc == selection.to.label.list$ER[[as.character(input$er.status)]]"
        }
    }
    if (input$pgr.status != ui.options$PgR[["Any"]]) {
        if (input$hr.cutoff == "hr.1perc") {
            filters[["PgR"]] = "PgR_1perc == selection.to.label.list$PgR[[as.character(input$pgr.status)]]"
        } else {
            filters[["PgR"]] = "PgR_10perc == selection.to.label.list$PgR[[as.character(input$pgr.status)]]"
        }
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
add.gene.mut.status <- function(input, sample.tbl, mut.tbl, gene.column.map, gene.status.map) {
    for (gene in input$gene.input) {
        mut.var = gene.column.map[[gene]]
        gene.status = gene.status.map[[gene]]
        if (!(mut.var %in% colnames(sample.tbl))) {
            sample.tbl <- mutate(sample.tbl, !!mut.var := as.factor(ifelse(SAMPLE %in% mut.tbl$SAMPLE[mut.tbl$gene.symbol == gene], gene.status[["abnormal"]], gene.status[["normal"]])))
        }
    }
    return(sample.tbl)
}

# Set mutation count and count per expressed MB depending on mutation set selection.
set.mutation.counts <- function(input, sample.tbl) {
    if (input$mutationSelection == "mutations.cosmic") {
        sample.tbl$current_mutation_count = sample.tbl$COSMIC_Mutation_Count
        sample.tbl$current_mutation_nonsynon_count = sample.tbl$COSMIC_Mutation_Nonsynon_Count
        sample.tbl$current_mutation_count_per_expressed_mb = sample.tbl$COSMIC_Mutation_Count_per_expressed_MB
        sample.tbl$current_mutation_nonsynon_count_per_expressed_mb = sample.tbl$COSMIC_Mutation_Nonsynon_Count_per_expressed_MB
    } else {
        sample.tbl$current_mutation_count = sample.tbl$Mutation_Count
        sample.tbl$current_mutation_nonsynon_count = sample.tbl$Mutation_Nonsynon_Count
        sample.tbl$current_mutation_count_per_expressed_mb = sample.tbl$Mutation_Count_per_expressed_MB
        sample.tbl$current_mutation_nonsynon_count_per_expressed_mb = sample.tbl$Mutation_Nonsynon_Count_per_expressed_MB
    }
    return(sample.tbl)
}

# Add mutational burden given a cutoff.
add.mutation.burden <- function(input, sample.tbl) {
    if (input$tmb.type == "tmb.absolute") {
        sample.tbl <- mutate(sample.tbl,
                             tumor_mutational_burden = as.factor(ifelse(current_mutation_nonsynon_count > input$tmb.cutoff, "High", "Low")))
    } else if (input$tmb.type == "tmb.normalized") {
        sample.tbl <- mutate(sample.tbl,
                             tumor_mutational_burden = as.factor(ifelse(current_mutation_nonsynon_count_per_expressed_mb > input$tmb.cutoff, "High", "Low")))
    } else {
        # shouldn't happen
    }
    return(sample.tbl)
}

# Determine pathway mutation status for each sample.
add.pathway.mut.status <- function(input, sample.tbl, mut.tbl) {
    if (input$pathwayType == "pathway.reactome") {
        for (pathway in input$pathway.input) {
            mut.var = paste0("mut.pathway.status.", pathway)
            if (!(mut.var %in% colnames(sample.tbl))) {
                sample.tbl <- mutate(sample.tbl, !!mut.var := as.factor(ifelse(SAMPLE %in% mut.tbl$SAMPLE[grepl(pathway, mut.tbl$Pathways.Reactome)], "mut", "wt")))
            }
        }
    } else if (input$pathwayType == "pathway.custom") {
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
    } else {
        # should not happen
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
    } else if (input$plotType == "mut.protein.plot") {
        mut.tbl <- filter(mut.tbl, gene.symbol %in% input$protein.plot.gene, TYPE == "SNV")
    } else if (input$plotType == "mut.pathway.plot") {
        if (input$pathwayType == "pathway.reactome") {
            # keep mutations present in any of the input pathways
            mut.tbl <- filter(mut.tbl, apply(sapply(input$pathway.input, function(pathway) grepl(pathway, Pathways.Reactome)), 1, any))
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
    } else if (input$plotType == "mut.waterfall.plot") {
        # count the occurrence of each mutation in our set
        mut_count <- plyr::count(plyr::count(mut.tbl, c('gene.symbol', 'SAMPLE'))[, 1:2], 'gene.symbol')

        # determine the top X most mutated genes
        topX.mut <- mut_count %>%
            arrange(desc(freq)) %>%
            head(input$waterfall.cutoff)
        topX.mut <- as.character(topX.mut$gene.symbol)

        # Restrict the mutation table to the genes we'll actually display.
        mut.tbl = dplyr::filter(mut.tbl, gene.symbol %in% topX.mut)
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

# Determine suitable row and column counts for grid plotting.
get.plot.grid.dimensions <- function(n.plots) {
    n.cols = n.rows = 1
    while (n.cols * n.rows < n.plots) {
        if (n.cols == n.rows) n.cols = n.cols + 1
        else n.rows = n.rows + 1
    }
    return(list(rows=n.rows, cols=n.cols))
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
        dplyr::select(Protein = 1, Gene = 2) %>%
        group_by(Gene) %>%
        filter(row_number() == 1)

    # Generate a mapping from gene to mutation status column names and add the custom columns.
    # Two types of "genes" are considered:
    # 1. Normal genes as defined by UCSC, Ensembl, etc. These have status mutated (mut), or wildtype (wt).
    # 2. Custom genes, with details and status defined in config.yaml. If no status is defined, "mut", and "wt"
    #    are used by default.
    # Currently only two statuses are supported, deemed "abnormal" and "normal".
    mutated.genes <- sort(unique(mutations$gene.symbol))
    mutated.gene.columns <- paste0("mut.status.", mutated.genes)
    names(mutated.gene.columns) = mutated.genes
    mutated.gene.status <- rep(list(list(abnormal="mut", normal="wt")), length(mutated.genes))
    names(mutated.gene.status) = mutated.genes
    for (val in config$custom_genes) {
        mutated.gene.columns[val$label] = val$column
        if (is.null(val$status)) {
            this.status = c("mut", "wt")
        } else {
            this.status = as.list(val$status)
        }
        names(this.status) = c("abnormal", "normal")
        mutated.gene.status[[val$column]] = this.status
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
        filtered.samples <- add.gene.mut.status(input, filtered.samples, mut.tbl(), mutated.gene.columns, mutated.gene.status)
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
    observeEvent(c(input$mutationSelection, input$tmb.type), {
        updateSliderInput(session, "tmb.cutoff",
                          min = round(min(switch(input$tmb.type,
                                                 tmb.absolute = sample.tbl()[["current_mutation_nonsynon_count"]],
                                                 tmb.normalized = sample.tbl()[["current_mutation_nonsynon_count_per_expressed_mb"]])),
                                      digits = 2),
                          max = round(max(switch(input$tmb.type,
                                                 tmb.absolute = sample.tbl()[["current_mutation_nonsynon_count"]],
                                                 tmb.normalized = sample.tbl()[["current_mutation_nonsynon_count_per_expressed_mb"]])),
                                      digits = 2),
                          value = round(median(switch(input$tmb.type,
                                                      tmb.absolute = sample.tbl()[["current_mutation_nonsynon_count"]],
                                                      tmb.normalized = sample.tbl()[["current_mutation_nonsynon_count_per_expressed_mb"]])),
                                        digits = 2),
                          step = ifelse(input$tmb.type == "tmb.absolute", 1, 0.01))
        updateNumericInput(session, "waterfall.cutoff",
                           max = length(mutated.genes))
    })

    # Hide the loading message when the rest of the server function has executed
    hideElement(id = "loading-content", anim = TRUE, animType = "fade")
    showElement(id = "app-content")

    # Plot using ggplot2
    output$plot <- renderPlot(
        height = function(x) plot.height(),
        width = function(x) plot.width(),
        {
            # Save plot object in case we need it for PDF download later.
            current.plot <<- create.plot()

            if (input$plotType == "mut.waterfall.plot") {
                grid::grid.draw(current.plot)
            } else {
                print(current.plot)
            }
        })

    create.plot <- function() {
        sample.data <- sample.tbl()
        treatment.label = treatment.tbl$plot.label[which(treatment.tbl$var == input$treatment.input)]

        if (input$plotType == "mut.burden.plot") {
            fit = survfit(Surv(OS_years, OS_event) ~ tumor_mutational_burden, data = sample.data)
            title = paste0("Mutation Burden (Cutoff ", input$tmb.cutoff, ") in ", treatment.label, " Treated Patients")
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else if (input$plotType == "mut.pathway.plot" & input$pathwayType == "pathway.reactome") {
            plot.list = list()

            # Determine plot grid dimensions
            grid.dims = get.plot.grid.dimensions(length(input$pathway.input))
            n.cols = grid.dims[["cols"]]
            n.rows = grid.dims[["rows"]]

            # scale plot dimensions to new settings
            # XXX currently resets user-specified dimensions
            updateNumericInput(session, "height.survival", value = round(500 + ((n.rows + log(n.rows)) * 100)))
            updateNumericInput(session, "width.survival", value = round(500 + ((n.cols + log(n.cols)) * 100)))

            for (pathway in input$pathway.input) {
                mut.var = paste0("mut.pathway.status.", pathway)

                # Call survfit with do.call to avoid a problem with ggsurvplot later on.
                # See: https://github.com/kassambara/survminer/issues/125
                fit <- do.call(survfit,
                               list(formula = Surv(OS_years, OS_event) ~ get(mut.var), data = sample.data))

                plot.list[[pathway]] = surv.plot(input, fit, data=sample.data, gene=pathway, title=pathway)
            }
            title.main = paste("Treatment Group: ", treatment.label)
            title.grob = text_grob(title.main, size = 23, face = "bold")
            plot = arrange_ggsurvplots(plot.list, nrow=n.rows, ncol=n.cols, byrow=TRUE, title=title.grob)
        } else if (input$plotType == "mut.pathway.plot" & input$pathwayType == "pathway.custom") {
            fit = survfit(Surv(OS_years, OS_event) ~ mut.pathway.status, data = sample.data)
            title = paste("Treatment Group:", treatment.label)
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else if (input$plotType == "mut.gene.plot") {
            plot.list = list()

            # Determine plot grid dimensions
            grid.dims = get.plot.grid.dimensions(length(input$gene.input))
            n.cols = grid.dims[["cols"]]
            n.rows = grid.dims[["rows"]]

            # scale plot dimensions to new settings
            # XXX currently resets user-specified dimensions
            updateNumericInput(session, "height.survival", value = round(500 + ((n.rows + log(n.rows)) * 100)))
            updateNumericInput(session, "width.survival", value = round(500 + ((n.cols + log(n.cols)) * 100)))

            for (gene in input$gene.input) {
                mut.var = mutated.gene.columns[[gene]]

                # Call survfit with do.call to avoid a problem with ggsurvplot later on.
                # See: https://github.com/kassambara/survminer/issues/125
                fit <- do.call(survfit,
                               list(formula = Surv(OS_years, OS_event) ~ get(mut.var), data = sample.data))

                plot.list[[gene]] = surv.plot(input, fit, data=sample.data, gene=gene, title=gene)
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
        display.mut.tbl = dplyr::select(mut.tbl(), -c(Pathways.Reactome))
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
        filename = function() { paste("mutation_explorer_plot", "pdf", sep='.') },
        content = function(file) {
            # if no plot, length(current.plot) == 0
            pdf(file, useDingbats = FALSE, width = plot.width() / 72, height = plot.height() / 72)
            if (input$plotType == "mut.waterfall.plot") {
                grid::grid.draw(current.plot)
            } else {
                print(current.plot)
            }
            dev.off()
        },
        contentType = "application/pdf"
    )

    output$downloadSamples <- downloadHandler(
        filename = function() { paste("samples", "zip", sep='.') },
        content = function(file) {
            tmpfile = gsub("(.+)\\..+", "\\1\\.tsv", file)
            cat(get.download.table.head(input, config$table_header), file=tmpfile)
            write.table(sample.tbl(), tmpfile, sep = "\t", quote = FALSE, na = "", row.names = FALSE, append = TRUE)
            zip(zipfile = file, files = tmpfile, flags = "-r9Xj")  # r9X is default; -j to trim input file names
        },
        contentType = "application/zip"
    )

    output$downloadMutations <- downloadHandler(
        filename = function() { paste("mutations", "zip", sep='.') },
        content = function(file) {
            tmpfile = gsub("(.+)\\..+", "\\1\\.tsv", file)
            cat(get.download.table.head(input, config$table_header), file=tmpfile)
            write.table(mut.tbl(), tmpfile, sep = "\t", quote = FALSE, na = "", row.names = FALSE, append = TRUE)
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
