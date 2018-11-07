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
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(reactome.db))
source("R/resources.R")
source("R/plot_survival_ggplot.R")


filter.sample.tbl <- function(input, master) {
    filters = list()

    # Filter by treatment
    if (input$treatment.input != "any") {
        filters[["Treatment"]] = "get(input$treatment.input) == 1"
    }
    # Filter down by biomarker selection
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
        master <- master %>% filter_(filter_str)
    }

    return(master)
}

# Determine mutation status of samples for selected genes
add.gene.mut.status <- function(input, sample.tbl, mut.tbl) {
    for (gene in input$gene.input) {
        mut.var = paste0("mut.status.", gene)
        sample.tbl <- sample.tbl %>%
            mutate(!!mut.var := as.factor(ifelse(rba %in% mut.tbl$SAMPLE[mut.tbl$gene.symbol == gene], "mut", "wt")))
    }
    return(sample.tbl)
}

# Add mutational burden given a cutoff.
add.mutation.burden <- function(input, sample.tbl) {
    sample.tbl <- sample.tbl %>%
        mutate(tumor_mutational_burden = as.factor(ifelse(mutation_count > input$tmb.cutoff, "High", "Low")))
    return(sample.tbl)
}

# Determine pathway mutation status for each sample.
add.pathway.mut.status <- function(input, sample.tbl, mut.tbl) {
    mut.tbl = filter.mut.tbl(input, sample.tbl, mut.tbl)

    if (input$pathwayType == "pathway.reactome") {
        sample.tbl <- sample.tbl %>%
            mutate(mut.pathway.status = as.factor(ifelse(sample.tbl$rba %in% mut.tbl$SAMPLE, "mut", "wt")))
    } else {  # pathway.custom
        pathway.genes <- input$custom.pathway.input

        if (!is.null(pathway.genes)) {
            mut.status = apply(sample.tbl, 1, function(sample) { pathway.genes %in% mut.tbl$gene.symbol[mut.tbl$SAMPLE %in% sample[["rba"]]] })
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

filter.mut.tbl <- function(input, sample.tbl, mut.tbl) {
    mut.tbl <- mut.tbl %>%
        filter(SAMPLE %in% sample.tbl$rba)

    if (input$plotType == "mut.gene.plot") {
        mut.tbl <- mut.tbl %>%
            filter(gene.symbol %in% input$gene.input)
    } else if (input$plotType == "mut.pathway.plot") {
        if (input$pathwayType == "pathway.reactome") {
            mut.tbl <- mut.tbl %>%
                filter(grepl(input$pathway.input, Pathways.Reactome))
        } else {  # pathway.custom
            pathway.genes <- input$custom.pathway.input

            if (!is.null(pathway.genes)) {
                # keep all genes in the selected custom pathway
                mut.tbl = mut.tbl %>% filter(gene.symbol %in% pathway.genes)
            } else {
                mut.tbl = mut.tbl %>% filter(FALSE)  # no genes defined -> empty table
            }
        }
    }  # No filtering for burden; we want all mutations in that case.

    return(mut.tbl)
}


############################################################################
#
# Server logic
#
############################################################################S

shinyServer(function(input, output, session) {

    master <- as.data.frame(get(load("data/sample_master_table_scanb.Rdata")))
    mutations <- as.data.frame(get(load("data/mutations_scanb.Rdata")))
    mutated.genes <- sort(unique(mutations$gene.symbol))
    n.mut = nrow(mutations)
    n.samples = nrow(master)

    # set default directory for help files
    observe_helpers(session, "helpfiles")

    # Filter the sample table down whenever an input control changes
    sample.tbl <- reactive({
        filtered.table <- filter.sample.tbl(input, master)
        filtered.table <- add.gene.mut.status(input, filtered.table, mutations)
        filtered.table <- add.mutation.burden(input, filtered.table)
        filtered.table <- add.pathway.mut.status(input, filtered.table, mutations)
        return(filtered.table)
    })

    ######################################################
    #
    # Update UI elements
    #
    ######################################################
    output$header_panel <- renderUI({
        n.samples = nrow(master)
        titlePanel(
            h1("SCAN-B Mutation Explorer",
               h3(paste(prettyNum(n.mut, big.mark=","), "mutations in", prettyNum(n.samples, big.mark=","), "Primary Breast Cancer Samples")))
        )
    })
    outputOptions(output, "header_panel", suspendWhenHidden=FALSE)  # make sure the header is shown before the loading screen is gone
    updateSelectInput(session, "gene.input",
                      choices = mutated.genes
    )
    updateSelectInput(session, "custom.pathway.input",
                      choices = mutated.genes
    )
    updateSliderInput(session, "tmb.cutoff", min = min(master$mutation_count),
                      max = max(master$mutation_count), value = 75)

    # Hide the loading message when the rest of the server function has executed
    hideElement(id = "loading-content", anim = TRUE, animType = "fade")
    showElement(id = "app-content")

    # Survival plot using ggplot2
    output$survplot <- renderPlot(
        height = function(x) input$height,
        width = function(x) input$width,
        {
            plot.surv <<- plot.survival()
            return(plot.surv)
        })

    plot.survival <- function() {
        sample.data <- sample.tbl()
        treatment.label = treatment.tbl$plot.label[which(treatment.tbl$var == input$treatment.input)]

        if (input$plotType == "mut.burden.plot") {
            fit = survfit(Surv(OS_years, OS_event) ~ tumor_mutational_burden, data = sample.data)
            title = paste0("Mutation Burden (Cutoff ", input$tmb.cutoff, ") in ", treatment.label, " treated patients")
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else if (input$plotType == "mut.pathway.plot") {
            fit = survfit(Surv(OS_years, OS_event) ~ mut.pathway.status, data = sample.data)
            title = paste("Pathway in", treatment.label, "treated patients")
            plot = surv.plot(input, fit, data=sample.data, title=title)
        } else {  # mut.gene.plot
            plot.list = list()
            title.main = paste(treatment.label, "treated patients")

            # brute-force determine the row/col counts
            n.cols = n.rows = 1
            while (n.cols * n.rows < length(input$gene.input)) {
                if (n.cols == n.rows) n.cols = n.cols + 1
                else n.rows = n.rows + 1
            }
            # scale plot dimensions to new settings
            # XXX currently resets user-specified dimensions
            updateNumericInput(session, "height", value = round(500 + ((n.rows + log(n.rows)) * 100)))
            updateNumericInput(session, "width", value = round(500 + ((n.cols + log(n.cols)) * 100)))

            for (gene in input$gene.input) {
                title.gene = paste(gene, "Mutation Status")

                mut.var = paste0("mut.status.", gene)

                # Call survfit with do.call to avoid a problem with ggsurvplot later on.
                # See: https://github.com/kassambara/survminer/issues/125
                fit <- do.call(survfit,
                               list(formula = Surv(OS_years, OS_event) ~ get(mut.var), data = sample.data))

                plot.list[[gene]] = surv.plot(input, fit, data=sample.data, gene=gene, title=title.gene)
            }
            plot = arrange_ggsurvplots(plot.list, nrow=n.rows, ncol=n.cols, byrow=TRUE, title=title.main)
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
        filtered.muts <- filter.mut.tbl(input, sample.tbl(), mutations)
        DT::datatable(filtered.muts,
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
            # if no plot, length(plot.surv) == 0
            pdf(file, useDingbats = FALSE, width = input$width / 72, height = input$height / 72)
            print(plot.surv, newpage = FALSE)
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
            filtered.muts <- filter.mut.tbl(input, sample.tbl(), mutations)
            tmpfile = gsub("(.+)\\..+", "\\1\\.tsv", file)
            write.table(filtered.muts, tmpfile, sep = "\t", quote = FALSE, na = "", row.names = FALSE)
            zip(zipfile = file, files = tmpfile, flags = "-r9Xj")  # r9X is default; -j to trim input file names
        },
        contentType = "application/zip"
    )

    output$appCiteAbout <- renderUI ({
        HTML("<h2>The Sweden Cancerome Analysis Network&mdash;Breast (SCAN-B)</h2>
<p>The Sweden Cancerome Analysis Network&mdash;Breast (SCAN-B) initiative (ClinicalTrials.gov identifier <a href='https://clinicaltrials.gov/ct2/show/NCT02306096'>NCT02306096</a>)
is a population-based, multicenter breast cancer study that started enrolling patients in 2010.  The study is described
in <a href='https://doi.org/10.1186/s13073-015-0131-9'>Saal <i>et al</i>, Genome Medicine (2015)</a>.
<h2>SCAN-B Mutation Explorer</h2>
<p>This software was developed as part of a PhD research project in the laboratory of Lao H. Saal, Translational Oncogenomics Unit, Department of Oncology and Pathology, Lund University, Sweden.</p>
<h3>Citation</h3>
<p>If you use any data or plots from this website in your publications, please cite the following paper:</p>
<p>Brueffer <i>et al</i>. Paper Name. Journal Name.</p>
<small>The source code for this software can be found on GitHub: <a href='https://github.com/cbrueffer/ShinyMutationExplorer'>https://github.com/cbrueffer/ShinyMutationExplorer</a></small>
")
    })

})
