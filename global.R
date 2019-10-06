############################################################################
#
# Global state, loaded before ui.R and server.R
#
############################################################################

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyhelper))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sessioninfo))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(reactome.db))
source("R/resources.R")
source("R/plot_survival_ggplot.R")
source("R/plot_waterfall.R")
source("R/plot_protein.R")


config = read_yaml("config.yaml")

con <- DBI::dbConnect(RSQLite::SQLite(), config$db_file)
samples <- collect(tbl(con, "samples"))
mutations <- collect(tbl(con, "mutations"))
DBI::dbDisconnect(con)

# Gene<->Protein map, only the first mapping for each gene is retained
gene_protein_mapping <- read.csv(config$gene_protein_map_file, sep='\t', header = F, stringsAsFactors = F)
gene_protein_mapping <- dplyr::select(gene_protein_mapping, Protein = 1, Gene = 2)
gene_protein_mapping <- dplyr::group_by(gene_protein_mapping, Gene)
gene_protein_mapping <- dplyr::filter(gene_protein_mapping, row_number() == 1)

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