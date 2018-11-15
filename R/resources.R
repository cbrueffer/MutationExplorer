suppressPackageStartupMessages(library(reactome.db))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))

#
# Treatment definitions
#
treatment.tbl <- as.tibble(rbind(c("any", "Any", "any"),
                                 c("adjv_endo", "Endocrine + any", "endo + any"),
                                 c("adjv_cyto", "Chemo + any", "chemo + any"),
                                 c("adjv_her2", "HER2 + any", "HER2 + any"),
                                 c("adjv_endo_only", "Endocrine only", "endo only"),
                                 c("adjv_cyto_only", "Chemo only", "chemo only"),
                                 c("adjv_her2_only", "HER2 only", "HER2 only"),
                                 c("adjv_endo_and_cyto", "Endo + Chemo + any", "endo + chemo + any"),
                                 c("adjv_endo_and_her2", "Endo + HER2 + any", "endo + HER2 + any"),
                                 c("adjv_cyto_and_her2", "Cyto + HER2 + any", "chemo + HER2 + any"),
                                 c("adjv_endo_and_cyto_only", "Endo + Chemo only", "endo + chemo only"),
                                 c("adjv_endo_and_her2_only", "Endo + HER2 only", "endo + HER2 only"),
                                 c("adjv_cyto_and_her2_only", "Cyto + HER2 only", "chemo + HER2 only"),
                                 c("adjv_endo_and_cyto_and_her2", "Endo + Chemo + HER2", "endo + chemo + HER2"),
                                 c("adjv_no_systemic", "No Systemic Treatment", "not systemically")
))
colnames(treatment.tbl) <- c("var", "label", "plot.label")

treatment.options = character()
for (i in 1:nrow(treatment.tbl)) {
    row = treatment.tbl[i, ]
    treatment.options[[row$label]] = row$var
}


#
# Biomarker definitions
#
biomarker.tbl <- as.tibble(rbind(c("ER", "ER_1perc", "POS", "Positive", "", 1),
                                 c("ER", "ER_1perc", "NEG", "Negative", "", 2),
                                 c("PgR", "PgR_1perc", "POS", "Positive", "", 1),
                                 c("PgR", "PgR_1perc", "NEG", "Negative", "", 2),
                                 c("HER2", "HER2", "POS", "Amplified", "red", 1),
                                 c("HER2", "HER2", "NEG", "Normal", "green", 2),
                                 c("Ki67", "Ki67", "POS", "High", "red", 1),
                                 c("Ki67", "Ki67", "NEG", "Low", "green", 2),
                                 c("NHG", "NHG", "G1", "G1", "green", 1),
                                 c("NHG", "NHG", "G2", "G2", "yellow", 2),
                                 c("NHG", "NHG", "G3", "G3", "red", 3),
                                 c("PAM50", "PAM50", "LumA", "Luminal A", "darkblue", 1),
                                 c("PAM50", "PAM50", "LumB", "Luminal B", "lightskyblue", 2),
                                 c("PAM50", "PAM50", "HER2", "HER2-enriched", "deeppink", 3),
                                 c("PAM50", "PAM50", "Basal", "Basal-like", "red", 4),
                                 c("PAM50", "PAM50", "Normal", "Normal-like", "green", 5),
                                 c("PAM50", "PAM50", "Unclassified", "Unclassified", "gray", 6)
))
colnames(biomarker.tbl) <- c("marker", "var", "level", "label", "color", "selection.value")


# get all Reactome Homo sapiens pathways
pathway.list = as.list(reactomePATHNAME2ID)[grepl("^Homo sapiens:", names(as.list(reactomePATHNAME2ID)))]

get.pathway.ui.options = function(show.pathway.tbl, pathways) {
    result.list = list()
    for (this.pathway.name in names(pathways)) {
        this.pathway.id.list <- pathways[[this.pathway.name]]

        # a few pathways have two IDs
        for (this.pathway.id in this.pathway.id.list) {
            this.pathway.label = paste0(gsub(".+: (.+)", "\\1", this.pathway.name), " (", this.pathway.id, ")")
            result.list[[this.pathway.label]] = this.pathway.id
        }
    }
    return(result.list)
}
pathway.ui.options <- get.pathway.ui.options(show.pathway.tbl, pathway.list)


# Generates a list of lists mapping each label to the corresponding
# UI selection values; one list per biomarker.
get.ui.options = function(biomarker.tbl) {
    result.list = list()
    for (this.marker in unique(biomarker.tbl$marker)) {
        marker.list = list()
        biomarker.subset = biomarker.tbl %>%
            filter(marker == this.marker)

        marker.list[["Any"]] <- 42  # every UI element has an "any" option
        for (this.level in biomarker.subset$level) {
            selection.value = biomarker.subset %>% filter(level == this.level) %>% dplyr::select(selection.value) %>% as.integer()
            label = biomarker.subset %>% filter(level == this.level) %>% dplyr::select(label) %>% as.character()
            marker.list[[label]] = selection.value
        }
        result.list[[this.marker]] = marker.list
    }
    return(result.list)
}
ui.options <- get.ui.options(biomarker.tbl)

# Generates a list of lists mapping each selection value to the corresponding
# biomarker levels; one list per biomarker.
get.selection.to.biomarker.levels = function(biomarker.tbl) {
    result.list = list()
    for (this.marker in unique(biomarker.tbl$marker)) {
        marker.list = list()
        biomarker.subset = biomarker.tbl %>%
            filter(marker == this.marker)

        marker.list[["42"]] <- "Any"  # every UI element has an "any" option
        for (this.level in biomarker.subset$level) {
            selection.value = biomarker.subset %>% filter(level == this.level) %>% dplyr::select(selection.value) %>% as.character()
            marker.list[[selection.value]] = this.level
        }
        result.list[[this.marker]] = marker.list
    }
    return(result.list)
}
selection.to.label.list <- get.selection.to.biomarker.levels(biomarker.tbl)


# Generates PAM50 HTML labels
get.pam50.html.labels <- function(biomarker.tbl) {
    biomarker.tbl.pam50 <- biomarker.tbl %>% filter(marker == "PAM50")
    result.list = list()
    result.list[["Any"]] = "<div style='background: white; color: black; font-size: 110%; padding-left: 5px;'>Any</div>"
    for (i in 1:nrow(biomarker.tbl.pam50)) {
        subtype = biomarker.tbl.pam50[i, ]$level
        label = biomarker.tbl.pam50[i, ]$label
        color = biomarker.tbl.pam50[i, ]$color
        result.list[[subtype]] = paste0("<div style='background: ", color, "; color: white; font-size: 110%; padding-left: 5px;'>", label, "</div>")
    }
    return(result.list)
}
pam50.html.labels = get.pam50.html.labels(biomarker.tbl)

plot.type.options = c("Mutated Genes" = "mut.gene.plot",
                      "Mutated Pathways" = "mut.pathway.plot",
                      "Tumor Mutational Burden" = "mut.burden.plot")

pathway.type.options = c("Reactome" = "pathway.reactome",
                         "Custom" = "pathway.custom")

plot.legend.loc.options = c("Top" = "top",
                            "Bottom" = "bottom",
                            "Left" = "left",
                            "Right" = "right",
                            "Custom" = "custom")

color.palette.options = c("NPG" = "npg",
                          "AAAS" = "aaas",
                          "NEJM" = "nejm",
                          "Lancet" = "lancet",
                          "JAMA" = "jama",
                          "JCO" = "jco",
                          "UCSCGB" = "ucscgb",
                          "D3" = "d3",
                          "LocusZoom" = "locuszoom",
                          "IGV" = "igv",
                          "UChicago" = "uchicago",
                          "Star Trek" = "startrek",
                          "Tron Legacy" = "tron",
                          "Futurama" = "futurama",
                          "Rick and Morty" = "rickandmorty",
                          "The Simpsons" = "simpsons")