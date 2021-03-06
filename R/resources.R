#
# Constants
#

# Waterfall plots
plot.waterfall.cutoff.min = 1
plot.waterfall.cutoff.max = 30
plot.waterfall.cutoff.default = 20
plot.waterfall.round.digits.min = 0
plot.waterfall.round.digits.max = 3
plot.waterfall.round.digits.default = 0
plot.waterfall.height.min = 100
plot.waterfall.height.max = 5000
plot.waterfall.height.default = 700
plot.waterfall.width.min = 100
plot.waterfall.width.max = 5000
plot.waterfall.width.default = 1200

# Survival plots (gene/pathway/TMB)
plot.survival.height.min = 100
plot.survival.height.max = 5000
plot.survival.height.default = 700
plot.survival.width.min = 100
plot.survival.width.max = 5000
plot.survival.width.default = 700

# Protein plots
plot.protein.cutoff.mutation.min = 0
plot.protein.cutoff.mutation.default = 0
plot.protein.cutoff.annotation.min = 0
plot.protein.cutoff.annotation.default = 10
plot.protein.height.min = 100
plot.protein.height.max = 1000
plot.protein.height.default = 300
plot.protein.width.min = 100
plot.protein.width.max = 5000
plot.protein.width.default = 1000

# Misc
pathway.custom.genes.max = 20
survival.gene.plots.max = 9
survival.pathway.plots.max = 9

#
# Treatment definitions
#
treatment.tbl <- as_tibble(rbind(c("any", "Any Treatment", "Any Treatment"),
                                 c("Adjuvant_Endo", "Endocrine + any", "Endo + any"),
                                 c("Adjuvant_Cyto", "Cyto + any", "Cyto + any"),
                                 c("Adjuvant_HER2", "HER2 + any", "HER2 + any"),
                                 c("Adjuvant_Endo_only", "Endocrine only", "Endo only"),
                                 c("Adjuvant_Cyto_only", "Cyto only", "Cyto only"),
                                 c("Adjuvant_HER2_only", "HER2 only", "HER2 only"),
                                 c("Adjuvant_Endo_and_Cyto", "Endo + Cyto + any", "Endo + Cyto + any"),
                                 c("Adjuvant_Endo_and_HER2", "Endo + HER2 + any", "Endo + HER2 + any"),
                                 c("Adjuvant_Cyto_and_HER2", "Cyto + HER2 + any", "Cyto + HER2 + any"),
                                 c("Adjuvant_Endo_and_Cyto_only", "Endo + Cyto only", "Endo + Cyto only"),
                                 c("Adjuvant_Endo_and_HER2_only", "Endo + HER2 only", "Endo + HER2 only"),
                                 c("Adjuvant_Cyto_and_HER2_only", "Cyto + HER2 only", "Cyto + HER2 only"),
                                 c("Adjuvant_Endo_and_Cyto_and_HER2_only", "Endo + Cyto + HER2", "Endo + Cyto + HER2"),
                                 c("Adjuvant_No_Systemic", "No Systemic Treatment", "Not Systemically Treated")
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
biomarker.tbl <- as_tibble(rbind(c("ER", "ER_1perc", "POS", "Positive", "", 1),
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
                                 c("PAM50", "PAM50", "Unclassified", "Unclassified", "gray", 6),
                                 c("HistType", "Histological_Type", "Ductal", "Ductal", "", 1),
                                 c("HistType", "Histological_Type", "Lobular", "Lobular", "", 2),
                                 c("HistType", "Histological_Type", "Other", "Other", "", 3)
))
colnames(biomarker.tbl) <- c("marker", "var", "level", "label", "color", "selection.value")

#
# Mutation effects
#
mut.effect.tbl <- as_tibble(rbind(c("Frameshift", "red"),
                                  c("Nonsense", "orange"),
                                  c("In-frame indel", "darkgreen"),
                                  c("Missense", "#bbebff"),
                                  c("Splicing", "#0092FF"),
                                  c("Synonymous", "#4900FF"),
                                  c("UTR", "#FF00DB"),
                                  c("Other", "gray")
))
colnames(mut.effect.tbl) <- c("effect", "color")


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
        biomarker.subset = filter(biomarker.tbl, marker == this.marker)

        marker.list[["Any"]] <- 42  # every UI element has an "any" option
        for (this.level in biomarker.subset$level) {
            this.biomarker = filter(biomarker.subset, level == this.level)
            this.selection.value = as.integer(dplyr::select(this.biomarker, selection.value))
            this.label = as.character(dplyr::select(this.biomarker, label))
            marker.list[[this.label]] = this.selection.value
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
        biomarker.subset = filter(biomarker.tbl, marker == this.marker)

        marker.list[["42"]] <- "Any"  # every UI element has an "any" option
        for (this.level in biomarker.subset$level) {
            this.biomarker = filter(biomarker.subset, level == this.level)
            this.selection.value = as.character(dplyr::select(this.biomarker, selection.value))
            marker.list[[this.selection.value]] = this.level
        }
        result.list[[this.marker]] = marker.list
    }
    return(result.list)
}
selection.to.label.list <- get.selection.to.biomarker.levels(biomarker.tbl)


# Generates PAM50 HTML labels
get.pam50.html.labels <- function(biomarker.tbl) {
    biomarker.tbl.pam50 <- filter(biomarker.tbl, marker == "PAM50")
    result.list = list()
    result.list[["Any"]] = "<div style='background: white; color: black; font-size: 110%; padding-left: 5px;'>Any</div>"
    for (i in 1:nrow(biomarker.tbl.pam50)) {
        subtype = biomarker.tbl.pam50$level[i]
        label = biomarker.tbl.pam50$label[i]
        color = biomarker.tbl.pam50$color[i]
        result.list[[subtype]] = paste0("<div style='background: ", color, "; color: white; font-size: 110%; padding-left: 5px;'>", label, "</div>")
    }
    return(result.list)
}
pam50.html.labels = get.pam50.html.labels(biomarker.tbl)

# Generates mutation effect HTML labels
get.mut.effect.html.labels <- function(effect.tbl) {
    result.list = list()
    for (i in 1:nrow(effect.tbl)) {
        effect = effect.tbl$effect[i]
        color = effect.tbl$color[i]
        result.list[[effect]] = sprintf('<span style="color: %s">%s</span>', color, effect)
    }
    return(result.list)
}
mut.effect.html.labels = get.mut.effect.html.labels(mut.effect.tbl)

mutation.selection.options = c("All" = "mutations.all",
                               "COSMIC only" = "mutations.cosmic")

plot.type.options = c("Mutational Landscape" = "mut.waterfall.plot",
                      "Survival by Mutated Genes" = "mut.gene.plot",
                      "Survival by Mutated Pathways" = "mut.pathway.plot",
                      "Survival by Mutational Burden" = "mut.burden.plot",
                      "Protein Plot" = "mut.protein.plot")

pathway.type.options = c("Reactome" = "pathway.reactome",
                         "Custom" = "pathway.custom")

tmb.type.options = c("Normalized TMB" = "tmb.normalized",
                     "Absolute TMB" = "tmb.absolute")

hr.cutoff.options = c("1%" = "hr.1perc",
                      "10%" = "hr.10perc")

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
