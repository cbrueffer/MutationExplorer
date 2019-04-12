#######################################################
#
# Functions for waterfall plots, based on GenVisR
#
#######################################################

suppressPackageStartupMessages(library(GenVisR))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))


plot.waterfall <- function(input, sample.tbl, mut.tbl) {

    #######################################################
    #
    # prepare data structures for waterfall()
    #
    #######################################################

    # create the clinical annotation table, in long format
    clin.anno <- sample.tbl %>%
        dplyr::select(sample = SAMPLE,
                      NHG,
                      Ki67,
                      HER2,
                      PgR = PgR_1perc,
                      ER = ER_1perc,
                      PAM50,
                      HistType = Histological_Type) %>%
        melt(id.vars = c("sample"))

    # ordering and coloring
    clinicColor <- c(LumA = "blue4", LumB = "deepskyblue", HER2 = "hotpink2", Basal = "firebrick2", Normal = "green4", Unclassified = "gray",
                     Ductal = "yellow", Lobular = "red", Other = "gray",
                     NEG = "white", POS = "black",
                     G1 = "#deebf7", G2 = "#9ecae1", G3 = "#3182bd",
                     `NA` = "gray")
    clinicOrder <- c("LumA", "LumB", "HER2", "Basal", "Normal", "Unclassified",
                     "Ductal", "Lobular", "Other",
                     "G1", "G2", "G3", "NEG", "POS", "NA")


    # mutation data, waterfall() needs the columns "sample", "gene", and "variant_class"
    mutDf <- mut.tbl %>%
        dplyr::select(sample = SAMPLE, gene = gene.symbol)

    # count the occurrence of each mutation in our set
    mut_count <- plyr::count(plyr::count(mutDf, c('gene', 'sample'))[,1:2], 'gene')

    mutDf <- mutDf %>%
        mutate(own_freq = mut_count$freq[match(gene, mut_count$gene)] / length(unique(sample))) %>%
        mutate(gene = paste0(gene, ' (', sprintf("%.0f%%", own_freq * 100), ')')) %>%
        mutate(variant_class = mut.tbl$ANN.effect.class)

    # nomenclature for alteration impact
    mut_impact_order <- c("Frameshift",
                          "Nonsense",
                          "In-frame indel",
                          "Missense",
                          "Splicing",
                          "Synonymous",
                          "UTR",
                          "Other"
    )

    mut_impact_palette <- c("red",
                            "orange",
                            "darkgreen",
                            "#bbebff",  # lightblue
                            "#0092FF",  # blue
                            "#4900FF",  # lila
                            "#FF00DB",  # fuchsia
                            "gray"
    )

    # determine the top X mutated genes
    topX.mut <- mut_count %>%
        arrange(desc(freq)) %>%
        dplyr::select(gene) %>%
        head(input$waterfall.cutoff)
    topX.mut <- as.character(topX.mut$gene)

    # Make sure the mutation status for the topX genes is present in the sample table for figure out sample ordering
    for (i in seq_along(topX.mut)) {
        var <- paste0("mut.status.", topX.mut[i])
        if (!(var %in% colnames(sample.tbl))) {
            sample.tbl[[var]] <- as.factor(ifelse(sample.tbl$SAMPLE %in% mut.tbl[mut.tbl$gene.symbol == topX.mut[i], ]$SAMPLE, "mut", "wt"))
        }
    }

    #
    # Sample arrangement
    #
    # order samples by histological type first, then by frequently mutated genes
    arrange_vars = c("Histological_Type", paste0("mut.status.", topX.mut))

    sample.order = sample.tbl %>%
        arrange_(.dots = arrange_vars) %>%
        dplyr::select(SAMPLE)
    sample.order = as.character(sample.order$SAMPLE)


    #######################################################
    #
    # plot
    #
    #######################################################

    samp_recur_layer = theme(
        axis.text.x = element_text(size = rel(1.3))
    )

    p = waterfall(mutDf,
                  fileType="Custom",
                  variant_class_order = mut_impact_order,
                  mainRecurCutoff = 0.00001,
                  mainDropMut = TRUE,
                  maxGenes = input$waterfall.cutoff,
                  mainPalette = mut_impact_palette,
                  clinDat = clin.anno,
                  clinLegCol = 2,
                  sampOrder = sample.order,
                  main_geneLabSize = 12,
                  clinVarCol = clinicColor,
                  clinVarOrder = clinicOrder,
                  section_heights = c(0.2, 1, 0.35),
                  plotMutBurden = FALSE,
                  mainGrid = F,
                  sampRecurLayer = samp_recur_layer,
                  out = "plot")

    return(p)
}