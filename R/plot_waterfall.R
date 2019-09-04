#######################################################
#
# Functions for waterfall plots, based on GenVisR
#
#######################################################

suppressPackageStartupMessages(library(GenVisR))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))


plot.waterfall <- function(input, sample.tbl, mut.tbl, gene.column.map) {

    #######################################################
    #
    # prepare data structures for waterfall()
    #
    #######################################################

    # XXX find a way to do this inside the pipe below
    if (input$hr.cutoff == "hr.1perc") {
        ER_select = sample.tbl$ER_1perc
        PgR_select = sample.tbl$PgR_1perc
    } else {
        ER_select = sample.tbl$ER_10perc
        PgR_select = sample.tbl$PgR_10perc
    }

    # create the clinical annotation table, in long format
    clin.anno <- sample.tbl %>%
        dplyr::mutate(PgR = PgR_select,
                      ER = ER_select) %>%
        dplyr::select(sample = SAMPLE,
                      NHG,
                      Ki67,
                      HER2,
                      PgR,
                      ER,
                      PAM50,
                      HistType = Histological_Type) %>%
        melt(id.vars = c("sample"))

    # ordering and coloring
    clinicColor <- c(LumA = "blue4", LumB = "deepskyblue", HER2 = "hotpink2", Basal = "firebrick2", Normal = "green4", Unclassified = "gray",
                     Ductal = "yellow", Lobular = "red", Other = "gray",
                     NEG = "white", POS = "black",
                     G1 = "#deebf7", G2 = "#9ecae1", G3 = "#3182bd",
                     `NA` = "gray")
    clinicOrder <- c("Ductal", "Lobular", "Other",
                     "LumA", "LumB", "HER2", "Basal", "Normal", "Unclassified",
                     "G1", "G2", "G3", "NEG", "POS", "NA")


    # mutation data, waterfall() needs the columns "sample", "gene", and "variant_class"
    mutDf <- mut.tbl %>%
        dplyr::select(sample = SAMPLE,
                      gene = gene.symbol,
                      variant_class = ANN.effect.class)

    # count the occurrence of each mutation in our set
    mut_count <- plyr::count(plyr::count(mutDf, c('gene', 'sample'))[, 1:2], 'gene')

    # determine the top X most mutated genes
    topX.mut <- mut_count %>%
        arrange(desc(freq)) %>%
        dplyr::select(gene) %>%
        head(input$waterfall.cutoff)
    topX.mut <- as.character(topX.mut$gene)

    # Make sure the mutation status for the topX genes is present in the sample table for figuring out sample ordering
    for (i in seq_along(topX.mut)) {
        var = gene.column.map[[topX.mut[i]]]
        if (!(var %in% colnames(sample.tbl))) {
            sample.tbl[[var]] <- as.factor(ifelse(sample.tbl$SAMPLE %in% mut.tbl$SAMPLE[mut.tbl$gene.symbol == topX.mut[i]], "mut", "wt"))
        }
    }

    # Restrict the mutation table to the genes we'll actually display.
    mutDf <- mutDf %>%
        mutate(own_freq = mut_count$freq[match(gene, mut_count$gene)] / length(unique(sample))) %>%
        dplyr::filter(gene %in% topX.mut)

    #
    # Effect/color settings
    #
    # Filter effects that don't occur in the displayed genes.
    this.mut.effect.tbl = dplyr::filter(mut.effect.tbl, effect %in% mutDf$variant_class[mutDf$gene %in% topX.mut])

    #
    # Sample arrangement
    #
    # order samples by histological type first, then by frequently mutated genes
    arrange_vars = c("Histological_Type", lapply(topX.mut, function(gene) gene.column.map[[gene]]))

    sample.order = sample.tbl %>%
        arrange_(.dots = arrange_vars)
    sample.order = as.character(sample.order$SAMPLE)

    # Add cohort frequency to the gene symbol for display purposes
    mutDf <- mutate(mutDf, gene = paste0(gene, ' (', sprintf("%.0f%%", own_freq * 100), ')'))

    #######################################################
    #
    # plot
    #
    #######################################################

    samp_recur_layer = theme(
        axis.text.x = element_text(size = rel(1.3))
    )

    p = waterfall(mutDf,
                  fileType = "Custom",
                  variant_class_order = this.mut.effect.tbl$effect,
                  mainRecurCutoff = 0.00001,
                  mainDropMut = FALSE,  # turned off, we do it manually above since this option doesn't drop the corresponding color
                  plotGenes = unique(mutDf$gene),
                  plotSamples = sample.tbl$SAMPLE,
                  mainPalette = this.mut.effect.tbl$color,
                  clinDat = clin.anno,
                  clinLegCol = 2,
                  sampOrder = sample.order,
                  main_geneLabSize = 12,
                  clinVarCol = clinicColor,
                  clinVarOrder = clinicOrder,
                  section_heights = c(0.2, 1, 0.35),
                  plotMutBurden = FALSE,
                  mainGrid = FALSE,
                  sampRecurLayer = samp_recur_layer,
                  out = "grob")

    return(p)
}