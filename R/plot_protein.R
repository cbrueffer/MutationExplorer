#######################################################
#
# Functions for protein plots
#
#######################################################

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))


plot.protein <- function(input, mutation_df, gene_protein_id_map) {

    # XXX make configurable later
    input_gene = input$protein.plot.gene
    cutOff_Mut = input$protein.plot.mutation.cutoff
    cutOff_Anno = input$protein.plot.anno.cutoff
    show_legend = TRUE
    break_divider = 4

    if (input_gene == "") return()

    # Retrieve the protein domain architecture.
    protein <- gene_protein_id_map$Protein[gene_protein_id_map$Gene == input_gene]
    pfamUrl <- paste0("http://pfam.xfam.org/protein/", protein, "/graphic")
    pfamGraphicsResponse <- GET(pfamUrl)

    # Alteration color definitions
    color.nonsense = "#7C7C7C"
    color.synonymous = "#4EF758"
    color.other = "black"

    # const colors vector
    constColors <- c("Nonsense" = color.nonsense, "Synonymous" = color.synonymous, "Other" = color.other)  # 5E5F5E
    # Color palette for missense mutations
    missense_palette <- "YlOrRd"  # was YlOrRd, RdYlBu

    pfamGraphics_json <- fromJSON(content(pfamGraphicsResponse, "text"))
    # extract text, color and positions for regions
    parsedGraphics <- pfamGraphics_json$regions[[1]] %>%
        dplyr::select(text, colour, start, end)

    proteinLen <- pfamGraphics_json$length

    # get aminoChanges for the gene
    aminoChanges <- mutation_df %>%
        filter(!is.na(ANN.prot.change)) %>%
        mutate(ANN.effect.class.lolli = ifelse(ANN.effect.class %in% c("Missense", "Nonsense", "Synonymous"), ANN.effect.class, "Other")) %>%
        mutate(AA_Change_s = gsub("p\\.(.+)", "\\1", ANN.prot.change)) %>%
        select(gene = gene.symbol,
               TYPE,
               AA_Change_s,
               changeType = ANN.effect.class.lolli,
               proteinPosition = ANN.prot.change.aa) %>%
        filter(gene == input_gene, TYPE == "SNV") %>%
        group_by(proteinPosition) %>%
        mutate(nMutPerPos = n()) %>%
        group_by(AA_Change_s, proteinPosition) %>%
        mutate(nMutPerPosCat = n()) %>%
        ungroup() %>%
        distinct() %>%
        group_by(proteinPosition) %>%
        arrange(proteinPosition, desc(nMutPerPosCat)) %>%
        mutate(y_end = cumsum(nMutPerPosCat), y_start = y_end - nMutPerPosCat,
               category_num = row_number(), category_n = n()) %>%
        ungroup() %>%
        # assign colors based on change type, and final labels only for highly mutated positions
        mutate(catColor = if_else(changeType == "Synonymous", constColors["Synonymous"],
                                  if_else(changeType == "Other", constColors["Other"],
                                          if_else(changeType == "Nonsense", constColors["Nonsense"],
                                                  rev(brewer.pal(7, missense_palette))[category_num]))),
               AA_Change_label = if_else(nMutPerPos >= cutOff_Anno, AA_Change_s, ""))

    # cut-off or different transcript comment
    if (cutOff_Mut > 0 | proteinLen < max(aminoChanges$proteinPosition)) {
        cutOffComment <- paste0(sum(aminoChanges$nMutPerPos <= cutOff_Mut | aminoChanges$proteinPosition > proteinLen),
                                " variants not shown due to cutoff or different transcript")
    } else {
        cutOffComment <- " "
    }

    # final filtering
    aminoChanges <- aminoChanges %>%
        filter(nMutPerPos >= cutOff_Mut, proteinPosition <= proteinLen)
    maxMutations <- max(aminoChanges$nMutPerPos)

    #
    #  Build the top part of the plot, containing mutation bars and labels.
    #
    p <- ggplot(aminoChanges, aes(x = proteinPosition, y = y_start, xend = proteinPosition, yend = y_end)) +
        geom_segment(size = 1,
                     lineend = "butt",
                     color = aminoChanges$catColor) +
        geom_text_repel(aes(x = proteinPosition, y = y_end, label = AA_Change_label),
                        size = 6,
                        segment.alpha = 0.5,
                        alpha = 0.8,
                        color = aminoChanges$catColor,
                        direction = "both",
                        ylim = c(2.1, max(aminoChanges$nMutPerPos))) +
        scale_x_continuous(expand = expand_scale(mult = c(0.01, 0.05)),
                           limits = c(0, proteinLen)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0)),
                           name = substitute(paste(italic(input_gene), " Mutations")),
                           breaks = c(0, maxMutations),
                           labels = c(0, maxMutations)) +
        theme(axis.title.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = rel(1.4)),
              panel.grid = element_blank(),
              panel.background = element_blank(),
              plot.margin=unit(c(7, 1, 0.5, 1), "mm")) +
        annotate("text",
                 Inf,
                 Inf,
                 label = cutOffComment,
                 hjust = 1,
                 vjust = -1,
                 color = "red")
    # for single character alignment
    if (maxMutations < 10) {
        p <- p + theme(axis.title.y = element_text(margin = margin(t = 1, r = 9, b = 1, l = 3), size = rel(1.4)))
    } else {
        p <- p + theme(axis.title.y = element_text(margin = margin(t = 1, r = 3, b = 1, l = 3), size = rel(1.4)))
    }
    # Disable clip-area.
    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"

    # add fake data for legend creation
    fakeD <- data.frame(const = c("Nonsense", "Synonymous", "Other"),
                        Missense = 1:3)

    # calculate x axis breaks and labels
    prot.breaks = c(seq(0, proteinLen, floor(proteinLen / 100) * 100 / break_divider), proteinLen)
    if (prot.breaks[length(prot.breaks)] - prot.breaks[length(prot.breaks) - 1] < 25) {
        prot.breaks.labels = prot.breaks[-c(length(prot.breaks) - 1)]
        prot.breaks.labels = append(prot.breaks.labels, "", after = length(prot.breaks) - 2)
    } else {
        prot.breaks.labels = prot.breaks
    }

    #
    # Build the bottom part of the plot containing the protein domains and legend.
    #
    pb <- ggplot(fakeD, aes(x = c(0, proteinLen), y = c(0, 1))) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.line.x = element_line(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size = rel(1.4)),
              panel.grid = element_blank(),
              panel.background = element_blank(),
              plot.margin = unit(c(0, 1, 0, 1), "mm")) +
        scale_x_continuous(expand = expand_scale(mult = c(0.01, 0.05)),
                           breaks = prot.breaks,
                           labels = prot.breaks.labels) +
        scale_y_continuous(limits = c(0.15, 0.9),
                           expand = c(0, 0)) +
        annotate("rect",
                 xmin = 0, xmax = proteinLen,
                 ymin = 0.45, ymax = 0.75, fill = "grey", alpha = 0.8)
    if (show_legend) {
        pb <- pb +
            geom_tile(aes(1:3, Missense, color = Missense, fill = factor(const, levels = c("Nonsense", "Synonymous", "Other")))) +
            scale_fill_manual(values = constColors) +
            scale_colour_gradient(low = brewer.pal(7, missense_palette)[5], high = brewer.pal(7, missense_palette)[1]) +
            guides(fill = guide_legend(title = element_blank(),
                                       label.theme = element_text(size = 16),
                                       order = 2),
                   color = guide_colourbar(barwidth = 1.5,
                                           barheight = 1,
                                           label = F,
                                           ticks = F,
                                           direction = "horizontal",
                                           title.position = "right",
                                           title.vjust = 0.75,
                                           title.theme = element_text(size = 16),
                                           order = 1)) +
            theme(legend.position = "bottom",
                  legend.margin = margin(t = 0, r = 0, b = 0.1, l = 0, unit = "mm"),
                  legend.justification = "center")
    }

    for (rowI in 1:nrow(parsedGraphics)) {
        pb <- pb +
            annotate("rect",
                     xmin = parsedGraphics[rowI, "start"],
                     xmax = parsedGraphics[rowI, "end"],
                     ymin = 0.3,
                     ymax = 0.9,
                     fill = parsedGraphics[rowI, "colour"],
                     alpha = 0.95) +
            annotate("text",
                     x = (parsedGraphics[rowI, "end"] + parsedGraphics[rowI, "start"])/2,
                     y = 0.6,
                     label = parsedGraphics[rowI, "text"],
                     color = "white",
                     size = 6)
    }

    if (show_legend) {
        fp <- plot_grid(gt, pb, nrow = 2, align = "v", rel_heights = c(2, 1))
    } else {
        fp <- plot_grid(gt, pb, nrow = 2, align = "v", rel_heights = c(4, 1))
    }
    return(fp)
}