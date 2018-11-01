suppressPackageStartupMessages(library(survminer))
source("R/utilities.R")
source("R/ggsurvplot_core.R")
source("R/arrange_ggsurvplots.R")

surv.plot <- function(input, fit, data, gene=NULL, legend.labs=NULL, title="") {

    if (is.null(legend.labs)) {
        # fix labels in indidual gene plots
        if (any(grepl("mut.var", names(fit$strata)))) {
            gene = ifelse(is.null(gene), "gene", gene)
            legend.labs = paste0(gene, gsub(".*(=.+)", "\\1", names(fit$strata)))
        } else if (any(grepl("tumor_mutational_burden", names(fit$strata)))) {
                legend.labs = paste0("TMB", gsub(".*(=.+)", "\\1", names(fit$strata)))
        } else {
            legend.labs = names(fit$strata)
        }
    }

    # configure legend location
    if (!isTRUE(input$showLegend)) {
        legend = "none"
    } else if (input$legendLoc == "custom") {
        legend = c(as.numeric(input$legend.coord.x), as.numeric(input$legend.coord.y))
    } else {
        legend = input$legendLoc
    }

    plot <- ggsurvplot(
        fit,                     # survfit object with calculated statistics.
        data = data,
        title = title,

        # risk table
        risk.table = input$show.risk.table,
        risk.table.y.text.col = TRUE, # colour risk table text annotations.
        risk.table.y.text = FALSE,    # show bars instead of names in text annotations
        risk.table.height = 0.17,
        tables.theme = theme_cleantable() +
            theme(
                plot.title = element_text(size=13, face="plain", color="black", hjust=0, margin=margin(b=-0.5, unit="pt"))
            ),
        tables.y.text = FALSE,
        fontsize = 4,            # risk table font size

        # p-value
        pval = input$show.pval,             # show p-value of log-rank test.
        pval.size = 4,
        pval.coord = c(0.13, 0.3),
        #pval.method = TRUE,
        #pval.method.size = 4,
        #pval.method.coord = c(0.13, 0.38),

        # censor details
        censor = input$show.censors,
        censor.shape="|",
        censor.size = 4,

        # confidence intervals
        conf.int = input$show.conf.int,
        #conf.int.style = "step",

        # point estimates of survival curves.
        xlab = "Time after Diagnosis (years)",
        ylab = "Overall Survival (proportion)",
        ggtheme = theme_classic() +
            theme(
                axis.title = element_text(size=13, face="plain", color="black"),
                axis.text = element_text(size=13, face="plain", colour="black"),
                axis.ticks.length=unit(.15, "cm"),
                legend.text=element_text(size=11, face="plain"),
                plot.title = element_text(size="16", face="bold", colour="black", hjust=0.5, margin=margin(b=10, unit="pt"))
            ),
        palette = input$color.palette,
        break.time.by = 1,

        # in legend of risk table
        legend = legend,
        legend.title = input$legendTitle,
        legend.labs = legend.labs
    )

    return(plot)
}