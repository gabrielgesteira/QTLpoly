#' QTLs with respective support interval plots
#'
#' Creates a plot where colored bars represent the support intervals for QTL peaks (black dots).
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @param main a character string with the main title; if \code{NULL}, no title will be shown.
#'
#' @param drop if \code{TRUE}, phenotypes with no QTL will be dropped; if \code{FALSE} (default), all phenotypes will be shown.
#'
#' @return A \pkg{ggplot2} with QTL bars for each linkage group.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}, \code{\link[qtlpoly]{profile_qtl}}
#'
#' @examples
#'   \donttest{
#'   # Estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   library(qtlpoly)
#'   genoprob4x = lapply(maps4x[c(5)], calc_genoprob)
#'   data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)
#'
#'   # Search for QTL
#'   remim.mod = remim(data = data, pheno.col = 1, w.size = 15, sig.fwd = 0.0011493379,
#' sig.bwd = 0.0002284465, d.sint = 1.5, n.clusters = 1)
#'
#'   # Plot support intervals
#'   plot_sint(data = data, model = remim.mod)
#'   }
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#' @export plot_sint
#' @import ggplot2

plot_sint <- function(data, model, pheno.col=NULL, main=NULL, drop=FALSE) {
  trait <- c(rep("LGs", data$nlgs))
  lg <- c(1:data$nlgs)
  lg.lab <- c(1:data$nlgs)
  lower <- c(rep(0, data$nlgs))
  pos <- c(rep(-10, data$nlgs))
  upper <- c(data$lgs.size)
  
  if(is.null(pheno.col)) pheno.col <- model$pheno.col
  nphen <- length(pheno.col)
  for(p in 1:nphen) {
    t <- which(model$pheno.col == pheno.col[p])
    nqtl <- dim(model$results[[t]]$qtls)[1]
    if(!is.null(nqtl)) {
      trait <- c(trait, rep(names(model$results)[[t]], nqtl), rep(names(model$results)[[t]], data$nlgs))
      lg <- c(lg, model$results[[t]]$qtls[,1], c(1:data$nlgs))
      lg.lab <- c(lg.lab, rep(" ", nqtl), rep(" ", data$nlgs))
      lower <- c(lower, model$results[[t]]$lower[,2], rep(-10,data$nlgs))
      pos <- c(pos, model$results[[t]]$qtls[,2], rep(-10,data$nlgs))
      upper <- c(upper, model$results[[t]]$upper[,2], rep(-10,data$nlgs))
    } else if(!drop) {
      trait <- c(trait, rep(names(model$results)[[t]], data$nlgs))
      lg <- c(lg, c(1:data$nlgs))
      lg.lab <- c(lg.lab, rep(" ", data$nlgs))
      lower <- c(lower, c(rep(-10, data$nlgs)))
      pos <- c(pos, c(rep(-10, data$nlgs)))
      upper <- c(upper, c(rep(-10, data$nlgs)))
    }
  }
  
  DF <- data.frame(trait=trait, lg=as.integer(lg), lg.lab=lg.lab, lower=lower, pos=pos, upper=upper)
  DF$trait <- factor(DF$trait, levels = unique(DF$trait))
  DF$lg <- factor(DF$lg, levels = unique(DF$lg))
  trait.color <- c("black", hcl(h = seq(15, 375, length = nlevels(DF$trait)), l = 65, c = 100)[1:nlevels(DF$trait)])
  bxp.table <- table(interaction(DF$trait, DF$lg))
  bxp.width <- c(bxp.table[unlist(lapply(interaction(DF$trait, DF$lg), function(x) which(names(bxp.table) %in% x)))])
  bxp.width <- ifelse(bxp.width == 1, 1, bxp.width-.25)
  
  p <- ggplot(data = DF) +
    facet_grid(~ lg, scales = "free_x", space = "free_x") +
    geom_boxplot(aes(x=trait, ymin = lower, lower = lower, middle = pos, upper = upper, ymax = upper, color=trait, fill = trait, group = interaction(trait, pos)), position = position_dodge(0), stat="identity", size = 0, width = bxp.width) +
    geom_crossbar(aes(x=trait, y=pos, ymin = -10, ymax = -10), width = 1, position=position_dodge(0)) +
    scale_fill_manual(values = trait.color) +
    scale_color_manual(values = trait.color) +
    coord_cartesian(ylim = c(max(upper), 5)) +
    labs(title=main, x = "Linkage Group", y = "Position (cM)")+
    scale_x_discrete(position="top") +
    scale_y_reverse(expand = c(0,8)) +
    theme_bw() +
    theme(legend.title=element_blank(), axis.text.x=element_blank(), legend.position="bottom",
          title=element_text(face="bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x=element_blank(),
          panel.spacing = unit(0.2, "lines"), strip.text.x = element_text(size = 10), plot.title = element_text(face="bold", hjust = 0.5))
  suppressWarnings(print(p))
}
