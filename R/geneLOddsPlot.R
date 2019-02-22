#' geneLOddsPlot
#'
#' This function visualizes the log odds of the individual genes.
#'
#' @param scObj A sparseClassification object as obtained from the macClassifier function.
#' @param samples (optional) A vector including the samples to be extracted. Uses the column names of the expression data.
#' @param phenotype (optional) A vector of length equal to the number of samples included in the classification containing the actual classes.
#'     Used for coloring the bars.
#' @param title (optional) A title for the plot.
#' @param hgnc A boolean argument whether you want the log odds with Entrez or HGNC symbols. Defaults to HGNC symbols.
#'
#' @import ggplot2
#' 
#' @export
#' 

geneLOddsPlot <- function(scObj, samples = NULL, phenotype = NULL, title = NULL, hgnc = T){
  if(is.null(scObj)) stop("No macrophage classifier data found!")
  if(class(scObj) != "sparseClassification") stop("scObj is not a sparseClassification object")
  if(!is.null(samples) & !is.null(phenotype) & length(samples) != length(phenotype)) stop("The number of requested samples and the associated phenotypes do not match.")
  if(is.null(samples) & !is.null(phenotype) & length(samples) != ncol(scObj)) stop("The number of phenotypes does not match the nuber of samples as found in the scObj.")
  if(is.null(samples)) samples <- names(scObj)
  
  GLO_df <- geneLOddsExtractor(scObj = scObj, hgnc = hgnc, intercept = F)
  GLO_df <- GLO_df[,c(samples, "entrez", "gene", "class")]
  
  GLO_melt <- melt(GLO_df)
  colnames(GLO_melt) <- c("Entrez", "Gene", "Class", "Sample", "Log_odds")
  
  if(!is.null(phenotype)){
    GLO_melt$Group <- rep(phenotype, each = unique(table(GLO_melt$Sample)))
    plot_obj <- ggplot(GLO_melt, aes(x = Gene, y = Log_odds, group = Group, col = Group))  
  } else{
    plot_obj <- ggplot(GLO_melt, aes(x = Gene, y = Log_odds))  
  }
  
  plot_obj <- plot_obj + geom_jitter(stat = "identity", position = position_dodge(width = 1)) +
    facet_wrap(~ Class, scales = "free_x", ncol = 4) +
    theme_bw() +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="bottom")
  
  if(!is.null(title)) plot_obj <- plot_obj + ggtitle(title)
  
  return(plot_obj)
}