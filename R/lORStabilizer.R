#' lORStabilizer
#'
#' This function stabilizes the coefficients as obtained from the repeated cross-validation data ("repcv_lor").
#'
#' @param threshold The proportion of the iterations that must be non-zero for a gene to be considered stable.
#' @param aggregation_method Which method for aggregation across the iterations should be used? Options are taking the median ("median")
#'     or taking the mean ("mean"). Defaults to "median".
#' @export

lORStabilizer <- function(threshold = 0.5, aggregation_method = c("median", "mean")){
  threshold <- as.numeric(abs(threshold))
  aggregation_method <- match.arg(aggregation_method)

  lapply(repcv_lor, function(mac_lors){
    rimportance <- apply(mac_lors, 1, function(gene){mean(abs(gene)>0)})
    stab_genes <- names(which(rimportance > threshold))

    stab_lors <- switch(aggregation_method,
                         "median" = apply(mac_lors[stab_genes,], 1, median),
                         "mean" = apply(mac_lors[stab_genes,], 1, mean))
    return(stab_lors)
  })
}
