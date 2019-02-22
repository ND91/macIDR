#' sparseOdds
#'
#' This function calculates the (log) odds for the individual genes composing each classifier.
#'
#' @param stab_lor_list A list of odds ratios annotated with Entrez IDs. Can simply use the output of the coef_stabilizer() function.
#' @param exprs_data A matrix with rows representing genes and columns representing samples.
#' @param expon A boolean argument whether or not to return the odds (TRUE) or the log odds (FALSE). Defaults to TRUE.
#' 
#' @export

sparseOdds <- function(stab_lor_list, exprs_data, expon = TRUE){
  or_list <- lapply(stab_lor_list, function(lor_mac){
    features_ol <- intersect(names(lor_mac[-1]), rownames(exprs_data))
    
    if(length(features_ol) != 0){
      if(length(features_ol) != 1){
        exprs_data_sub <- rbind(1, as.matrix(exprs_data[features_ol,]))
      } else{
        exprs_data_sub <- t(cbind(1, as.matrix(exprs_data[features_ol,]))) #Don't ask...
      }
      
      lor_mac_sub <- c(lor_mac[1], lor_mac[features_ol])
      lodds <- apply(exprs_data_sub, 2, function(i){i*lor_mac_sub})
    } else{
      lodds <- NA
    }
    
    ifelse(expon, odds <- exp(lodds), odds <- lodds)
    return(odds)
  })
  return(or_list)
}
