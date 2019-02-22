#' sparseSumLOdds
#'
#' This function calculates the sum of the log odds per class, even if not all classifier genes are present.
#' @param stab_lor_list A list of odds ratios annotated with Entrez IDs. Can simply use the output of the coef_stabilizer() function.
#' @param exprs_data A matrix with rows representing genes and columns representing samples.
#' @param intercept A boolean whether to include the intercept
#'
#' @export

sparseSumLOdds <- function(stab_lor_list, exprs_data, intercept){
  lapply(stab_lor_list, function(mac_lor){
    mac_lor_sub <- mac_lor[-which(names(mac_lor) %in% "(Intercept)")]

    ol_genes <- intersect(names(mac_lor_sub), rownames(exprs_data))

    #Classifier genes present
    in_genes <- names(mac_lor_sub)[names(mac_lor_sub) %in% rownames(exprs_data)]
    out_genes <- names(mac_lor_sub)[!names(mac_lor_sub) %in% rownames(exprs_data)]

    presence_stats <- data.frame(logOR_in = sum(abs(mac_lor_sub[in_genes])),
                                 logOR_total = sum(abs(mac_lor_sub)),
                                 logOR_relative = sum(abs(mac_lor_sub[in_genes]))/sum(abs(mac_lor_sub)))

    if(intercept){
      mac_lor_ol <- mac_lor[c("(Intercept)", ol_genes)]
      if(is.null(dim(exprs_data))) lodds <- mac_lor_ol %*% rbind(1, matrix(exprs_data[ol_genes,]))
      else lodds <- mac_lor_ol %*% rbind(1, exprs_data[ol_genes,])
    } else{
      mac_lor_ol <- mac_lor[ol_genes]
      if(is.null(dim(exprs_data))) lodds <- mac_lor_ol %*% matrix(exprs_data[ol_genes,])
      else lodds <- mac_lor_ol %*% exprs_data[ol_genes,]
    }

    return(list(lodds = lodds, in_genes = in_genes, out_genes = out_genes, presence_stats = presence_stats))
  })
}
