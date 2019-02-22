#' macClassifier
#'
#' This function takes the expression data in matrix format and the stabilized log odds ratios (as calculated using the lORStabilizer
#' function) and calculates the response for different macrophage subsets. #' Currently, the expression data can be classified as
#' "CTRL", "LPS_Early", "LPS_Late", "LPS_IFNg", "IFNg", "IL4", "IL10", or "DEX" macrophage. The results are then stored in a
#' sparseClassification object, which can be queried using specific accessor methods.
#'
#' @param exprs_data A matrix with rows representing genes and columns representing samples.
#' @param lor_list A list of log odds ratios per macrophage subset.
#' @param intercept A boolean whether to include the intercept
#' @examples
#'
#' #Prepare the test data
#' data(testset_exprdata)
#' head(testset_exprdata)
#'
#' #Prepare the stabilized odds ratios (defaults to 0.5 median stabilization)
#' stab_lor_list <- lORStabilizer()
#'
#' #Calculate the predictions for the test data per class
#' testset_spobj <- macClassifier(testset_exprdata, stab_lor_list)
#'
#' #Extract the predictions [0,1] (highest value is the predicted class)
#' testset_props <- getProportions(testset_spobj)
#'
#' @export

macClassifier <- function(exprs_data, lor_list, intercept = T){
  if(is.null(exprs_data)) stop("No expression data found!")
  if(is.null(lor_list)) stop("No log odds ratio list found! Please provide run lORStabilizer first.")

  lodd_list <- sparseSumLOdds(stab_lor_list = lor_list, exprs_data = exprs_data, intercept = intercept)

  #Gene-wise log odds
  gene_lodds <- sparseOdds(stab_lor_list = lor_list, exprs_data = exprs_data, expon = F)

  #Sum log odds
  lodds <- t(data.frame(do.call(rbind, lapply(lodd_list, function(lodds){
    lodds$lodds
  })), row.names = names(lodd_list)))

  #Odds
  odds <- exp(lodds)

  #Proportions
  proportions <- odds/rowSums(odds)

  #Predictions
  predictions <- colnames(proportions)[apply(proportions, 1, which.max)]
  names(predictions) <- rownames(proportions)

  #Presence absence statistics
  present_genes <- lapply(lodd_list, function(lodds){lodds$in_genes})
  absent_genes <- lapply(lodd_list, function(lodds){lodds$out_genes})
  presence_stats <- data.frame(do.call(rbind, lapply(lodd_list, function(lodds){
    lodds$presence_stats
  })), row.names = names(lodd_list))

  sparse_classes <- new("sparseClassification",
                        gene_lodds = gene_lodds,
                        lodds = lodds,
                        odds = odds,
                        proportions = proportions,
                        predictions = predictions,
                        presence_stats = presence_stats,
                        present_genes = present_genes,
                        absent_genes = absent_genes)

  return(sparse_classes)
}
