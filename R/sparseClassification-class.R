#' An S4 class to represent the sparseClassification object.
#'
#' @slot gene_lodds A list of vectors, where each vector represents the log odds ratios for the classifier gene per macrophage subset.
#' @slot lodds A matrix of summed log odds per macrophage subset.
#' @slot odds A matrix of summed odds per macrophage subset.
#' @slot proportions A matrix of responses [0-1] per macrophage subset.
#' @slot predictions A data frame of predictions per sample.
#' @slot presence_stats A data frame containing the presence and absence statistics of the classifier genes per macrophage subset.
#' @slot present_genes A list of genes present per macrophage subset.
#' @slot absent_genes A list of genes absent per macrophage subset.
#'
#' @export
#'
#' @include allGenerics.R

setClass(Class = "sparseClassification",
         representation = representation(
           gene_lodds = "list",
           lodds = "matrix",
           odds = "matrix",
           proportions = "matrix",
           predictions = "character",
           presence_stats = "data.frame",
           present_genes = "list",
           absent_genes = "list")
         )

setMethod(f = "length", signature = "sparseClassification", definition = function(x){
  return(nrow(x@lodds))
})

setMethod(f = "dim", signature = "sparseClassification", definition = function(x){
  return(dim(x@lodds))
})

setMethod(f = "names", signature = "sparseClassification", definition = function(x){
  return(rownames(x@lodds))
})

setMethod(f = "getGeneLOdds", signature = "sparseClassification", definition = function(object){
  return(object@gene_lodds)
})

setMethod(f = "getLOdds", signature = "sparseClassification", definition = function(object){
  return(object@lodds)
})

setMethod(f = "getOdds", signature = "sparseClassification", definition = function(object) {
  return(object@odds)
})

setMethod(f = "getProportions", signature = "sparseClassification", definition = function(object){
  return(object@proportions)
})

setMethod(f = "getPredictions", signature = "sparseClassification", definition = function(object){
  return(object@predictions)
})

setMethod(f = "getPresent", signature = "sparseClassification", definition = function(object){
  object@present_genes
})

setMethod(f = "getAbsent", signature = "sparseClassification", definition = function(object){
  object@absent_genes
})

setMethod(f = "getPresenceStats", signature = "sparseClassification", definition = function(object){
  object@presence_stats
})
