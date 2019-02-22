#' The repeated log Odds Ratios
#'
#' A list containing the 500 times repeated 10-fold cross-validation log Odds Ratios for M0, MLPSearly, MLPSlate, MLPSIFNg, MIFNg, MIL4, MIL10, and Mdex macrophages.
#'
#' @import Matrix
#' @format A list of 8 sparse matrices each of which consists of 5987 rows and 500 columns representing the genes (Entrez ID) and iterations, respectively.
#' \describe{
#'   \item{M0}{log OR for M0 macrophages}
#'   \item{MLPSearly}{log OR for LPS early (2-4 hours) macrophages}
#'   \item{MLPSlate}{log OR for LPS late (>24 hours) macrophages}
#'   \item{MLPSIFNg}{log OR for LPS late with IFNg macrophages}
#'   \item{MIFNg}{log OR for IFNg macrophages}
#'   \item{MIL4}{log OR for IL4 macrophages}
#'   \item{MIL10}{log OR for IL10 macrophages}
#'   \item{Mdex}{log OR for dex macrophages}
#'   ...
#' }
#' @examples
#' data(repcv_lor)
#' @source \url{https://github.com/ND91/PRJ0000004_MACMETA}
"repcv_lor"
