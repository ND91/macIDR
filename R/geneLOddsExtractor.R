#' geneLOddsExtractor
#'
#' This function extracts the gene log odds and stores it in a dataframe
#'
#' @param scObj A sparseClassification object as obtained from the macClassifier function.
#' @param hgnc A boolean argument whether you want the log odds with Entrez (FALSE) or HGNC symbols (TRUE). Defaults to TRUE.
#' @param intercept A boolean argument whether you want the intercept to be returned as well. Defaults to TRUE.
#' @param na.rm A boolean argument whether you want to remove absent classifier genes. Defaults to TRUE.
#' @param loddsum A boolean argument on whether you want to sum the log odds per activation state. Defaults to FALSE.
#'
#' @import org.Hs.eg.db
#' 
#' @export
#' 

geneLOddsExtractor <- function(scObj, hgnc = T, intercept = T, na.rm = T, loddsum = F){
  if(is.null(scObj)) stop("No macrophage classifier data found!")
  if(class(scObj) != "sparseClassification") stop("scObj is not a sparseClassification object")
  
  GLO <- getGeneLOdds(scObj)
  
  if(!na.rm){
    GLO_absent <- getAbsent(scObj)
    GLO <- mapply(FUN = function(GLO_entry, GLO_entry_absent){
      GLO_absent_df <- data.frame(matrix(NA, nrow = length(GLO_entry_absent), ncol = ncol(GLO_entry), dimnames = list(GLO_entry_absent, colnames(GLO_entry))))
      GLO_entry <- rbind(GLO_entry, GLO_absent_df)
      return(GLO_entry)
    }, GLO_entry = GLO, GLO_entry_absent = GLO_absent, SIMPLIFY = F)
  }
  
  GLO_rn <- lapply(GLO, function(mac_lodd){
    rownames(mac_lodd)[1] <- "(Intercept)"
    if(loddsum) mac_lodd <- rbind(mac_lodd, Sum = colSums(mac_lodd, na.rm = T))
    mac_lodd <- data.frame(mac_lodd, 
                           entrez = rownames(mac_lodd)) 
    
    if(!intercept) mac_lodd <- mac_lodd[-1,]
    
    mac_lodd$gene <- mapIds(org.Hs.eg.db,
                            keys = rownames(mac_lodd),
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
    
    return(mac_lodd)
  })
  
  GLO_df <- data.frame(do.call(rbind, GLO_rn))
  GLO_df$class <- gsub("^(.+)\\..+$", "\\1", rownames(GLO_df))
  
  return(GLO_df)
}