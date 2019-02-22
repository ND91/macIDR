#' responsePlot
#'
#' This function acts as a wrapper to plot the responses.
#'
#' @param scObj A sparseClassification object as obtained from the macClassifier function.
#' @param actual_class An optional vector of length equal to the number of samples included in the classification containing the actual classes.
#'     Used for coloring the bars.
#' @param type_response Should the response be a proportion ("proportion"; [0;1]), odds ("odds"; [-Inf;Inf]), log odds ("lodds"; [0;Inf]),
#'     or should all aforementioned values be returned as separate list elements ("all"). Defaults to "proportion".
#'
#' @import ggplot2
#' @export

responsePlot <- function(scObj, actual_class = NULL, type_response = c("proportion", "odds", "lodds")){
  if(is.null(scObj)) stop("No macrophage classifier data found!")
  if(class(scObj) != "sparseClassification") stop("scObj is not a sparseClassification object")
  if(!is.null(actual_class) & length(actual_class) != length(scObj)) stop("scObj and actual_class have to be equal in length")

  type_response <- match.arg(type_response)

  response <- switch(type_response,
                     "proportion" = getProportion(scObj),
                     "odds" = getOdds(scObj),
                     "lodds" = getLodds(scObj))
  response <- data.frame(response,
                         Predicted = unlist(colnames(response)[apply(response, 1, which.max)]),
                         SampleID = rownames(response))

  if(!is.null(actual_class)){
    response <- data.frame(response,
                           Actual = actual_class)
    response_melt <- reshape2::melt(response, id = c("Predicted", "SampleID", "Actual"), variable.name = "Class", value.name = "Response")
  } else{
    response_melt <- reshape2::melt(response, id = c("Predicted", "SampleID"), variable.name = "Class", value.name = "Response")
  }
  
  response_melt$pCol <- NA
  response_melt$pCol[as.character(response_melt$Predicted) == as.character(response_melt$Class)] <- "Predicted"

  if(!is.null(actual_class)){
    response_melt$aCol <- NA
    response_melt$aCol[as.character(response_melt$Actual) == as.character(response_melt$Class)] <- as.character(response_melt$Actual)[as.character(response_melt$Actual) == as.character(response_melt$Class)]
  }

  gplot_obj <- ggplot(response_melt, aes(x = Class, y = Response)) +
    geom_bar(stat = "identity", position = "dodge", aes(fill = pCol))
  
  if(!is.null(actual_class)){
    gplot_obj <- gplot_obj + geom_point(aes(x = aCol, y = 0.5))
  } 
  gplot_obj <- gplot_obj + 
    facet_wrap(~SampleID, ncol = 1, strip.position="right") +
    ylab(type_response) +
    theme_bw() +
    ylim(0,1) +
    ylab(NULL) +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          legend.position="bottom")


  print(gplot_obj)

  return(gplot_obj)
}
