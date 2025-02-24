#' 
#' Get the stripped cell line name and the tumors ID of the user selected combined matrix
#' 
#' @import ggplot2
#' @import dplyr
#' @param mat_m3 dataframe of lineage frequencies
#' @return lineages heatmap
#' @export



get_confusion_matrix_shiny <- function(mat_m3) {
  
  ccle_heatmap_m3_dist <- ggplot2::ggplot(mat_m3, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
    geom_tile(color = "black") +
    coord_fixed() +
    theme_dark() +
    scale_fill_gradient(low = "white", high = "red") + 
    labs(x = "CCLE", y = "TCGA", title = "Proportion of Agreement") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 12, face = "italic")) + 
    guides(fill = guide_colourbar(barwidth = 0.5,barheight = 8)) 
  
  return(ccle_heatmap_m3_dist)
  
}