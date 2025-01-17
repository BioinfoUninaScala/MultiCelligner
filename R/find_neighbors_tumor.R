#' 
#' Generate an interactive plot that highlighting the tumor k nearest neighbors for the user query
#'
#' @import plotly
#' @import crosstalk
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @import BiocNeighbors
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param input_sample samples of TCGA or CCLE choosed by the user 
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param BNindex a BiocNeighborIndex object containing precomputed index information
#' @param sample_order index of the original combined matrix
#' @return an interactive plot that highlighting the tumor k nearest neighbors found by queryKNN
#' @export

find_neighbors_tumor <- function(combined_mat, reduced_mat, input_sample, k, ann, BNindex, sample_order) {

  query <- matrix(combined_mat[which(rownames(combined_mat) %in% input_sample), ], nrow = 1)
  out <- queryKNN(BNINDEX = BNindex,query = query, k = k + 100)
  
  combined_mat_1 <- combined_mat[sample_order,]
  sample_neighbors <- rownames(combined_mat_1[out$index,])
  sample_neighbors_1 <- sample_neighbors[grepl('TCGA', sample_neighbors)]
  sample_neighbors_2 <- sample_neighbors_1[1:k]

  top_k_tumors_1 <- as.character(sample_neighbors_2)
  
  data_res <- reduced_mat %>% t() %>% 
    as.data.frame() %>% 
    mutate("sampleID" = colnames(reduced_mat)) %>% 
    left_join(., ann, by = "sampleID")
  
  colnames(data_res)[c(1,2)] <- c("UMAP_1", "UMAP_2")
  
  data_res_1 <- data_res %>% mutate('show_it' = ifelse(data_res$sampleID %in% c(input_sample,top_k_tumors_1), 'show', 'not'))
  
  data_res_2 <- data_res_1 %>% mutate('size' = if_else(show_it == 'show', 16, 5))
  #data_res_2 <- data_res_2 %>% mutate('show_it_CL' = if_else(which(data_res_2$sampleID == input$both_sample,),
  
  shared <- SharedData$new(data_res_2)
  
  row_1 <- crosstalk::bscols(
    widths = c(2, 10), 
    list(
      crosstalk::filter_checkbox("Type", 
                                 label = "type",
                                 sharedData = shared, 
                                 group = ~type),
      crosstalk::filter_checkbox("Lineage", 
                                 label = "lineage",
                                 sharedData = shared, 
                                 group = ~lineage),
      crosstalk::filter_select("subtype",
                               label = "subtype",
                               sharedData = shared, 
                               group = ~subtype)), 
    
    plot_ly(
      data = shared,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = 'scatter',
      mode = 'markers',
      color = ~lineage,  
      symbol = ~type, 
      symbols = c('circle',"x"),
      stroke = ~show_it,
      strokes = c('show' = "red"),
      size = ~size,
      sizes = c(5,9),
      marker = list(
        line = list(
          width = 3)))
    %>%
      layout(
        xaxis = list(title = "UMAP 1", zeroline  = F),
        yaxis = list(title = "UMAP 2", zeroline = F),
        legend = list(
          title = list(text = 'Lineage'),
          traceorder = 'normal'),
        height = 600))
  
  
  row_2 <- crosstalk::bscols(
    widths = c(10,2),
    
    htmltools::browsable(
      tagList(
        tags$button(
          tagList(fontawesome::fa("download"), "Download"),
          onclick = "Reactable.downloadDataCSV('alignment-download-table', 'alignment.csv')"),
        
        reactable(shared$origData()[shared$origData()$show_it == 'show',], searchable = TRUE, minRows = 3, 
                  showPageSizeOptions = TRUE,
                  pageSizeOptions = c(10, 20, 30),
                  defaultPageSize = 10,
                  resizable = TRUE, highlight = TRUE, 
                  selection = "multiple",
                  onClick = "select",
                  theme = reactableTheme(
                    headerStyle = list(
                      "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                      "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"),
                      borderColor = "#555")),
                  bordered = TRUE,
                  striped = TRUE,
                  filterable = TRUE,
                  elementId = "alignment-download-table",
                  columns = list(
                    UMAP_1 = colDef(show = FALSE),
                    UMAP_2 = colDef(show = FALSE),
                    size = colDef(show = FALSE),
                    show_it = colDef(show = FALSE),
                    sampleID = colDef(name = "sampleID"),
                    lineage = colDef(name = "lineage"),
                    subtype = colDef(name = "subtype"),
                    type = colDef(name = "type")),))))
  
  x <- htmltools::browsable(
    htmltools::tagList(row_1, row_2))
  
  return(x)
  
}

