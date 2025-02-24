#' 
#' Generate an interactive plot that highlighting the CL and tumor k nearest neighbors for the user query
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
#' @param omics_name name of the shiny selected omics
#' @return an interactive plot that highlighting the CL and tumor k nearest neighbors found by queryKNN
#' @export


find_neighbors_both <- function(combined_mat, reduced_mat, input_sample, k, ann, BNindex, sample_order, omics_name) {
  
  query <- matrix(combined_mat[input_sample, ], nrow = 1)
  rownames(query)[1] <- rownames(combined_mat)[which(rownames(combined_mat) %in% input_sample)]
  out <- queryKNN(BNINDEX = BNindex,query = query, k = k + 1)
  
  combined_mat_1 <- combined_mat[sample_order,]
  sample_neighbors <- rownames(combined_mat_1[out$index,])
  
  top_k_tumors_1 <- as.character(sample_neighbors)
  
  data_res <- reduced_mat %>% t() %>% 
    as.data.frame() %>% 
    mutate("sampleID" = colnames(reduced_mat)) %>% 
    left_join(., ann, by = "sampleID")
  
  colnames(data_res)[c(1,2)] <- c("UMAP_1", "UMAP_2")
  
  data_res_1 <- data_res %>% mutate('show_it' = ifelse(data_res$sampleID %in% c(input_sample,top_k_tumors_1), 'show', 'not'))
  
  data_res_2 <- data_res_1 %>% mutate('size' = if_else(show_it == 'show', 16, 5))
  
  dist <- SNFtool::dist2(matrix(combined_mat[which(rownames(combined_mat) %in% input_sample), ], nrow = 1), combined_mat)
  dist_1 <- reshape::melt(dist) 
  colnames(dist_1)[c(2,3)] <- c('sampleID','dist')
  dist_2 <- dist_1 %>% dplyr::select(sampleID, dist)
  
  data_res_3 <- data_res_2 %>% left_join(dist_2, by = 'sampleID') %>% arrange(dist)
  data_res_3$dist <- round(data_res_3$dist, 3)
  
  data_res_3 <- data_res_3 %>% select(UMAP_1,UMAP_2,stripped_cell_line_name,sampleID,lineage,
                                      subtype,subtype_1,type,dist,show_it,size)
  
  shared <- SharedData$new(data_res_3, key = ~sampleID)
  
  if(omics_name %in% c('MoNETA multiomics ', 'MOFA multiomics ')) {
    
    row_1 <- crosstalk::bscols(
      widths = c(2, 10), 
      list(
        crosstalk::filter_checkbox("Type", 
                                   label = "Select Type",
                                   sharedData = shared, 
                                   group = ~type),
        crosstalk::filter_checkbox("Lineage", 
                                   label = "Select Lineage",
                                   sharedData = shared, 
                                   group = ~lineage),
        crosstalk::filter_select("subtype",
                                 label = "Select Subtype",
                                 sharedData = shared, 
                                 group = ~subtype)), 
      
      plot_ly(
        data = shared,
        x = ~UMAP_1,
        y = ~UMAP_2,
        key = ~sampleID,
        type = 'scatter',
        mode = 'markers',
        color = ~lineage,  
        symbol = ~type, 
        symbols = c('circle',"x"),
        stroke = ~show_it,
        strokes = c('show' = "red"),
        size = ~size,
        sizes = c(5,9),
        hoverinfo = "text",
        hovertext = ~paste("ID:", sampleID,
                           '\nName:', stripped_cell_line_name, 
                           '\nLineage:', lineage,
                           '\nSubtype:', subtype_1,
                           '\nType:', type),
        marker = list(
          line = list(
            width = 3)))
      %>%
        layout(
          title = list(
            text = paste('tSNE projection of', omics_name, 'alignment'), 
            font = list(size = 21, family = "Candara", color = "black", weight = "bold"), 
            x = 0.3,          
            xanchor = "center",  
            yanchor = "top"
          ),
          xaxis = list(title = "UMAP 1", zeroline  = F),
          yaxis = list(title = "UMAP 2", zeroline = F),
          legend = list(
            title = list(text = 'Select lineage-type pair'),
            traceorder = 'normal'),
          height = 600) %>% 
        event_register("plotly_selected") %>% 
        highlight(on = "plotly_selected", off = "plotly_doubleclick",color = 'green', persistent = FALSE)
      )
    
    
    row_2 <- crosstalk::bscols(
      widths = c(10,2),
      
      htmltools::browsable(
        tagList(
          tags$button(
            tagList(fontawesome::fa("download"), "Download neighbors list"),
            onclick = "Reactable.downloadDataCSV('alignment-download-table', 'alignment.csv')"),
          
          reactable(shared$origData()[shared$origData()$show_it == 'show',], searchable = TRUE, minRows = 3, 
                    showPageSizeOptions = TRUE,
                    pageSizeOptions = c(25,30,40,50),
                    defaultPageSize = 25,
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
                      stripped_cell_line_name = colDef(name = 'Name'), 
                      sampleID = colDef(name = "ID"),
                      type = colDef(name = "Type"),
                      lineage = colDef(name = "Lineage"),
                      subtype = colDef(name = "Subtype"),
                      subtype_1 = colDef(name = 'Subtype code'),
                      dist = colDef(name = "Distance"),
                      UMAP_1 = colDef(show = FALSE),
                      UMAP_2 = colDef(show = FALSE),
                      size = colDef(show = FALSE),
                      show_it = colDef(show = FALSE)),
                    
                    style = list(
                      height = "400px",  
                      overflowY = "auto",
                      overflowX = "hidden"
                    )
          ))))
    
    x <- htmltools::browsable(
      htmltools::tagList(row_1, row_2))
    
  } else {
    
    row_1 <- crosstalk::bscols(
      widths = c(2, 10), 
      list(
        crosstalk::filter_checkbox("Type", 
                                   label = "Select Type",
                                   sharedData = shared, 
                                   group = ~type),
        crosstalk::filter_checkbox("Lineage", 
                                   label = "Select Lineage",
                                   sharedData = shared, 
                                   group = ~lineage),
        crosstalk::filter_select("subtype",
                                 label = "Select Subtype",
                                 sharedData = shared, 
                                 group = ~subtype)), 
      
      plot_ly(
        data = shared,
        x = ~UMAP_1,
        y = ~UMAP_2,
        key = ~sampleID,
        type = 'scatter',
        mode = 'markers',
        color = ~lineage,  
        symbol = ~type, 
        symbols = c('circle',"x"),
        stroke = ~show_it,
        strokes = c('show' = "red"),
        size = ~size,
        sizes = c(5,9),
        hoverinfo = "text",
        hovertext = ~paste("ID:", sampleID,
                           '\nName:', stripped_cell_line_name, 
                           '\nLineage:', lineage,
                           '\nSubtype:', subtype_1,
                           '\nType:', type),
        marker = list(
          line = list(
            width = 3)))
      %>%
        layout(
          title = list(
            text = paste('UMAP projection of', omics_name, 'alignment'), 
            font = list(size = 21, family = "Candara", color = "black", weight = "bold"), 
            x = 0.3,          
            xanchor = "center",  
            yanchor = "top"
          ),
          xaxis = list(title = "UMAP 1", zeroline  = F),
          yaxis = list(title = "UMAP 2", zeroline = F),
          legend = list(
            title = list(text = 'Select lineage-type pair'),
            traceorder = 'normal'),
          height = 600) %>%
        event_register("plotly_selected") %>% 
        highlight(on = "plotly_selected", off = "plotly_doubleclick",color = 'green', persistent = FALSE)
      )
    
    
    row_2 <- crosstalk::bscols(
      widths = c(10,2),
      
      htmltools::browsable(
        tagList(
          tags$button(
            tagList(fontawesome::fa("download"), "Download neighbors list"),
            onclick = "Reactable.downloadDataCSV('alignment-download-table', 'alignment.csv')"),
          
          reactable(shared$origData()[shared$origData()$show_it == 'show',], searchable = TRUE, minRows = 3, 
                    showPageSizeOptions = TRUE,
                    pageSizeOptions = c(25,30,40,50),
                    defaultPageSize = 25,
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
                      stripped_cell_line_name = colDef(name = 'Name'), 
                      sampleID = colDef(name = "ID"),
                      type = colDef(name = "Type"),
                      lineage = colDef(name = "Lineage"),
                      subtype = colDef(name = "Subtype"),
                      subtype_1 = colDef(name = 'Subtype code'),
                      dist = colDef(name = "Distance"),
                      UMAP_1 = colDef(show = FALSE),
                      UMAP_2 = colDef(show = FALSE),
                      size = colDef(show = FALSE),
                      show_it = colDef(show = FALSE)),
                    
                    style = list(
                      height = "400px",  
                      overflowY = "auto",
                      overflowX = "hidden"
                    )
          ))))
    
    x <- htmltools::browsable(
      htmltools::tagList(row_1, row_2))
    
  }
  
  return(x)
  
}