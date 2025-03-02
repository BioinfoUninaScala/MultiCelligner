#' 
#' Create an interactive plot that highlighting k nearest neighbors for the choosed query
#'
#' @import plotly
#' @import crosstalk
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @import reshape
#' @import reshape2 
#' @import SNFtool
#' @import tidyverse
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param omics_name name of the shiny selected omics
#' @return an interactive plot that highlighting the tumor k nearest neighbors 
#' @export

find_neighbors <- function(combined_mat, reduced_mat, input_sample = NULL, selected_samples = NULL, type, k, ann, omics_name) {
  
  if(is.null(selected_samples)) {
  
  x_1 <- combined_mat[which(rownames(combined_mat) %in% input_sample),]
  x_2 <- as.matrix(x_1)
  x_3 <- t(x_2)
  rownames(x_3) <- input_sample
  
  dist <- SNFtool::dist2(x_3, combined_mat)
  dist_1 <- reshape2::melt(dist)
  colnames(dist_1)[c(1,2,3)] <- c('ref_ID', 'sampleID', 'dist')
  
  if ("Tumors" %in% type) { 
  
  dist_2 <- dist_1 %>% 
    arrange(dist) %>%
    mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
    filter(grepl('TCGA', x = sampleID) & priority == 2) 
  }
  
  if ("Cell lines" %in% type) {
    
    dist_2 <- dist_1 %>% 
      arrange(dist) %>%
      mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
      filter(grepl('TCGA', x = sampleID) & priority == 2) 
  }
  
  if (all(c("Cell lines", "Tumors") %in% type)) {
    
    dist_2 <- dist_1 %>% 
      arrange(dist) %>%
      mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
      filter(priority == 2)
  }
  
  dist_top_n <- dist_2 %>% 
    group_by(ref_ID) %>%                  
    arrange(dist) %>%  
    slice_head(n = k) %>% 
    ungroup()
  
  top_k_tumors_1 <- as.character(dist_top_n$sampleID)
  
  } else if (is.null(input_sample)) {
    
    x_1 <- combined_mat[which(rownames(combined_mat) %in% selected_samples),]
    x_2 <- as.data.frame(colMeans(x_1))
    x_3 <- t(x_2)
    rownames(x_3) <- 'metasample'
    mat_metasample <- as.matrix(x_3)
    
    dist_metasample <- SNFtool::dist2(mat_metasample, combined_mat)
    dist_1 <- reshape2::melt(dist_metasample)
    colnames(dist_1)[c(1,2,3)] <- c('metasample','sampleID','dist')
    
    if ("Tumors" %in% type) { 
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(grepl('TCGA', x = sampleID) & priority == 2)
    } 
    
    if ("Cell lines" %in% type) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(grepl('ACH-00', x = sampleID) & priority == 2)
    }
    
    if (all(c("Cell lines", "Tumors") %in% type)) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(priority == 2)
    }
    
    
    dist_top_n <- dist_2 %>% 
      group_by(metasample) %>%                  
      arrange(dist) %>%  
      slice_head(n = k) %>% 
      ungroup()
    
    top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
    
  }
  
  data_res <- reduced_mat %>% t() %>% 
    as.data.frame() %>% 
    mutate("sampleID" = colnames(reduced_mat)) %>% 
    left_join(., ann, by = "sampleID")
  
  colnames(data_res)[c(1,2)] <- c("UMAP_1", "UMAP_2")
  
  data_res_1 <- data_res %>% mutate('show_it' = ifelse(data_res$sampleID %in% top_k_tumors_1, 'show', 'not'))
  
  data_res_2 <- data_res_1 %>% mutate('size' = if_else(show_it == 'show', 16, 5))
  
  dist_metasample_2 <- dist_1 %>% dplyr::select(sampleID, dist)
  data_res_3 <- data_res_2 %>% left_join(dist_metasample_2, by = 'sampleID')  %>% arrange(dist)
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
        source = 'A',
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
            font = list(size = 21, family = "Arial", color = "black", weight = "bold"), 
            x = 0.3,          
            xanchor = "center",  
            yanchor = "top"
          ),
          xaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
          yaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
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
                    pageSizeOptions = c(25,50,75,100,150,200),
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
        source = 'A',
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
            font = list(size = 21, family = "Arial", color = "black", weight = "bold"), 
            x = 0.3,          
            xanchor = "center",  
            yanchor = "top"
          ),
          xaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
          yaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
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
                    pageSizeOptions = c(25,50,75,100,150,200),
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