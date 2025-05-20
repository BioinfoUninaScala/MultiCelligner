#' 
#' Create a pie chart that illustrates the percentage of lineage or subtype belong to the k nearest neighbors 
#' 
#' @import dplyr
#' @import ggplot2
#' @import reshape
#' @import reshape2
#' @import magrittr
#' @import SNFtool
#' @import stringr
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param dist_top_n 
#' @param value string value to select to determine if the pie chart will be based on subtype percentage or lineage percentage
#' @return pie chart illustrating the distribution of k nearest neighbors subtype to the query 
#' @export

get_piechart <- function(combined_mat, input_sample = NULL, selected_samples = NULL , type, k, ann, value, dist_top_n) {
  
  top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
  dist_df <- data.frame(sampleID = top_k_tumors_1)
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','lineage_tcga')
  } else if ("Cell lines" %in% type & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Lineage_CCLE')
  } else if ("Tumors" %in% type & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Lineage_TCGA')
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_tcga')
  } else if ("Cell lines" %in% type & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_CCLE')
  } else if ("Tumors" %in% type & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_TCGA')
  }
  
  if(value %in% "lineage") {
    dist_top25_2 <- dist_top25_1 %>%
      mutate(sampleID = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                            length(sample_1))) %>%
      left_join(., ann[, c(1,2)], by = 'sampleID')
  } else {
    dist_top25_2 <- dist_top25_1 %>%
      mutate(sampleID = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                            length(sample_1))) %>%
      left_join(., ann[, c(1,4)], by = 'sampleID')
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, lineage_tcga, lineage_ccle)
    colnames(dist_top25_3)[3] <- 'Lineage'
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_ccle, Lineage) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage"
  } else if ("Cell lines" %in% type & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_tcga')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Lineage_CCLE, lineage_tcga)
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_tcga, Lineage_CCLE) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage_CCLE"
  } else if ("Tumors" %in% type & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Lineage_TCGA, lineage_ccle)
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_ccle, Lineage_TCGA) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage_TCGA"
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_tcga, subtype_ccle)
    colnames(dist_top25_3)[3] <- 'Subtype'
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_ccle, Subtype) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype"
  } else if ("Cell lines" %in% type & length(type) == 1 & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_tcga')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_CCLE, subtype_tcga)
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_tcga, Subtype_CCLE) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype_CCLE"
  } else if ("Tumors" %in% type & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_TCGA, subtype_ccle)
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_ccle, Subtype_TCGA) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype_TCGA"
  }
  
  dist_top25_4 <- dist_top25_4 %>%
    mutate(percentage = Freq / sum(Freq) * 100,
           label = paste0(round(percentage, 1), "%")) %>%
    dplyr::filter(Freq > 0)
  
  if(all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_wrap = stringr::str_wrap(Subtype, width = 14))
    fill_var <- "Subtype_wrap"
  }
  
  else if ("Cell lines" %in% type & length(type) == 1 & value %in% "subtype") {
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_CCLE_wrap = stringr::str_wrap(Subtype_CCLE, width = 14))
    fill_var <- "Subtype_CCLE_wrap"
  }
  
  y <- ggplot(dist_top25_4, aes(x = "", y = Freq, fill = .data[[fill_var]])) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Spectral") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3) +
    labs(if(value %in% 'lineage') title = "Neighbors lineage distribution" else title = "Subtype lineage distribution")
  
  return(y)
  
}

#' 
#' Create an interactive plot that highlighting k nearest neighbors for the choosed query
#'
#' @import plotly
#' @import crosstalk
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @import magrittr
#' @import reshape
#' @import SNFtool
#' @import reshape2
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param ann annotation file of tumors and cell lines 
#' @param dist_top_n neighbors dataframe
#' @return an interactive plot that highlighting the tumor k nearest neighbors 

get_alignment_plot <- function(reduced_mat, ann, dist_top_n = NULL) {
  
  top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
  
  data_res <- reduced_mat %>% t() %>% 
    as.data.frame() %>% 
    mutate("sampleID" = colnames(reduced_mat)) %>% 
    left_join(., ann, by = "sampleID")
  
  colnames(data_res)[c(1,2)] <- c("UMAP_1", "UMAP_2")
  
  if(!is.null(dist_top_n) || !missing(dist_top_n)) {
    
    data_res_1 <- data_res %>% mutate('show_it' = ifelse(data_res$sampleID %in% top_k_tumors_1, 'show', 'not'))
    
    data_res_2 <- data_res_1 %>% mutate('size' = if_else(show_it == 'show', 16, 5))
    
    
    
    dist_metasample_2 <- dist_top_n %>% dplyr::select("sampleID", "dist")
    data_res_3 <- data_res_2 %>% left_join(dist_metasample_2, by = 'sampleID')  %>% arrange(dist)
    data_res_3$dist <- round(data_res_3$dist, 3)
    
    data_res_3 <- data_res_3 %>% select("UMAP_1","UMAP_2","stripped_cell_line_name","sampleID","lineage",
                                        "subtype","subtype_1", "link", "type","dist","show_it","size")
    
    shared <- SharedData$new(data_res_3, key = ~sampleID)
    
    row_1 <- bscols(
      widths = c(2, 10),
      list(
        div(style = "height: 40px;"),
        filter_checkbox("Type", 
                        label = "View Model",
                        sharedData = shared, 
                        group = ~type),
        filter_select("Lineage", 
                      label = "View Lineage",
                      sharedData = shared, 
                      group = ~lineage),
        filter_select("Subtype",
                      label = "View Subtype",
                      sharedData = shared, 
                      group = ~subtype)
      ),
      
      div(
        style = "height: 600px; width: 100%;",  
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
            dragmode = "zoom",
            xaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
            yaxis = list(zeroline  = F, showticklabels = FALSE, showgrid = FALSE, title = ''),
            legend = list(
              title = list(text = 'Select lineage-type pair'),
              traceorder = 'normal'),
            height = 600) %>% 
          event_register("plotly_selected") %>% 
          highlight(on = "plotly_selected", off = "plotly_doubleclick",color = 'green', persistent = FALSE)
      ))
    
    
    row_2 <- crosstalk::bscols(
      widths = c(10,2),
      
      htmltools::browsable(
        tagList(
          tags$button(
            tagList(fontawesome::fa("download"), "Download Table"),
            id = "download-table-btn"
          ),
          tags$script(HTML("
  document.getElementById('download-table-btn').onclick = function() {
    const table = Reactable.getState('alignment-download-table');
    if (!table) {
      alert('Table not ready yet!');
      return;
    }

    const data = table.data.map(row => ({
      stripped_cell_line_name: row.stripped_cell_line_name ?? '',
      sampleID: row.sampleID ?? '',
      lineage: row.lineage ?? '',
      subtype: row.subtype ?? '',
      subtype_1: row.subtype_1 ?? '',
      type: row.type ?? ''
    }));

    const csvHeader = ['stripped_cell_line_name','sampleID', 'lineage', 'subtype', 'subtype_1', 'type'];
    const csvRows = data.map(row => [row.stripped_cell_line_name, row.sampleID, row.lineage, row.subtype, row.subtype_1, row.type]);

    const csvContent = [csvHeader, ...csvRows]
      .map(e => e.join(','))
      .join('\\n');

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.setAttribute('download', 'alignment.csv');
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };
")),
          
          tags$div(
            style = "margin-bottom: 10px;",
            tags$button(
              tagList(fontawesome::fa("download"), "Download plot as SVG"), 
              id = "download-svg-btn")
          ),
          tags$script(HTML(
            "
  document.getElementById('download-svg-btn').onclick = function() {
    var plot = document.getElementsByClassName('plotly')[0];
    Plotly.downloadImage(plot, {
      format: 'svg',
      filename: 'alignment_plot',
      width: 1500,
      height: 700,
      scale: 1
    });
  };
  "
          )),
          
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
                      sampleID = colDef(
                        cell = function(value, index) {
                          if (is.na(data_res_3$link[index])) {
                            value 
                          } else {
                            url <- data_res_3$link[index]
                            htmltools::tags$a(href = url, target = "_blank", value)
                          }
                        } 
                      ),
                      type = colDef(name = "Type"),
                      lineage = colDef(name = "Lineage"),
                      subtype = colDef(name = "Subtype"),
                      subtype_1 = colDef(name = 'Subtype code'),
                      dist = colDef(name = "Distance"),
                      UMAP_1 = colDef(show = FALSE),
                      UMAP_2 = colDef(show = FALSE),
                      link = colDef(show = FALSE),
                      size = colDef(show = FALSE),
                      show_it = colDef(show = FALSE)),
                    
                    style = list(
                      height = "400px",  
                      overflowY = "auto",
                      overflowX = "hidden"
                    )
          ))))
    
    
  } else {
    
    data_res_1 <- data_res %>% mutate(
      size = if_else(data_res$type == "tumor", 5, 9),
      width = if_else(data_res$type == "tumor", 0.2, 1.3),
      opacity = if_else(data_res$type == "tumor", 1, 1)
    )
    
    data_res_1 <- data_res_1 %>% select(
      "UMAP_1", "UMAP_2", "stripped_cell_line_name", "sampleID", "lineage",
      "subtype", "subtype_1", "link", "type", "size", "width", "opacity"
    )
    
    shared <- SharedData$new(data_res_1, key = ~sampleID)
    
    row_1 <- bscols(
      widths = c(2, 10),
      list(
        div(style = "height: 40px;"),
        filter_checkbox("Type",
                        label = "View Model",
                        sharedData = shared,
                        group = ~type
        ),
        filter_select("Lineage",
                      label = "View Lineage",
                      sharedData = shared,
                      group = ~lineage
        ),
        filter_select("Subtype",
                      label = "View Subtype",
                      sharedData = shared,
                      group = ~subtype
        )
      ),
      div(
        style = "height: 600px; width: 100%;",
        plot_ly(
          data = shared,
          x = ~UMAP_1,
          y = ~UMAP_2,
          key = ~sampleID,
          source = "A",
          type = "scatter",
          mode = "markers",
          color = ~lineage,
          symbol = ~type,
          symbols = c("circle", "x"),
          stroke = ~type,
          strokes = c("CL" = "black"),
          size = ~size,
          sizes = c(5, 9),
          hoverinfo = "text",
          hovertext = ~ paste(
            "ID:", sampleID,
            "\nName:", stripped_cell_line_name,
            "\nLineage:", lineage,
            "\nSubtype:", subtype_1,
            "\nType:", type
          ),
          marker = list(
            line = list(width = 1.3)
          )
        ) %>%
          layout(
            dragmode = "zoom",
            autosize = TRUE,
            xaxis = list(title = "", zeroline = F, showticklabels = FALSE, showgrid = FALSE),
            yaxis = list(title = "", zeroline = F, showticklabels = FALSE, showgrid = FALSE),
            legend = list(
              title = list(text = "Select lineage-type pair"),
              traceorder = "normal"
            ),
            height = 600
          ) %>%
          event_register(event = "plotly_selected") %>%
          highlight(on = "plotly_selected", off = "plotly_doubleclick", color = "green", persistent = FALSE)
      )
    )
    
    
    row_2 <- crosstalk::bscols(
      widths = c(10, 2),
      htmltools::browsable(
        tagList(
          tags$button(
            tagList(fontawesome::fa("download"), "Download Table"),
            id = "download-table-btn"
          ),
          tags$script(HTML("
  document.getElementById('download-table-btn').onclick = function() {
    const table = Reactable.getState('alignment-download-table');
    if (!table) {
      alert('Table not ready yet!');
      return;
    }

    const data = table.data.map(row => ({
      stripped_cell_line_name: row.stripped_cell_line_name ?? '',
      sampleID: row.sampleID ?? '',
      lineage: row.lineage ?? '',
      subtype: row.subtype ?? '',
      subtype_1: row.subtype_1 ?? '',
      type: row.type ?? ''
    }));

    const csvHeader = ['stripped_cell_line_name','sampleID', 'lineage', 'subtype', 'subtype_1', 'type'];
    const csvRows = data.map(row => [row.stripped_cell_line_name, row.sampleID, row.lineage, row.subtype, row.subtype_1, row.type]);

    const csvContent = [csvHeader, ...csvRows]
      .map(e => e.join(','))
      .join('\\n');

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.setAttribute('download', 'alignment.csv');
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };
")),
          
          tags$div(
            style = "margin-bottom: 10px;",
            tags$button(
              tagList(fontawesome::fa("download"), "Download plot as SVG"),
              id = "download-svg-btn"
            )
          ),
          tags$script(HTML(
            "
                document.getElementById('download-svg-btn').onclick = function() {
                  var plot = document.getElementsByClassName('plotly')[0];
                  Plotly.downloadImage(plot, {
                    format: 'svg',
                    filename: 'alignment_plot',
                    width: 1500,
                    height: 700,
                    scale: 1
                  });
                };
            "
          )),
          reactable(
            shared,
            searchable = TRUE,
            minRows = 3,
            showPageSizeOptions = TRUE,
            pageSizeOptions = c(25, 50, 75, 100, 150, 200, 300, 500),
            defaultPageSize = 25,
            resizable = TRUE, #highlight = TRUE,
            #selection = "multiple",
            #onClick = "select",
            theme = reactableTheme(
              headerStyle = list(
                "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"),
                borderColor = "#555"
              )
            ),
            bordered = TRUE,
            striped = TRUE,
            filterable = TRUE,
            elementId = "alignment-download-table",
            columns = list(
              stripped_cell_line_name = colDef(name = "Name"),
              sampleID = colDef(
                cell = function(value, index) {
                  if (is.na(data_res_1$link[index])) {
                    value
                  } else {
                    url <- data_res_1$link[index]
                    htmltools::tags$a(href = url, target = "_blank", value)
                  }
                }
              ),
              type = colDef(name = "Type"),
              lineage = colDef(name = "Lineage"),
              subtype = colDef(name = "Subtype"),
              subtype_1 = colDef(name = "Subtype code"),
              link = colDef(show = FALSE),
              UMAP_1 = colDef(show = FALSE),
              UMAP_2 = colDef(show = FALSE),
              size = colDef(show = FALSE),
              width = colDef(show = FALSE),
              opacity = colDef(show = FALSE)
            ),
            style = list(
              height = "400px",
              overflowY = "auto",
              overflowX = "hidden"
            )))))
  }
  
  x <- htmltools::browsable(
    htmltools::tagList(row_1, row_2))
  
  return(x)
  
}



#' 
#' Get the stripped cell line name of the user selected combined matrix
#' 
#' @import dplyr
#' @param combined_mat matrix samples x genes of corrected data by MNN: 
#' @param ann annotation file of tumors and cell lines 
#' @return a vector with the stripped cell line name


get_CL_strp_names <- function(combined_mat ,ann) {
  
  depmap_id <- rownames(combined_mat)
  depmap_id_df <- as.data.frame(depmap_id)
  colnames(depmap_id_df)[1] <- 'sampleID'
  df <- left_join(depmap_id_df, ann[,c(1,6)], by = 'sampleID')
  nm <- df$stripped_cell_line_name[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', df$stripped_cell_line_name)]
  
  return(nm)
  
}

#' 
#' Generate a lineage and subtype distance distribution plot
#'
#' @import plotly
#' @import dplyr
#' @import scales
#' @param ann annotation file of tumors and cell lines 
#' @param dist_top_n neighbors dataframe
#' @return lineage and subtype distance distribution plot

c_distribution <- function(dist_top_n, ann) {
  
  df <- dist_top_n %>%
    arrange(dist)  
  
  df$position <- 1:nrow(df)
  
  df <- df %>%
    left_join(., ann, by = join_by(sampleID))
  
  df$to_sum <- c(0, df$dist[-length(df$dist)])
  df$cum_dist <- df$dist + df$to_sum
  
  p1 <- ggplot(df, aes(x = position, y = dist)) +
    geom_line(linewidth = 1.2, color = "steelblue") +
    theme_minimal() +
    ylab("Cumulative Distance") +
    xlab(NULL) 
  
  p2 <- ggplot(df, aes(x = position, y = 1, fill = lineage)) +
    geom_tile() +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    scale_fill_brewer(palette = "Dark2")
  #theme(legend.position = "bottom") 
  
  
  # p2 <- ggplot(df, aes(x = position, y = 1, fill = lineage, text= paste(
  #   "Position:", position,
  #   "\nID:", sampleID,
  #   "\nName:", stripped_cell_line_name,
  #   "\nLineage:", lineage,
  #   "\nSubtype:", subtype_1,))) +
  #   geom_tile() +
  #   scale_y_continuous(expand = c(0, 0)) +
  #   theme_void() +
  #   scale_fill_brewer(palette = "Dark2")
  # #theme(legend.position = "bottom") 
  
  
  
  p3 <- ggplot(df, aes(x = position, y = 1, fill = subtype_1)) +
    geom_tile() +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Dark2")
  
  p1_plotly <- ggplotly(p1)
  p2_plotly <- ggplotly(p2, tooltip = "text")
  p3_plotly <- ggplotly(p3)
  
  x <- subplot(p1_plotly, p3_plotly, nrows = 2, heights = c(0.75, 0.25), shareX = TRUE)
  y <- subplot(p1_plotly, p2_plotly, nrows = 2, heights = c(0.75, 0.25), shareX = TRUE)
  
  return(list(lineage_distribution = x,
              subtype_distribution = y))
  
}
