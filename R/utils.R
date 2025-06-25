#' 
#' Create a pie chart that illustrates the percentage of lineage or subtype belong to the k nearest neighbors 
#' 
#' @import dplyr
#' @import magrittr
#' @import SNFtool
#' @import stringr
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample a vector of one or multiple TCGA and/or CCLE IDs. 
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param ann annotation file of tumors and cell lines 
#' @param dist_top_n neighbors dataframe
#' @param value string value to select to determine if the pie chart will be based on subtype percentage or lineage percentage
#' @return pie chart illustrating the distribution of k nearest neighbors subtype to the query 
#' @export

get_piechart <- function(combined_mat, input_sample, type, ann, value, dist_top_n) {
  
  annot_col <- ifelse(value == 'lineage', 'lineage', 'subtype_1')
  
  dist_topk_1 <- dist_top_n %>% dplyr::select(sampleID) %>% dplyr::distinct() %>% 
    dplyr::left_join(., ann[, c('sampleID', annot_col)], by = 'sampleID')
  colnames(dist_topk_1) <- c('neighbor', 'annot_neighbor')
  
  dist_topk_2 <- 
    dist_topk_1 %>% dplyr::group_by(annot_neighbor) %>% 
    dplyr::summarize(Freq = n()) %>%
    dplyr::mutate(percentage = Freq / sum(Freq) * 100,
                  label = paste0(round(percentage, 1), "%")) %>%
    dplyr::filter(Freq > 0)
  
    
  dist_topk_2 <- dist_topk_2 %>%
    dplyr::mutate(annot_neighbor = stringr::str_wrap(annot_neighbor, width = 14))
  
  
  ann <- ann %>% dplyr::mutate(annot_neighbor = stringr::str_wrap(ann[[annot_col]], width = 14))
  annot_neighbor_levels <- sort(unique(ann$annot_neighbor))  
  annot_neighbor_colors <- stats::setNames(brewer_recycled("Dark2", length(annot_neighbor_levels)), annot_neighbor_levels)
  
  y <- dist_topk_2 %>% 
    ggplot2::ggplot(., ggplot2::aes(x = "", y = Freq, fill = annot_neighbor)) +
    ggplot2::geom_bar(width = 1, stat = "identity") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::scale_fill_manual(values = annot_neighbor_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5)) +
    ggplot2::geom_text(ggplot2::aes(label = label),
              position = ggplot2::position_stack(vjust = 0.5),
              size = 3) +
    ggplot2::labs(title = paste("Neighbors", value, "distribution")) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = stringr::str_to_title(value)))
  
  return(y)
  
}



#' 
#' Create an interactive plot that highlighting k nearest neighbors for the choosed query
#'
#' @import plotly
#' @importFrom crosstalk filter_checkbox filter_select filter_select bscols
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @import magrittr
#' @import SNFtool
#' @importFrom tibble rownames_to_column
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param ann annotation file of tumors and cell lines 
#' @param input_sample a vector of one or multiple TCGA and/or CCLE IDs. 
#' @param dist_top_n neighbors dataframe
#' @param annot_value indicates the annotation column from the annotation table 'ann'. Values acceptes are: 'lineage' or 'subtype'
#' @return an interactive plot that highlighting the tumor k nearest neighbors 

get_alignment_plot <- function(reduced_mat, ann, input_sample = NULL, dist_top_n = NULL, annot_value = 'lineage') {
  
  annot_col <- ifelse(annot_value == 'lineage', 'lineage', 'subtype_1')
  ann <- ann %>% mutate(str_wrap_annot_col = stringr::str_wrap(ann[[annot_col]], width = 14))
  annot_levels <- sort(unique(ann$str_wrap_annot_col))  
  annot_colors <- stats::setNames(brewer_recycled("Dark2", length(annot_levels)), annot_levels)
  
  data_res <- reduced_mat %>% t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("sampleID") %>% 
    dplyr::left_join(., ann, by = "sampleID")
  colnames(data_res)[2:3] <- c('COMP_1', 'COMP_2')
  
  data_res_1 <- data_res %>% dplyr::mutate(
    size = dplyr::if_else(type == "Tumor", 5, 9),
    width = dplyr::if_else(type == "Tumor", 0.2, 1.3),
    opacity = dplyr::if_else(type == "Tumor", 1, 1),
    show_it = 'not'
  ) %>% dplyr::select(
    "COMP_1", "COMP_2", "stripped_cell_line_name", "sampleID", "lineage",
    "subtype", "subtype_1", "link", 'str_wrap_annot_col', "type", "size", "width", "opacity", "show_it"
  )
  
  if(!is.null(dist_top_n) || !missing(dist_top_n)) {
    print('Generating projection in "Neighbors search" mode...')
    
    top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
    
    data_res_2 <- data_res_1 %>%
      dplyr::mutate(
        'show_it' = dplyr::case_when(
          sampleID %in% top_k_tumors_1 ~ "show",
          sampleID %in% input_sample ~ "input",
          .default = 'not'
        ),
        'size' = dplyr::if_else(show_it %in% c('show', 'input'), 9, 
                                dplyr::if_else(type == "Tumor", 5, 9))
      )
    
    data_res_3 <- data_res_2 %>% dplyr::left_join(dist_top_n, by = 'sampleID')  %>%
      dplyr::arrange(dist) %>% dplyr::mutate(dist = round(dist, 3)) %>% 
      dplyr::select("COMP_1","COMP_2","stripped_cell_line_name","sampleID","lineage",
                    "subtype","subtype_1", "link", 'str_wrap_annot_col', "type",
                    "size", "width", "opacity", "show_it", "dist")
    
  } else {
    data_res_3 <- data_res_1 %>% dplyr::mutate(dist = NA)
  }
  
  shared <- plotly::highlight_key(data_res_3,~stripped_cell_line_name,"Highligth a sample")
  
  if(!is.null(dist_top_n) || !missing(dist_top_n)) {
    react_table <- shared$origData()[shared$origData()$show_it == 'show',]
  } else {
    react_table <- shared$origData()
  }
  
  
  x_range <- range(shared$data()$COMP_1, na.rm = TRUE)
  y_range <- range(shared$data()$COMP_2, na.rm = TRUE)

  alignPlot <- withCallingHandlers({
    
    plot_ly(data = shared,
            x = ~COMP_1,
            y = ~COMP_2,
            source = 'A',
            type = 'scatter',
            mode = 'markers',
            color = ~factor(str_wrap_annot_col, levels = annot_levels),
            colors = annot_colors,
            symbol = ~type, 
            symbols = c('circle',"x"),
            stroke = ~show_it,
            strokes = c('show' = "red", 'input' = "green"),
            size = ~size,
            sizes = c(10, 25),
            hoverinfo = "text",
            hovertext = ~paste("ID:", sampleID,
                               '\nName:', stripped_cell_line_name, 
                               '\nLineage:', lineage,
                               '\nSubtype:', subtype_1,
                               '\nType:', type),
            marker = list(line = list(width = 3)),
            height = 600
    ) %>%
      layout(
        dragmode = "zoom",
        autosize = TRUE,
        xaxis = list(title = "", zeroline = FALSE, showticklabels = FALSE, showgrid = FALSE, range = x_range),
        yaxis = list(title = "", zeroline = FALSE, showticklabels = FALSE, showgrid = FALSE, range = y_range),
        legend = list(title = list(text = "Select lineage-type pair"), traceorder = "normal")
      ) %>% 
      plotly::event_register("plotly_selected") %>% 
      plotly::highlight(on = "plotly_selected", off = "plotly_doubleclick", color = 'green', persistent = FALSE) %>% 
      plotly::highlight(on = "plotly_click", selectize = TRUE, persistent = TRUE, off = "plotly_doubleclick", color = "blue")
    
  }, message = function(m) {
    if (grepl("We recommend setting `persistent` to `FALSE`", conditionMessage(m))) {
      invokeRestart("muffleMessage")  # sopprime solo quel messaggio
    }
  })
  
  row_1 <- bscols(
    widths = c(2, 10),
    list(
      div(style = "height: 40px;"),
      crosstalk::filter_checkbox("Type", 
                                 label = "View Model",
                                 sharedData = shared, 
                                 group = ~type),
      crosstalk::filter_select("Lineage", 
                               label = "View Lineage",
                               sharedData = shared, 
                               group = ~lineage),
      crosstalk::filter_select("Subtype",
                               label = "View Subtype",
                               sharedData = shared, 
                               group = ~subtype)
    ),
    
    div(
      style = "height: 600px; width: 100%;",  
      suppressWarnings(plotly_build(alignPlot))
    )
  )
  
  
  row_2 <- crosstalk::bscols(
    widths = c(12),
    
    htmltools::browsable(
      tagList(
        
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
                  type: row.type ?? '',
                  dist: row.dist ?? ''
                }));
            
                const csvHeader = ['stripped_cell_line_name','sampleID', 'lineage', 'subtype', 'subtype_1', 'type', 'dist'];
                const csvRows = data.map(row => [row.stripped_cell_line_name, row.sampleID, row.lineage, row.subtype, row.subtype_1, row.type, row.dist]);
            
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
        
        reactable::reactable(react_table, 
                             searchable = TRUE, minRows = 3, 
                             showPageSizeOptions = TRUE,
                             pageSizeOptions = c(25, 50, 75, 100, 150, 200, 300, 500),
                             defaultPageSize = 25,
                             resizable = TRUE, #highlight = TRUE, 
                             # selection = "multiple",
                             # onClick = "select",
                             theme = reactable::reactableTheme(
                               headerStyle = list(
                                 "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                 "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"),
                                 borderColor = "#555")),
                             bordered = TRUE,
                             striped = TRUE,
                             filterable = TRUE,
                             elementId = "alignment-download-table",
                             columns = list(
                               stripped_cell_line_name = reactable::colDef(name = 'Name'), 
                               sampleID = reactable::colDef(
                                 cell = function(value, index) {
                                   if (is.na(data_res_3$link[index])) {
                                     value 
                                   } else {
                                     url <- data_res_3$link[index]
                                     htmltools::tags$a(href = url, target = "_blank", value)
                                   }
                                 } 
                               ),
                               type = reactable::colDef(name = "Type"),
                               lineage = reactable::colDef(name = "Lineage"),
                               subtype = reactable::colDef(name = "Subtype"),
                               subtype_1 = reactable::colDef(name = 'Subtype code'),
                               dist = reactable::colDef(name = "Distance"),
                               str_wrap_annot_col = reactable::colDef(show = FALSE),
                               COMP_1 = reactable::colDef(show = FALSE),
                               COMP_2 = reactable::colDef(show = FALSE),
                               link = reactable::colDef(show = FALSE),
                               size = reactable::colDef(show = FALSE),
                               show_it = reactable::colDef(show = FALSE),
                               size = reactable::colDef(show = FALSE),
                               width = reactable::colDef(show = FALSE),
                               opacity = reactable::colDef(show = FALSE)
                             ),
                             
                             style = list(
                               height = "400px",  
                               overflowY = "auto",
                               overflowX = "hidden"
                             )
        )
      )
    )
  )
  
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
  
  if (is.vector(combined_mat)) {
    depmap_id <- combined_mat
  } else {
    depmap_id <- rownames(combined_mat)
  }
  
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
#' @param ann annotation file of tumors and cell lines 
#' @param dist_top_n neighbors dataframe
#' @return lineage and subtype distance distribution plot

c_distribution <- function(dist_top_n, ann) {
  
  df <- dist_top_n %>%
    dplyr::arrange(dist)
  df$position <- 1:nrow(df)
  
  df <- df %>%
    dplyr::left_join(., ann, by = join_by(sampleID))
  df$to_sum <- c(0, df$dist[-length(df$dist)])
  df$cum_dist <- df$dist + df$to_sum
  
  lineage_levels <- sort(unique(ann$lineage))  
  lineage_colors <- stats::setNames(brewer_recycled("Dark2", length(lineage_levels)), lineage_levels)
  
  subtype_levels <- sort(unique(ann$subtype_1))  
  subtype_colors <- stats::setNames(brewer_recycled("Dark2", length(subtype_levels)), subtype_levels)
  
  main_plt <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = dist)) +
    ggplot2::geom_line(linewidth = 1.2, color = "steelblue") +
    ggplot2::geom_point(colour = 'black') +
    ggplot2::theme_minimal() +
    ggplot2::ylab("Distance") +
    ggplot2::xlab(NULL)
  
  lin_plt <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = 1, fill = lineage, text = paste(
    "Position:", position,
    "\nID:", sampleID,
    "\nName:", stripped_cell_line_name,
    "\nLineage:", lineage,
    "\nSubtype:", subtype_1
  ))) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_manual(values = lineage_colors) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = 'Lineage'))
  
  subtype_plt <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = 1, fill = subtype_1, text = paste(
    "Position:", position,
    "\nID:", sampleID,
    "\nName:", stripped_cell_line_name,
    "\nLineage:", lineage,
    "\nSubtype:", subtype_1
  ))) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_manual(values = subtype_colors) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = 'Subtype'))
  
  main_plotly <- ggplotly(main_plt)
  lin_plotly <- ggplotly(lin_plt, tooltip = "text")
  subtype_plotly <- ggplotly(subtype_plt, tooltip = "text")
  
  lin_dist <- subplot(main_plotly, lin_plotly, nrows = 2, heights = c(0.75, 0.25), shareX = TRUE)
  sub_dist <- subplot(main_plotly, subtype_plotly, nrows = 2, heights = c(0.75, 0.25), shareX = TRUE)
  
  return(list(
    lineage_distribution = lin_dist,
    subtype_distribution = sub_dist
  ))
}



#' 
#' Function to create a new color palette of size n starting from an existing RColorBrewer palette.
#' 
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grDevices colorRampPalette
#' @param palette name of a RColorBrewer palette (e.g., "Dark2")
#' @param n numeric integer indicating the number of colors in the new palette
#' @return a new color palette of size equal to n

brewer_recycled <- function(palette = "Dark2", n) {
  max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(max_colors, palette)
  grDevices::colorRampPalette(base_colors)(n)
}
