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
#' @export

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
                      label = "Select Model",
                      sharedData = shared, 
                      group = ~type),
      filter_select("Lineage", 
                    label = "Select Lineage",
                    sharedData = shared, 
                    group = ~lineage),
      filter_select("Subtype",
                    label = "Select Subtype",
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
                        label = "Select Model",
                        sharedData = shared,
                        group = ~type
        ),
        filter_select("Lineage",
                      label = "Select Lineage",
                      sharedData = shared,
                      group = ~lineage
        ),
        filter_select("Subtype",
                      label = "Select Subtype",
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
  
  
  
  
  
  