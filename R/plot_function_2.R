#'
#' Generate an interactive plot from the dimensionally reduced matrix
#'
#' @import plotly
#' @import crosstalk
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @param reduced_mat dimensionally reduced matrix (UMAP): sample x features
#' @param ann annotation file of tumors and cell lines
#' @param omics_name name of the shiny selected omics
#' @return a interactive plot
#' @export


my_plotting <- function(reduced_mat, ann, omics_name) {
  data_res <- reduced_mat %>%
    t() %>%
    as.data.frame() %>%
    mutate("sampleID" = colnames(reduced_mat)) %>%
    left_join(., ann, by = "sampleID")

  colnames(data_res)[c(1, 2)] <- c("UMAP_1", "UMAP_2")

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
          onclick = "Reactable.downloadDataCSV('alignment-download-table', 'alignment.csv')"
        ),
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
          resizable = TRUE, highlight = TRUE,
          selection = "multiple",
          onClick = "select",
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
          )
        )
      )
    )
  )


  x <- htmltools::browsable(
    htmltools::tagList(row_1, row_2)
  )


  return(x)
}
