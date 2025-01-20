#' 
#' MultiCelligner Shiny App
#'
#' @import dplyr
#' @import plotly
#' @import shiny
#' @import tidyverse
#' @import crosstalk
#' @import htmltools
#' @import fontawesome
#' @return Shiny app
#' @export

MultiCellignerShiny = shiny::shinyApp(ui, server)

ui <- fluidPage(
  
  titlePanel("MultiCelligner"),
  
  sidebarLayout( 
    sidebarPanel(
      fluidRow(
        column(9,
               selectInput("omics_plot", "Select the omics plot:",
                           choices = c("Expression", 
                                       "Methylation",
                                       "Mutational_process",
                                       "MoNETA_multiomics_tSNE",
                                       "MoNETA_multiomics_UMAP",
                                       "MOFA_multiomics_tSNE",
                                       "MOFA_multiomics_UMAP",
                                       "SNF_multiomics_UMAP"))),
        column(3,
               actionButton("get_plot", "show"))),
      
      hr(),
      
      fluidRow(
        column(6,
               checkboxGroupInput(
                 "df_selection_input", 
                 "query sample type:", 
                 choices = c("CL", "tumor"), 
                 selected = NULL)),
        column(6,
               checkboxGroupInput(
                 "df_selection_output", 
                 "neighbors type:", 
                 choices = c("CL", "tumor"), 
                 selected = NULL))),
      hr(), 
      
      fluidRow(
        column(6, 
               selectizeInput("both_sample", 
                              "get neighbors: write or select the sample", 
                              choices = NULL, 
                              multiple = FALSE)),
        column(3, 
               numericInput("num_neighbors", 
                            "number of neighbors:", 
                            value = 25, min = 1)),
        column(3, 
               actionButton("subset_btn", "show"))),
      
      
      hr(),
      
      fluidRow(
        column(6,
               selectizeInput("multiple_samples", 
                              "find neighbors for multiple_samples", 
                              choices = NULL, 
                              multiple = TRUE)),
        column(3,
               numericInput("num_neighbors_multiple", 
                            "number of neighbors:", 
                            value = 25, min = 1)),
        column(3,
               actionButton("get_centroid", "show"))),
      
      div(style = "max-height: 500px; overflow-y: scroll; overflow-x: hidden;",
          tableOutput("top_k_tumors")),
      hr(),
      plotOutput("piechart"),
      
      hr(),
      
      div(style = "max-height: 500px; overflow-y: scroll; overflow-x: hidden;",
          tableOutput("dist_top_n_1")),
      hr(),
      plotOutput("piechart_2")),
    
    mainPanel(
      
      uiOutput("plot"),
      plotOutput("heatmap", height = "600px", width = "800px") 
    )))

server <- function(input, output, session) {
  
  selected_plot <- reactive({
    switch(input$omics_plot,
           "Expression" = my_plotting(exp_umap, ann_multiomics_v6),
           "Methylation" = my_plotting(meth_umap, ann_multiomics_v6),
           "Mutational_process" = my_plotting(mut_umap, ann_multiomics_v6),
           "MoNETA_multiomics_tSNE" = my_plotting(tsne_MoNETA, ann_multiomics_v6),
           "MoNETA_multiomics_UMAP" = my_plotting(umap_MoNETA, ann_multiomics_v6),
           "MOFA_multiomics_tSNE" = my_plotting(tsne_MOFA, ann_multiomics_v6),
           "MOFA_multiomics_UMAP" = my_plotting(umap_MOFA, ann_multiomics_v6),
           "SNF_multiomics_UMAP" = my_plotting(SNF_umap, ann_multiomics_v6))
  })
  
  output$plot <- renderUI({
    selected_plot()  
  })
  
  observeEvent(input$get_plot, {
    output$plot <- renderUI({
      selected_plot()  
    })
  })
  
  selected_combined_mat <- reactive({
    switch(input$omics_plot,
           "Expression" = t(exp_index@data),
           "Methylation" = t(meth_index@data),
           "Mutational_process" = t(mut_index@data),
           "MoNETA_multiomics_tSNE" = t(MoNETA_index@data),
           "MoNETA_multiomics_UMAP" = t(MoNETA_index@data),
           "MOFA_multiomics_tSNE" = t(MOFA_index@data),
           "MOFA_multiomics_UMAP" = t(MOFA_index@data),
           "SNF_multiomics_UMAP" = t(SNF_index@data)
    )
  })
  
  selected_BNindex <- reactive({
    switch(input$omics_plot,
           "Expression" = exp_index,
           "Methylation" = meth_index,
           "Mutational_process" = mut_index,
           "MoNETA_multiomics_tSNE" = MoNETA_index,
           "MoNETA_multiomics_UMAP" = MoNETA_index,
           "MOFA_multiomics_tSNE" = MOFA_index,
           "MOFA_multiomics_UMAP" = MOFA_index,
           "SNF_multiomics_UMAP" = SNF_index
    )
  })
  
  selected_reduced_mat <- reactive({
    switch(input$omics_plot,
           "Expression" = exp_umap, 
           "Methylation" = meth_umap,
           "Mutational_process" = mut_umap,
           "MoNETA_multiomics_tSNE" = tsne_MoNETA,
           "MoNETA_multiomics_UMAP" = umap_MoNETA,
           "MOFA_multiomics_tSNE" = tsne_MOFA,
           "MOFA_multiomics_UMAP" = umap_MOFA,
           "SNF_multiomics_UMAP" = SNF_umap
    )
  })
  
  selected_sample_order <- reactive({
    switch(input$omics_plot,
           "Expression" = sample_exp_order, 
           "Methylation" = sample_meth_order,
           "Mutational_process" = sample_mut_order,
           "MoNETA_multiomics_tSNE" = sample_MoNETA_order,
           "MoNETA_multiomics_UMAP" = sample_MoNETA_order,
           "MOFA_multiomics_tSNE" = sample_MOFA_order,
           "MOFA_multiomics_UMAP" = sample_MOFA_order,
           "SNF_multiomics_UMAP" = sample_SNF_order
    )
  })
  
  observe({
    choices <- list() 
    
    if ("CL" %in% input$df_selection_input) {
      cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v6)
      cl_values <- rownames(selected_combined_mat())[!grepl('TCGA', rownames(selected_combined_mat()))] 
      choices <- c(choices, setNames(as.list(cl_values), cl_names))
    }
    
    if ("tumor" %in% input$df_selection_input) {
      tumor_names <- rownames(selected_combined_mat())[grepl('TCGA', rownames(selected_combined_mat()))]
      tumor_values <- rownames(selected_combined_mat())[grepl('TCGA', rownames(selected_combined_mat()))]
      choices <- c(choices, setNames(as.list(tumor_values), tumor_names))
    }
    
    if (all(c("CL", "tumor") %in% input$df_selection_input)) {
      all_names <- get_all_name_strpp(combined_mat = selected_combined_mat(), 
                                      ann = ann_multiomics_v6, 
                                      sample_info = sample_info)
      all_values <- rownames(selected_combined_mat()) 
      choices <- c(choices, setNames(as.list(all_values), all_names))
    }
    
    updateSelectizeInput(session, "both_sample", choices = choices)
    updateSelectizeInput(session, "multiple_samples", choices = choices)
  })
  
  selected_samples <- reactive({
    input$multiple_samples
  })
  
  observeEvent(input$subset_btn, {
    
    if ("CL" %in% input$df_selection_output) {
      top_k_tumors <- get_neighbors_table_CL(combined_mat = selected_combined_mat(),
                                             input_sample = input$both_sample,
                                             k = input$num_neighbors,
                                             ann = ann_multiomics_v6,
                                             BNindex = selected_BNindex(),
                                             sample_order = selected_sample_order())
    }
    
    if ("tumor" %in% input$df_selection_output) {
      
      top_k_tumors <- get_neighbors_table_tumor(combined_mat = selected_combined_mat(),
                                                input_sample = input$both_sample,
                                                k = input$num_neighbors,
                                                ann = ann_multiomics_v6,
                                                BNindex = selected_BNindex(),
                                                sample_order = selected_sample_order())
      
    }
    
    
    if(all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      top_k_tumors <- get_neighbors_table_both(combined_mat = selected_combined_mat(),
                                               input_sample = input$both_sample,
                                               k = input$num_neighbors,
                                               ann = ann_multiomics_v6,
                                               BNindex = selected_BNindex(),
                                               sample_order = selected_sample_order())
    }
    
    
    output$top_k_tumors <- renderTable({
      top_k_tumors
    })
  })
  
  observeEvent(input$subset_btn, {
    
    if ("CL" %in% input$df_selection_output) {
      
      x <- find_neighbors_CL(combined_mat = selected_combined_mat(), 
                             reduced_mat = selected_reduced_mat(),
                             input_sample = input$both_sample,
                             k = input$num_neighbors,
                             ann = ann_multiomics_v6,
                             BNindex = selected_BNindex(), 
                             sample_order = selected_sample_order()) 
    }
    
    if ("tumor" %in% input$df_selection_output) {
      x <- find_neighbors_tumor(combined_mat = selected_combined_mat(), 
                                reduced_mat = selected_reduced_mat(),
                                input_sample = input$both_sample,
                                k = input$num_neighbors,
                                ann = ann_multiomics_v6,
                                BNindex = selected_BNindex(), 
                                sample_order = selected_sample_order())
    }
    
    if (all(c("CL", "tumor") %in% input$df_selection_output)) {
      x <- find_neighbors_both(combined_mat = selected_combined_mat(), 
                               reduced_mat = selected_reduced_mat(),
                               input_sample = input$both_sample,
                               k = input$num_neighbors,
                               ann = ann_multiomics_v6,
                               BNindex = selected_BNindex(),
                               sample_order = selected_sample_order())
    }
    
    output$plot <- renderUI({
      x  
    })
  }) 
  
  observeEvent(input$get_centroid, {
    
    selected_samples <- reactive({
      input$multiple_samples
    })
    
    if("CL" %in% input$df_selection_output) {
      
      x <- omics_metasample_both_CL(combined_mat = selected_combined_mat(), 
                                    reduced_mat = selected_reduced_mat(), 
                                    selected_samples = selected_samples(), 
                                    n = input$num_neighbors_multiple, 
                                    ann = ann_multiomics_v6)
    }
    
    if("tumor" %in% input$df_selection_output) {
      
      x <- omics_metasample_both_tumor(combined_mat = selected_combined_mat(), 
                                       reduced_mat = selected_reduced_mat(), 
                                       selected_samples = selected_samples(), 
                                       n = input$num_neighbors_multiple, 
                                       ann = ann_multiomics_v6)
    }
    
    
    if(all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      x <- omics_metasample_both_both(combined_mat = selected_combined_mat(), 
                                      reduced_mat = selected_reduced_mat(), 
                                      selected_samples = selected_samples(), 
                                      n = input$num_neighbors_multiple, 
                                      ann = ann_multiomics_v6)
    }
    
    output$plot <- renderUI({
      x  
    })
  }) 
  
  observeEvent(input$get_centroid, {
    
    
    if("CL" %in% input$df_selection_output) {
      
      x_2 <- get_metasample_ann_both_CL(combined_mat = selected_combined_mat(),
                                        reduced_mat = selected_reduced_mat(),
                                        selected_samples = selected_samples(),
                                        n = input$num_neighbors_multiple,
                                        ann = ann_multiomics_v6)
    }
    
    if("tumor" %in% input$df_selection_output) {
      
      x_2 <- get_metasample_ann_both_tumor(combined_mat = selected_combined_mat(),
                                           reduced_mat = selected_reduced_mat(),
                                           selected_samples = selected_samples(),
                                           n = input$num_neighbors_multiple,
                                           ann = ann_multiomics_v6)
    }
    
    if(all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      x_2 <-  get_metasample_ann_both_both(combined_mat = selected_combined_mat(),
                                           reduced_mat = selected_reduced_mat(),
                                           selected_samples = selected_samples(),
                                           n = input$num_neighbors_multiple,
                                           ann = ann_multiomics_v6)
    }
    
    
    output$dist_top_n_1 <- renderTable({
      x_2
    })
  })
  
  observeEvent(input$subset_btn, {
    
    if ("CL" %in% input$df_selection_output) {
      
      piechart <- get_piechart_CL(combined_mat = selected_combined_mat(),
                                  input_sample = input$both_sample,
                                  k = input$num_neighbors,
                                  ann = ann_multiomics_v6,
                                  BNindex = selected_BNindex(),
                                  sample_order = selected_sample_order())
      
    }
    
    if ("tumor" %in% input$df_selection_output) {
      
      piechart <- get_piechart_tumor(combined_mat = selected_combined_mat(),
                                     input_sample = input$both_sample,
                                     k = input$num_neighbors,
                                     ann = ann_multiomics_v6,
                                     BNindex = selected_BNindex(),
                                     sample_order = selected_sample_order())
    }
    
    if (all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      piechart <- get_piechart_both(combined_mat = selected_combined_mat(),
                                    input_sample = input$both_sample,
                                    k = input$num_neighbors,
                                    ann = ann_multiomics_v6,
                                    BNindex = selected_BNindex(),
                                    sample_order = selected_sample_order())
    }
    
    output$piechart <- renderPlot({
      piechart
    })
  })
  
  observeEvent(input$get_centroid, {
    
    if (all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      piechart_2 <- get_piechart_both_2(combined_mat = selected_combined_mat(),
                                        selected_samples = selected_samples(),
                                        n = input$num_neighbors_multiple,
                                        ann = ann_multiomics_v6)
    }
    
    
    if ("CL" %in% input$df_selection_output & !'tumor' %in% input$df_selection_output) {
      
      piechart_2 <- get_piechart_CL_2(combined_mat = selected_combined_mat(),
                                      selected_samples = selected_samples(),
                                      n = input$num_neighbors_multiple,
                                      ann = ann_multiomics_v6)
      
    }
    
    if ("tumor" %in% input$df_selection_output & !"CL" %in% input$df_selection_output) {
      
      piechart_2 <- get_piechart_tumor_2(combined_mat = selected_combined_mat(),
                                         selected_samples = selected_samples(),
                                         n = input$num_neighbors_multiple,
                                         ann = ann_multiomics_v6)
    }
    
    output$piechart_2 <- renderPlot({
      piechart_2
    })
  })
  
  selected_heatmap <- reactive({
    switch(input$omics_plot,
           "Expression" = get_confusion_matrix_shiny(exp_m3),
           "Methylation" = get_confusion_matrix_shiny(meth_m3),
           "Mutational_process" = get_confusion_matrix_shiny(mut_m3),
           "MoNETA_multiomics_tSNE" = get_confusion_matrix_shiny(MoNETA_m3),
           "MoNETA_multiomics_UMAP" = get_confusion_matrix_shiny(MoNETA_m3),
           "MOFA_multiomics_tSNE" = get_confusion_matrix_shiny(MOFA_m3),
           "MOFA_multiomics_UMAP" = get_confusion_matrix_shiny(MOFA_m3),
           "SNF_multiomics_UMAP" = get_confusion_matrix_shiny(SNF_m3))
  })
  
  output$heatmap <- renderPlot({
    heatmap_plot <- selected_heatmap()
    heatmap_plot
  })
  
}
