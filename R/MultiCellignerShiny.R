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

MultiCellignerShiny <- function() {shiny::shinyApp(ui, server)}

ui <- fluidPage(
  
  titlePanel("MultiCelligner"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        type = "pills",
        
        tabPanel("neighbors", 
                 
                 hr(),
                 
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
                          actionButton("get_plot", "show"))
                 ),
                 
                 hr(),
                 
                 fluidRow(
                   column(6,
                          checkboxGroupInput("df_selection_input", "query sample type:", 
                                             choices = c("CL", "tumor"), selected = 'CL')),
                   column(6,
                          checkboxGroupInput("df_selection_output", "neighbors type:", 
                                             choices = c("CL", "tumor"), selected = 'tumor'))
                 ),
                 
                 hr(), 
                 
                 fluidRow(
                   column(6, 
                          selectizeInput("both_sample", "get neighbors: select the sample/s", 
                                         choices = NULL, multiple = TRUE)),
                   column(3, 
                          numericInput("num_neighbors", "number of neighbors:", 
                                       value = 25, min = 1)),
                   column(3, 
                          actionButton("subset_btn", "show"))
                 ),
                 
                 plotOutput("piechart")),
        
        tabPanel("heatmap", 
                 h3("Quality score of the aligment:", style = "font-size: 14px; text-align: center;"),
                 hr(),
                 plotOutput("heatmap", height = "400px", width = "430px")))),
    
    mainPanel(
      uiOutput("plot"))))

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
           "Expression" = t(index_exp@data),
           "Methylation" = t(index_meth@data),
           "Mutational_process" = t(index_mut@data),
           "MoNETA_multiomics_tSNE" = t(index_MoNETA@data),
           "MoNETA_multiomics_UMAP" = t(index_MoNETA@data),
           "MOFA_multiomics_tSNE" = t(index_MOFA@data),
           "MOFA_multiomics_UMAP" = t(index_MOFA@data),
           "SNF_multiomics_UMAP" = t(index_SNF@data)
    )
  })
  
  selected_BNindex <- reactive({
    switch(input$omics_plot,
           "Expression" = index_exp,
           "Methylation" = index_meth,
           "Mutational_process" = index_mut,
           "MoNETA_multiomics_tSNE" = index_MoNETA,
           "MoNETA_multiomics_UMAP" = index_MoNETA,
           "MOFA_multiomics_tSNE" = index_MOFA,
           "MOFA_multiomics_UMAP" = index_MOFA,
           "SNF_multiomics_UMAP" = index_SNF
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
  })
  
  selected_samples <- reactive({
    input$both_sample
  })
  
  observeEvent(input$subset_btn, {
    
    if (is.null(input$both_sample) || length(input$both_sample) == 0) {
      showNotification("Please select at least one sample before proceeding.", type = "error", duration = 3)
      return()
    }
    
    if(length(input$both_sample) == 1) {
    
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
    
    } else {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      if("CL" %in% input$df_selection_output) {
        
        x <- omics_metasample_both_CL(combined_mat = selected_combined_mat(), 
                                      reduced_mat = selected_reduced_mat(), 
                                      selected_samples = selected_samples(), 
                                      n = input$num_neighbors, 
                                      ann = ann_multiomics_v6)
      }
      
      if("tumor" %in% input$df_selection_output) {
        
        x <- omics_metasample_both_tumor(combined_mat = selected_combined_mat(), 
                                         reduced_mat = selected_reduced_mat(), 
                                         selected_samples = selected_samples(), 
                                         n = input$num_neighbors, 
                                         ann = ann_multiomics_v6)
      }
      
      
      if(all(c("CL", "tumor") %in% input$df_selection_output)) {
        
        x <- omics_metasample_both_both(combined_mat = selected_combined_mat(), 
                                        reduced_mat = selected_reduced_mat(), 
                                        selected_samples = selected_samples(), 
                                        n = input$num_neighbors, 
                                        ann = ann_multiomics_v6)
      }
      
      output$plot <- renderUI({
        x  
      })
    }
  })
  
  observeEvent(input$subset_btn, {
    
    if(length(input$both_sample) == 1) {
    
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
    
    } else {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      if (all(c("CL", "tumor") %in% input$df_selection_output)) {
      
      piechart_2 <- get_piechart_both_2(combined_mat = selected_combined_mat(),
                                        selected_samples = selected_samples(),
                                        n = input$num_neighbors,
                                        ann = ann_multiomics_v6)
    }
    
    
    if ("CL" %in% input$df_selection_output & !'tumor' %in% input$df_selection_output) {
      
      piechart_2 <- get_piechart_CL_2(combined_mat = selected_combined_mat(),
                                      selected_samples = selected_samples(),
                                      n = input$num_neighbors,
                                      ann = ann_multiomics_v6)
      
    }
    
    if ("tumor" %in% input$df_selection_output & !"CL" %in% input$df_selection_output) {
      
      piechart_2 <- get_piechart_tumor_2(combined_mat = selected_combined_mat(),
                                         selected_samples = selected_samples(),
                                         n = input$num_neighbors,
                                         ann = ann_multiomics_v6)
    }
    
    output$piechart <- renderPlot({
      piechart_2
    })
  }
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
