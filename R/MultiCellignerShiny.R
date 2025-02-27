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
#' @import shinyjs
#' @return Shiny app
#' @export

MultiCellignerShiny <- function() {shiny::shinyApp(ui, server)}

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      tags$style(HTML("hr { margin-top: 4px !important; margin-bottom: 4px !important; padding: 4px !important; }")),
      
      div(h3("MultiCelligner"), class = "text-center"),
      hr(),
      
      
      selectInput("omics_plot", "Select an omics alignment:",
                  choices = c("Methylation",
                              "Mutational signature (COSMIC)", 
                              "Expression",
                              "MoNETA multiomics tSNE",  
                              "MoNETA multiomics UMAP", 
                              "MOFA multiomics tSNE",
                              "MOFA multiomics UMAP",
                              "SNF multiomics UMAP")),
      
      hr(),
      
      tags$p("Find neighbors", style = "text-align: center; font-size: 18px; font-weight: bold;"),
      
      fluidRow(
        column(6,
               selectizeInput("sel_type",
                              'Select type',
                              choices = c('Cell lines', 'Tumors'), 
                              multiple = TRUE)),
        column(6,
               selectizeInput('sel_lineage',
                              'Select lineage',
                              choices = NULL,
                              multiple = FALSE))
      ),
      
      
      selectizeInput("both_sample", 
                     "Select one or more reference tumor/cell line:", 
                     choices = NULL, 
                     multiple = TRUE),
      
      hr(),
      
      fluidRow(
        column(8,
               tags$strong("or load from map selection:")),
        column(4,
               actionButton("load_selection", "Load", style = "text-align: center;") %>%
                 tagAppendChild(tags$script(HTML('
    $(document).ready(function(){
      $("#load_selection").tooltip({
        title: "Possibility to drop in the samples directly from the plot using Lasso Select",
        placement: "right",
        trigger: "hover",
        html: true
      });
    });
  ')))
        )),
      
      
      hr(),
      
      fluidRow(
        column(7,
               checkboxGroupInput("df_selection_output", 
                                  "Search among the closest:", 
                                  choices = c("Cell lines", "Tumors"), 
                                  selected = 'Tumors', inline = TRUE)),
        column(5,
               div("Number of neighbors:", style = "text-align: left; font-weight: bold;"),
               numericInput("num_neighbors", NULL, value = 25, min = 1, width = "70px"),
        )),
      
      hr(),
      
      div(style = "text-align: center;",
          actionButton("subset_btn", "Show", #Show
                       style = "font-size: 14px; padding:10px 70px;")
      ),
      
      
      hr(),
      
      # CSS x tab
      tags$style(HTML("
    .nav-pills .nav-link {
      text-align: center;  /* Centra il testo dei titoli dei tab */
      font-size: 12px;
    }
    
    .nav-pills {
      font-size: 14px;
      display: flex;      /* Imposta flexbox per una distribuzione uniforme */
      justify-content: center; /* Centra i tab orizzontalmente */
      align-items: center; /* Centra i tab verticalmente */
    }

    /* CSS per allineare il contenuto dei tab */
    .tab-pane {
      display: flex;
      justify-content: center;  /* Centra i contenuti all'interno di ogni tab */
      align-items: center;
      text-align: center;
    }
  ")),
      
      tabsetPanel(
        type = "pills",
        selected = 'Distribution of neighbors lineages',
        
        hr(),
        
        tabPanel("Agreement matrix", 
                 #h3("Quality score of the aligment:", style = "font-size: 14px; text-align: center;"),
                 #hr(),
                 plotOutput("heatmap", dblclick = 'heatmap_click' ,height = "280px", width = "310px"),
                 tags$script(HTML("
    Shiny.addCustomMessageHandler('open_heatmap', function(url) {
      window.open(url, '_blank');
    });
  ")),
                 
                 # CSS x tooltip
                 tags$style(HTML('
  /* Personalizzazione del tooltip */
  .tooltip-inner {
    font-size: 13px;       /* Dimensione del testo */
    padding: 10px 20px;    /* Spazio interno */
    max-width: none;       /* Disabilita il limite di larghezza predefinito */
    width: 280px;          /* Imposta una larghezza fissa, puoi regolare questo valore */
  }
')),
                 
                 tags$script(HTML('
    $(document).ready(function(){
      // Tooltip per la heatmap
      $("#heatmap").tooltip({
        title: "Agreement matrix based on % of CL and tumor samples correctly aligned<br><br>Double click to open up in a new window", 
        placement: "right",
        html: true
      });
    });
  ')),
        ),
        
        hr(),
        
        tabPanel("Neighbors lineages", 
                 plotOutput("piechart",height = "220px") 
        ),
        
      )),
    
    mainPanel(
      hr(),
      uiOutput("plot"),
      width = 9)
  )
)

server <- function(input, output, session) {
  
  selected_omics_name <- reactive({
    switch(input$omics_plot,
           "Methylation" = 'Methylation',
           "Mutational signature (COSMIC)" = 'Mutational signature',
           "Expression" = 'Expression',
           "MoNETA multiomics tSNE" = 'MoNETA multiomics ',
           "MoNETA multiomics UMAP" = 'MoNETA multiomics',
           "MOFA multiomics tSNE" = 'MOFA multiomics ',
           "MOFA multiomics UMAP" = 'MOFA multiomics',
           "SNF multiomics UMAP" = 'SNF multiomics')
  })
  
  selected_plot <- reactive({
    switch(input$omics_plot,
           "Methylation" = my_plotting(meth_umap, ann_multiomics_v6, selected_omics_name()),
           "Mutational signature (COSMIC)" = my_plotting(mut_umap, ann_multiomics_v6, selected_omics_name()),
           "Expression" = my_plotting(exp_umap, ann_multiomics_v6, selected_omics_name()),
           "MoNETA multiomics tSNE" = my_plotting_tSNE(tsne_MoNETA, ann_multiomics_v6, selected_omics_name()),
           "MoNETA multiomics UMAP" = my_plotting(umap_MoNETA, ann_multiomics_v6, selected_omics_name()),
           "MOFA multiomics tSNE" = my_plotting_tSNE(tsne_MOFA, ann_multiomics_v6, selected_omics_name()),
           "MOFA multiomics UMAP" = my_plotting(umap_MOFA, ann_multiomics_v6, selected_omics_name()),
           "SNF multiomics UMAP" = my_plotting_2(SNF_umap, ann_multiomics_v6, selected_omics_name()))
  })
  
  output$plot <- renderUI({
    selected_plot()  
  })
  
  observeEvent(input$omics_plot, {
    output$plot <- renderUI({
      selected_plot() 
    })
  })
  
  
  selected_combined_mat <- reactive({
    switch(input$omics_plot,
           "Expression" = t(index_exp@data),
           "Methylation" = t(index_meth@data),
           "Mutational signature (COSMIC)" = t(index_mut@data),
           "MoNETA multiomics tSNE" = t(index_MoNETA@data),
           "MoNETA multiomics UMAP" = t(index_MoNETA@data),
           "MOFA multiomics tSNE" = t(index_MOFA@data),
           "MOFA multiomics UMAP" = t(index_MOFA@data),
           "SNF multiomics UMAP" = t(index_SNF@data)
    )
  })
  
  selected_BNindex <- reactive({
    switch(input$omics_plot,
           "Expression" = index_exp,
           "Methylation" = index_meth,
           "Mutational signature (COSMIC)" = index_mut,
           "MoNETA multiomics tSNE" = index_MoNETA,
           "MoNETA multiomics UMAP" = index_MoNETA,
           "MOFA multiomics tSNE" = index_MOFA,
           "MOFA multiomics UMAP" = index_MOFA,
           "SNF multiomics UMAP" = index_SNF
    )
  })
  
  selected_reduced_mat <- reactive({
    switch(input$omics_plot,
           "Expression" = exp_umap, 
           "Methylation" = meth_umap,
           "Mutational signature (COSMIC)" = mut_umap,
           "MoNETA multiomics tSNE" = tsne_MoNETA,
           "MoNETA multiomics UMAP" = umap_MoNETA,
           "MOFA multiomics tSNE" = tsne_MOFA,
           "MOFA multiomics UMAP" = umap_MOFA,
           "SNF multiomics UMAP" = SNF_umap
    )
  })
  
  selected_sample_order <- reactive({
    switch(input$omics_plot,
           "Expression" = sample_exp_order, 
           "Methylation" = sample_meth_order,
           "Mutational signature (COSMIC)" = sample_mut_order,
           "MoNETA multiomics tSNE" = sample_MoNETA_order,
           "MoNETA multiomics UMAP" = sample_MoNETA_order,
           "MOFA multiomics tSNE" = sample_MOFA_order,
           "MOFA multiomics UMAP" = sample_MOFA_order,
           "SNF multiomics UMAP" = sample_SNF_order
    )
  })
  
  updateSelectizeInput(session, "sel_lineage", choices = ann_multiomics_v6$lineage, server = TRUE)
  
  observeEvent(input$subset_btn, {
    if(is.null(input$both_sample) || length(input$both_sample) == 0) {
      
      output$plot <- renderUI({
        selected_plot()
      })
    }
  })
  
  observe({
    choices <- list()
    subset_1 <- ann_multiomics_v6
    
    if (input$sel_lineage == "") {
      
      if ("Cell lines" %in% input$sel_type) {
        
        cl_values <- rownames(selected_combined_mat())
        
        if (!is.null(cl_values)) {
          
          cl_values <- cl_values[!grepl('TCGA', cl_values)]
          cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v6)
          
          if (!is.null(cl_names)) {
            
            choices <- c(choices, setNames(as.list(cl_values), cl_names))
          }
        }
      }
      
      if ("Tumors" %in% input$sel_type) {
        
        tumor_values <- rownames(selected_combined_mat())
        
        if (!is.null(tumor_values)) {
          
          tumor_values <- tumor_values[grepl('TCGA', tumor_values)]
          choices <- c(choices, setNames(as.list(tumor_values), tumor_values))
        }
      }
      
    } else {
      
      if (length(input$sel_lineage) >= 1) {
        
        subset_1 <- subset_1[ann_multiomics_v6$lineage == input$sel_lineage, ]
        
        subset <- subset_1$sampleID
        all_names <- get_all_name_strpp(combined_mat = selected_combined_mat(),
                                        ann = ann_multiomics_v6,
                                        sample_info = sample_info)
        
        all_values <- rownames(selected_combined_mat())
        
        if (!is.null(all_values) && !is.null(all_names)) {
          
          all_values_x <- all_values[all_values %in% subset]
          all_names_x <- all_names[all_names %in% subset_1$stripped_cell_line_name]
          
        } else {
          
          all_values_x <- character(0)
          all_names_x <- character(0)
        }
        
        if ("Cell lines" %in% input$sel_type) {
          
          cl_values <- rownames(selected_combined_mat())
          
          if (!is.null(cl_values)) {
            
            cl_values <- cl_values[!grepl('TCGA', cl_values)]
            cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v6)
            
            if (!is.null(cl_names)) {
              
              cl_filtered <- intersect(cl_values, all_values_x)
              cl_names_filtered <- intersect(cl_names, all_names_x)
              choices <- c(choices, setNames(as.list(cl_filtered), cl_names_filtered))
            }
          }
        }
        
        if ("Tumors" %in% input$sel_type) {
          
          tumor_values <- rownames(selected_combined_mat())
          
          if (!is.null(tumor_values)) {
            
            tumor_values <- tumor_values[grepl('TCGA', tumor_values)]
            tumor_filtered <- intersect(tumor_values, all_values_x)
            choices <- c(choices, setNames(as.list(tumor_filtered), tumor_filtered))
          }
        }
        
        if (all(c("Cell lines", "Tumors") %in% input$sel_type)) {
          
          choices <- c(choices, setNames(as.list(all_values_x), all_names_x))
        }
      }
    }
    
    updateSelectizeInput(session, "both_sample", choices = choices, server = TRUE)
  })
  
  
  
  selected_samples <- reactive({
    input$both_sample
  })
  
  
  observeEvent(input$subset_btn, {
    
    if(length(input$both_sample) == 1) {
      
      if ("Cell lines" %in% input$df_selection_output) {
        
        x <- find_neighbors_CL(combined_mat = selected_combined_mat(), 
                               reduced_mat = selected_reduced_mat(),
                               input_sample = input$both_sample,
                               k = input$num_neighbors,
                               ann = ann_multiomics_v6,
                               BNindex = selected_BNindex(), 
                               sample_order = selected_sample_order(),
                               omics_name = selected_omics_name()) 
      }
      
      if ("Tumors" %in% input$df_selection_output) {
        
        x <- find_neighbors_tumor(combined_mat = selected_combined_mat(), 
                                  reduced_mat = selected_reduced_mat(),
                                  input_sample = input$both_sample,
                                  k = input$num_neighbors,
                                  ann = ann_multiomics_v6,
                                  BNindex = selected_BNindex(), 
                                  sample_order = selected_sample_order(),
                                  omics_name = selected_omics_name())
      }
      
      if (all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
        x <- find_neighbors_both(combined_mat = selected_combined_mat(), 
                                 reduced_mat = selected_reduced_mat(),
                                 input_sample = input$both_sample,
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v6,
                                 BNindex = selected_BNindex(),
                                 sample_order = selected_sample_order(),
                                 omics_name = selected_omics_name())
      }
      
      output$plot <- renderUI({
        x  
      })
      
    } else if (length(input$both_sample) > 1) {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      if("Cell lines" %in% input$df_selection_output) {
        
        x <- omics_metasample_both_CL(combined_mat = selected_combined_mat(), 
                                      reduced_mat = selected_reduced_mat(), 
                                      selected_samples = selected_samples(), 
                                      n = input$num_neighbors, 
                                      ann = ann_multiomics_v6,
                                      omics_name = selected_omics_name())
      }
      
      if("Tumors" %in% input$df_selection_output) {
        
        x <- omics_metasample_both_tumor(combined_mat = selected_combined_mat(), 
                                         reduced_mat = selected_reduced_mat(), 
                                         selected_samples = selected_samples(), 
                                         n = input$num_neighbors, 
                                         ann = ann_multiomics_v6,
                                         omics_name = selected_omics_name())
      }
      
      
      if(all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
        
        x <- omics_metasample_both_both(combined_mat = selected_combined_mat(), 
                                        reduced_mat = selected_reduced_mat(), 
                                        selected_samples = selected_samples(), 
                                        n = input$num_neighbors, 
                                        ann = ann_multiomics_v6,
                                        omics_name = selected_omics_name())
      }
      
      output$plot <- renderUI({
        x  
      })
    }
  })
  
  observeEvent(input$subset_btn, {
    
    if(length(input$both_sample) == 1) {
      
      if ("Cell lines" %in% input$df_selection_output) {
        
        piechart <- get_piechart_CL(combined_mat = selected_combined_mat(),
                                    input_sample = input$both_sample,
                                    k = input$num_neighbors,
                                    ann = ann_multiomics_v6,
                                    BNindex = selected_BNindex(),
                                    sample_order = selected_sample_order())
        
      }
      
      if ("Tumors" %in% input$df_selection_output) {
        
        piechart <- get_piechart_tumor(combined_mat = selected_combined_mat(),
                                       input_sample = input$both_sample,
                                       k = input$num_neighbors,
                                       ann = ann_multiomics_v6,
                                       BNindex = selected_BNindex(),
                                       sample_order = selected_sample_order())
      }
      
      if (all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
        
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
      
    } else if (length(input$both_sample) > 1) {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      if (all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
        
        piechart_2 <- get_piechart_both_2(combined_mat = selected_combined_mat(),
                                          selected_samples = selected_samples(),
                                          n = input$num_neighbors,
                                          ann = ann_multiomics_v6)
      }
      
      
      if ("Cell lines" %in% input$df_selection_output & !'Tumors' %in% input$df_selection_output) {
        
        piechart_2 <- get_piechart_CL_2(combined_mat = selected_combined_mat(),
                                        selected_samples = selected_samples(),
                                        n = input$num_neighbors,
                                        ann = ann_multiomics_v6)
        
      }
      
      if ("Tumors" %in% input$df_selection_output & !"CL" %in% input$df_selection_output) {
        
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
  
  
  
  observeEvent(input$omics_plot, {
    
    updateSelectizeInput(session, "both_sample", choices = selected_samples(), server = TRUE, selected = selected_samples())
    
  }, ignoreInit = TRUE) 
  
  observeEvent(input$omics_plot, {
    
    if (all(input$both_sample %in% selected_samples())) {
      
      if(length(input$both_sample) == 1) {
        
        if ("Cell lines" %in% input$df_selection_output) {
          
          x <- find_neighbors_CL(combined_mat = selected_combined_mat(),
                                 reduced_mat = selected_reduced_mat(),
                                 input_sample = input$both_sample,
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v6,
                                 BNindex = selected_BNindex(),
                                 sample_order = selected_sample_order(),
                                 omics_name = selected_omics_name())
        }
        
        if ("Tumors" %in% input$df_selection_output) {
          
          x <- find_neighbors_tumor(combined_mat = selected_combined_mat(),
                                    reduced_mat = selected_reduced_mat(),
                                    input_sample = input$both_sample,
                                    k = input$num_neighbors,
                                    ann = ann_multiomics_v6,
                                    BNindex = selected_BNindex(),
                                    sample_order = selected_sample_order(),
                                    omics_name = selected_omics_name())
        }
        
        if (all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
          x <- find_neighbors_both(combined_mat = selected_combined_mat(),
                                   reduced_mat = selected_reduced_mat(),
                                   input_sample = input$both_sample,
                                   k = input$num_neighbors,
                                   ann = ann_multiomics_v6,
                                   BNindex = selected_BNindex(),
                                   sample_order = selected_sample_order(),
                                   omics_name = selected_omics_name())
        }
        
        output$plot <- renderUI({
          x
        })
        
      } else if (length(input$both_sample) > 1) {
        
        selected_samples <- reactive({
          input$both_sample
        })
        
        if("Cell lines" %in% input$df_selection_output) {
          
          x <- omics_metasample_both_CL(combined_mat = selected_combined_mat(),
                                        reduced_mat = selected_reduced_mat(),
                                        selected_samples = selected_samples(),
                                        n = input$num_neighbors,
                                        ann = ann_multiomics_v6,
                                        omics_name = selected_omics_name())
        }
        
        if("Tumors" %in% input$df_selection_output) {
          
          x <- omics_metasample_both_tumor(combined_mat = selected_combined_mat(),
                                           reduced_mat = selected_reduced_mat(),
                                           selected_samples = selected_samples(),
                                           n = input$num_neighbors,
                                           ann = ann_multiomics_v6,
                                           omics_name = selected_omics_name())
        }
        
        
        if(all(c("Cell lines", "Tumors") %in% input$df_selection_output)) {
          
          x <- omics_metasample_both_both(combined_mat = selected_combined_mat(),
                                          reduced_mat = selected_reduced_mat(),
                                          selected_samples = selected_samples(),
                                          n = input$num_neighbors,
                                          ann = ann_multiomics_v6,
                                          omics_name = selected_omics_name())
        }
        
        output$plot <- renderUI({
          x
        })
      }
    }
  }, ignoreInit = TRUE)
  
  
  
  selected_heatmap <- reactive({
    switch(input$omics_plot,
           "Expression" = get_confusion_matrix_shiny(exp_m3),
           "Methylation" = get_confusion_matrix_shiny(meth_m3),
           "Mutational signature (COSMIC)" = get_confusion_matrix_shiny(mut_m3),
           "MoNETA multiomics tSNE" = get_confusion_matrix_shiny(MoNETA_m3),
           "MoNETA multiomics UMAP" = get_confusion_matrix_shiny(MoNETA_m3),
           "MOFA multiomics tSNE" = get_confusion_matrix_shiny(MOFA_m3),
           "MOFA multiomics UMAP" = get_confusion_matrix_shiny(MOFA_m3),
           "SNF multiomics UMAP" = get_confusion_matrix_shiny(SNF_m3))
  })
  
  output$heatmap <- renderPlot({
    heatmap_plot <- selected_heatmap()
    heatmap_plot
  })
  
  observeEvent(input$heatmap_click, {
    req(selected_heatmap())
    file_path <- file.path(tempdir(), "heatmap.png")
    
    if (!dir.exists("www")) {
      dir.create("www")
    }
    
    ggsave(file_path, plot = selected_heatmap(), width = 10, height = 7, dpi = 300)
    
    file.copy(file_path, "www/heatmap.png", overwrite = TRUE)
    
    session$sendCustomMessage("open_heatmap", "heatmap.png")
  })
  
  filtered_data <- reactiveVal()
  
  observeEvent(event_data("plotly_selected"), {
    selected_data <- event_data("plotly_selected")
    
    if (!is.null(selected_data)) {
      
      selected_samples <- selected_data$key
      
      x_1 <- ann_multiomics_v6 %>%
        filter(sampleID %in% selected_samples)
      
      filtered_data(x_1$stripped_cell_line_name)
    }
  })
  
  lasso_selected_samples <- reactive({
    req(filtered_data()) 
    
    filtered_data()
  })
  
  observeEvent(input$load_selection, {
    
    updateSelectizeInput(session, "both_sample",
                         choices = lasso_selected_samples(),
                         selected = lasso_selected_samples())
  })
  
}

shiny::shinyApp(ui, server)



