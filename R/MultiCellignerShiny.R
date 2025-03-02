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
      
      tags$style(HTML("
  /* Forza altezza fissa del box di input solo per 'both_sample' */
  #both_sample + .selectize-control.multi .selectize-input {
    min-height: 190px !important; 
    max-height: 190px !important; 
    overflow-y: auto !important;  /* Attiva lo scroll verticale */
    overflow-x: hidden !important; /* Disabilita lo scroll orizzontale */
    display: flex !important;
    align-items: flex-start !important; /* Allinea il testo in cima */
    flex-wrap: wrap !important; /* Permette ai campioni di andare a capo */
  }

  /* Imposta altezza massima per il menu dropdown di 'both_sample' */
  #both_sample + .selectize-control.multi .selectize-dropdown-content {
    max-height: 200px !important;  /* Altezza massima del menu a tendina */
    overflow-y: auto !important;   /* Attiva lo scroll verticale */
    overflow-x: hidden !important; /* Disabilita lo scroll orizzontale */
  }

  /* Adatta la larghezza del menu dropdown per 'both_sample' */
  #both_sample + .selectize-control.multi .selectize-dropdown {
    width: 100% !important; /* Adatta la larghezza */
  }
")),
      
      
      
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
      
      tags$script(HTML('
        $(document).ready(function(){
          // Tooltip per il bottone "Show"
          $("#subset_btn").tooltip({
            title: "Click on Show or switch omics without sample/s in search bar to get the basic plot without neighbors",
            placement: "right",
            trigger: "hover",
            html: true
          });
        });
      ')),
      
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
        selected = 'Neighbors lineages',
        
        hr(),
        
        tabPanel("Neighbors lineages", 
                  plotOutput("piechart",height = "220px") , 

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
                 ),
        
        hr(),
        
        tabPanel("Neighbors subtypes", 
                 plotOutput("piechart_subtype",height = "220px"),
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
  
  selected_plot <- reactive({
    switch(input$omics_plot,
           "Methylation" = my_plotting(meth_umap, ann_multiomics_v6, selected_omics_name()),
           "Mutational signature (COSMIC)" = my_plotting(mut_umap, ann_multiomics_v6, selected_omics_name()),
           "Expression" = my_plotting(exp_umap, ann_multiomics_v6, selected_omics_name()),
           "MoNETA multiomics tSNE" = my_plotting_tSNE(tsne_MoNETA, ann_multiomics_v6, selected_omics_name()),
           "MoNETA multiomics UMAP" = my_plotting(umap_MoNETA, ann_multiomics_v6, selected_omics_name()),
           "MOFA multiomics tSNE" = my_plotting_tSNE(tsne_MOFA, ann_multiomics_v6, selected_omics_name()),
           "MOFA multiomics UMAP" = my_plotting(umap_MOFA, ann_multiomics_v6, selected_omics_name()),
           "SNF multiomics UMAP" = my_plotting(SNF_umap, ann_multiomics_v6, selected_omics_name()))
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
           "Expression" = pca_exp_1,
           "Methylation" = pca_meth_1,
           "Mutational signature (COSMIC)" = pca_mut_1,
           "MoNETA multiomics tSNE" = pca_MoNETA,
           "MoNETA multiomics UMAP" = pca_MoNETA,
           "MOFA multiomics tSNE" = pca_MOFA,
           "MOFA multiomics UMAP" = pca_MOFA,
           "SNF multiomics UMAP" = pca_SNF
    )
  })
  
  ### get the lineage in the menu
  updateSelectizeInput(session, "sel_lineage", choices = ann_multiomics_v6$lineage, selected = character(0))
  
  ### when there is no sample/s in search bar, reolad the omics base plot
  observeEvent(input$subset_btn, {
    if(is.null(input$both_sample) || length(input$both_sample) == 0) {
      
      output$plot <- renderUI({
        selected_plot()
      })
    }
  })
  
  
  ### update the sample/s for each selected omics; using depmap code but show nameID
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
  
  #### click show in the shiny to find neighbors
  observeEvent(input$subset_btn, {
    
    if(length(input$both_sample) == 1) {
      
      x <- find_neighbors(combined_mat = selected_combined_mat(), 
                          reduced_mat = selected_reduced_mat(),
                          input_sample = input$both_sample,
                          k = input$num_neighbors,
                          ann = ann_multiomics_v6,
                          type = input$df_selection_output,
                          omics_name = selected_omics_name())

      output$plot <- renderUI({
        x  
      })
      
    } else if (length(input$both_sample) > 1) {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      x <- find_neighbors(combined_mat = selected_combined_mat(), 
                          reduced_mat = selected_reduced_mat(),
                          selected_samples = selected_samples(),
                          k = input$num_neighbors,
                          ann = ann_multiomics_v6,
                          type = input$df_selection_output,
                          omics_name = selected_omics_name())
      
     
      output$plot <- renderUI({
        x  
      })
    }
  })
  
  #### click show in the shiny to get both kind of piechart
  observeEvent(input$subset_btn, {
    
    if(length(input$both_sample) == 1) {
      
     piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                input_sample = input$both_sample,
                                k = input$num_neighbors,
                                ann = ann_multiomics_v6,
                                type = input$df_selection_output)
     
     piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(), 
                                              input_sample = input$both_sample,
                                              k = input$num_neighbors,
                                              ann = ann_multiomics_v6,
                                              type = input$df_selection_output)
     
     output$piechart_subtype <- renderPlot({
       piechart_subtype
     })

      output$piechart <- renderPlot({
        piechart
      })
      
    } else if (length(input$both_sample) > 1) {
      
      selected_samples <- reactive({
        input$both_sample
      })
      
      piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                               selected_samples = selected_samples(),
                               k = input$num_neighbors,
                               ann = ann_multiomics_v6,
                               type = input$df_selection_output)
      
      piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(), 
                                               selected_samples = selected_samples(),
                                               k = input$num_neighbors,
                                               ann = ann_multiomics_v6,
                                               type = input$df_selection_output)
      
      
      output$piechart_subtype <- renderPlot({
        piechart_subtype
      })
      
      output$piechart <- renderPlot({
        piechart
      })
    }
  })
  
  ### keep the sample/s in the search bar when you switch omics
  observeEvent(input$omics_plot, {
    
    updateSelectizeInput(session, "both_sample", choices = selected_samples(), server = TRUE, selected = selected_samples())
    
  }, ignoreInit = TRUE) 
  
  ### when you already find neighbors for a sample/s just switching plot you will find neighbors for other omics
  observeEvent(input$omics_plot, { 
    
    if (all(input$both_sample %in% selected_samples())) {
      
       if(length(input$both_sample) == 1) {
          
         x <- find_neighbors(combined_mat = selected_combined_mat(), 
                             reduced_mat = selected_reduced_mat(),
                             input_sample = input$both_sample,
                             k = input$num_neighbors,
                             ann = ann_multiomics_v6,
                             type = input$df_selection_output,
                             omics_name = selected_omics_name())
         
         piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                  input_sample = input$both_sample,
                                  k = input$num_neighbors,
                                  ann = ann_multiomics_v6,
                                  type = input$df_selection_output)

        
        output$plot <- renderUI({
          x  
        })
        
        output$piechart <- renderPlot({
          piechart
        })
        
      } else if (length(input$both_sample) > 1) {
        
        selected_samples <- reactive({
          input$both_sample
        })
        
        x <- find_neighbors(combined_mat = selected_combined_mat(), 
                            reduced_mat = selected_reduced_mat(),
                            selected_samples = selected_samples(),
                            k = input$num_neighbors,
                            ann = ann_multiomics_v6,
                            type = input$df_selection_output,
                            omics_name = selected_omics_name())
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                 selected_samples = selected_samples(),
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v6,
                                 type = input$df_selection_output)
        
        output$plot <- renderUI({
          x  
        })
        
        output$piechart <- renderPlot({
          piechart
        })
      }
    }
  }, ignoreInit = TRUE)
  
  #### lasso select code
  
  filtered_data <- reactiveVal()
  
  observeEvent(event_data("plotly_selected"), {
    selected_data <- event_data("plotly_selected")
    
    if (!is.null(selected_data)) {
      selected_samples <- selected_data$key
      
      x_1 <- ann_multiomics_v6 %>%
        filter(sampleID %in% selected_samples)
      
      filtered_data(x_1$sampleID)
    }
  })
  
  lasso_selected_samples <- reactive({
    req(filtered_data()) 
    
    filtered_data()
  })
  
  observeEvent(input$load_selection, {
    
    choices <- list()
    
    all_names <- ann_multiomics_v6$stripped_cell_line_name
    all_values <- ann_multiomics_v6$sampleID
    
    sub_values <- all_values[all_values %in% lasso_selected_samples()]
    sub_names <- ann_multiomics_v6 %>% filter(sampleID %in% sub_values)
    sub_names <- sub_names$stripped_cell_line_name
    
    choices <- c(choices, setNames(as.list(sub_values), sub_names))
  
    updateSelectizeInput(session, "both_sample",
                         choices = choices,
                         selected = choices)
  })
  
}

shiny::shinyApp(ui, server)



