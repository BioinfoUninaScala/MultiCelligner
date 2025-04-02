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
#'

MultiCellignerShiny <- function() {shiny::shinyApp(ui, server)}

ui <- fluidPage(              ######## DEVI RIMUOVERE DALLE UMAP DI MONETA I CAMPIONI NON ANNOTATI!!!!!!!!!!
  
  shinyjs::useShinyjs(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      tags$style(HTML("hr { margin-top: 4px !important; margin-bottom: 4px !important; padding: 4px !important; }")),
      
      div(h3("MultiCelligner"), class = "text-center"),
      hr(),
      
      fluidRow(
        column(12, 
               div(style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   selectInput('omics_plot', 'Select an omics alignment:',
                               choices = c('Methylation',
                                           "Mutational signature (COSMIC)",
                                           "Expression"),
                               multiple = TRUE,
                               selected = 'Methylation',
                               selectize = TRUE))
        )),
      
      fluidRow(
        column(6, 
               div(style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   selectInput('reduction_method', HTML("Select<br>reduction:"),
                               choices = c('UMAP', "tSNE"),
                               selected = 'UMAP',
                               width = "150px"))
        ),
        
        column(6,
               div(style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   selectInput('multiomics_method', 'Select integration method:',
                               choices = NULL))
        )),
      
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
    min-height: 130px !important; 
    max-height: 130px !important; 
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
      
      fluidRow(
        column(12,
               selectizeInput("lin_output", 
                              "Limit query lineage/s to:", 
                              choices = NULL, 
                              multiple = TRUE)
        )
      ),
      
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
    
    if(length(input$omics_plot) == 1) {
      
      switch(input$omics_plot,
             "Methylation" = 'Methylation',
             "Mutational signature (COSMIC)" = 'Mutational signature',
             "Expression" = 'Expression')
    }
    
    if(length(input$omics_plot) > 1) {
      
      switch(input$multiomics_method,
             "MoNETA" = 'MoNETA multiomics',
             "MOFA" = 'MOFA multiomics',
             "SNF" = 'SNF multiomics')
    }
    
  })
  
  selected_reduced_mat <- reactive({
    
    if(length(input$omics_plot) == 1) {
      
      return(switch(input$omics_plot,
                    "Expression" = exp_umap, 
                    "Methylation" = meth_umap,
                    "Mutational signature (COSMIC)" = mut_umap))
      
    } else if (length(input$omics_plot) > 1) {
      
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(umap_exp_meth_mut)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(umap_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(umap_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(umap_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(tsne_exp_meth_mut)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(tsne_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(tsne_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(tsne_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(mofa_umap_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(mofa_umap_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(mofa_umap_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(mofa_umap_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(mofa_tsne_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(mofa_tsne_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(mofa_tsne_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(mofa_tsne_meth_mut)
        }
      }
    }
    
    return(NULL)  
  }) 
  
  
  selected_plot <- reactive({
    
    if(length(input$omics_plot) == 1) {
      
      return(switch(input$omics_plot,
                    "Expression" = my_plotting(exp_umap, ann_multiomics_v8, selected_omics_name()), 
                    "Methylation" = my_plotting(meth_umap, ann_multiomics_v8, selected_omics_name()),
                    "Mutational signature (COSMIC)" = my_plotting(mut_umap, ann_multiomics_v8, selected_omics_name())))
      
    } else if (length(input$omics_plot) > 1) {
      
      # --- MoNETA - UMAP ---
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting(umap_exp_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(my_plotting(umap_exp_meth, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting(umap_exp_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(my_plotting(umap_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
      }
      
      # --- MoNETA - tSNE ---
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(tsne_exp_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(my_plotting_tSNE(tsne_exp_meth, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(tsne_exp_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(tsne_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
      }
      
      # --- MOFA - UMAP ---
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting(mofa_umap_all, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(my_plotting(mofa_umap_exp_meth, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting(mofa_umap_exp_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(my_plotting(mofa_umap_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
      }
      
      # --- MOFA - tSNE ---
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(mofa_tsne_all, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(my_plotting_tSNE(mofa_tsne_exp_meth, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(mofa_tsne_exp_mut, ann_multiomics_v8, selected_omics_name()))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(my_plotting_tSNE(mofa_tsne_meth_mut, ann_multiomics_v8, selected_omics_name()))
        }
      }
    }
    
    return(NULL)  
  })
  
  
  selected_combined_mat <- reactive({
    
    if(length(input$omics_plot) == 1) {
      
      switch(input$omics_plot,
             "Expression" = pca_exp, 
             "Methylation" = pca_meth_1,
             "Mutational signature (COSMIC)" = combined_mat_mut
      )
      
    } else if (length(input$omics_plot) > 1) {
      
      if (input$multiomics_method == 'MoNETA') {
        
        if (setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(pca_rwr_exp_meth_mut)
        } else if (setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(pca_rwr_exp_meth)
        } else if (setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(pca_rwr_exp_mut)
        } else if (setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(pca_rwr_meth_mut)
        }
      }
      
      if (input$multiomics_method == 'MOFA') {
        
        if (setequal(input$omics_plot, c("Methylation","Expression","Mutational signature (COSMIC)"))) {
          return(MOFA_mat_all)
        } else if (setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(MOFA_mat_exp_meth)
        } else if (setequal(input$omics_plot, c("Expression","Mutational signature (COSMIC)"))) {
          return(MOFA_mat_exp_mut)
        } else if (setequal(input$omics_plot, c("Methylation", "Mutational signature (COSMIC)"))) {
          return(MOFA_mat_meth_mut)
        }
        
      }
    }
  })
  
  observeEvent(list(input$reduction_method, input$multiomics_method, input$omics_plot), {
    output$plot <- renderUI({
      selected_plot() 
    })
  })
  
  ### get the multiomics data intergration method in the menu
  updateSelectizeInput(session, "multiomics_method", choices = c('MoNETA','MOFA','SNF'), selected = character(0))
  
  ### select the linage based on the mat that the user choose
  selected_linages <- reactiveVal(NULL)
  
  observeEvent(list(input$multiomics_method, input$omics_plot), {
    lin <- ann_multiomics_v8$lineage[ann_multiomics_v8$sampleID %in% rownames(selected_combined_mat())]
    selected_linages(lin)
    updateSelectizeInput(session, "sel_lineage", choices = lin, selected = character(0))
  })
  
  ### select the Query lineage/s for neighbors search
  query_linages <- reactiveVal(NULL)
  
  observe({
    lin_out <- ann_multiomics_v8$lineage[ann_multiomics_v8$sampleID %in% rownames(selected_combined_mat())]
    query_linages(lin_out)
    updateSelectizeInput(session, "lin_output", choices = c('All',lin_out), selected = 'All' ,server = TRUE)
  }) 
  
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
    subset_1 <- ann_multiomics_v8
    
    if (input$sel_lineage == "") {
      
      if ("Cell lines" %in% input$sel_type) {
        
        cl_values <- rownames(selected_combined_mat())
        
        if (!is.null(cl_values)) {
          
          cl_values <- cl_values[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', cl_values)]
          cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v8)
          
          if (!is.null(cl_names)) {
            
            choices <- c(choices, setNames(as.list(cl_values), cl_names))
          }
        }
      }
      
      if ("Tumors" %in% input$sel_type) {
        
        tumor_values <- rownames(selected_combined_mat())
        
        if (!is.null(tumor_values)) {
          
          tumor_values <- tumor_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', tumor_values)]
          choices <- c(choices, setNames(as.list(tumor_values), tumor_values))
        }
      }
      
      if (all(c("Cell lines", "Tumors") %in% input$sel_type)) {
        
        all_values <- rownames(selected_combined_mat())
        
        if (!is.null(all_values)) {
          
          cl_values <- all_values[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', all_values)]
          tumor_values <- all_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', all_values)]
          
          cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v8)
          
          if (!is.null(cl_values) && !is.null(cl_names)) {
            choices <- c(choices, setNames(as.list(cl_values), cl_names))
          }
          
          if (!is.null(tumor_values)) {
            choices <- c(choices, setNames(as.list(tumor_values), tumor_values))
          }
        }
      }
      
    } else {
      
      if (length(input$sel_lineage) >= 1) {
        
        subset_1 <- subset_1[ann_multiomics_v8$lineage == input$sel_lineage, ]
        subset_1 <- subset_1[subset_1$sampleID %in% rownames(selected_combined_mat()),]
        subset <- subset_1$sampleID
        
        all_names <- subset_1$stripped_cell_line_name[subset_1$sampleID %in% rownames(selected_combined_mat())]
        
        all_values <- rownames(selected_combined_mat())[rownames(selected_combined_mat()) %in% subset]
        
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
            
            cl_values <- cl_values[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', cl_values)]
            cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v8)
            
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
            
            tumor_values <- tumor_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', tumor_values)]
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
                          ann = ann_multiomics_v8,
                          type = input$df_selection_output,
                          omics_name = selected_omics_name(),
                          red_method = input$reduction_method,
                          query_lineage = input$lin_output)
      
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
                          ann = ann_multiomics_v8,
                          type = input$df_selection_output,
                          omics_name = selected_omics_name(),
                          red_method = input$reduction_method,
                          query_lineage = input$lin_output)
      
      
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
                               ann = ann_multiomics_v8,
                               type = input$df_selection_output,
                               query_lineage = input$lin_output)
      
      piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(), 
                                               input_sample = input$both_sample,
                                               k = input$num_neighbors,
                                               ann = ann_multiomics_v8,
                                               type = input$df_selection_output,
                                               query_lineage = input$lin_output)
      
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
                               ann = ann_multiomics_v8,
                               type = input$df_selection_output,
                               query_lineage = input$lin_output)
      
      piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(), 
                                               selected_samples = selected_samples(),
                                               k = input$num_neighbors,
                                               ann = ann_multiomics_v8,
                                               type = input$df_selection_output,
                                               query_lineage = input$lin_output)
      
      
      output$piechart_subtype <- renderPlot({
        piechart_subtype
      })
      
      output$piechart <- renderPlot({
        piechart
      })
    }
  })
  
  ### keep the sample/s in the search bar when you switch omics
  # observeEvent(list(input$reduction_method, input$multiomics_method, input$omics_plot), {
  #   
  #   updateSelectizeInput(session, "both_sample", choices = selected_samples(), server = TRUE, selected = selected_samples())
  #   
  # }, ignoreInit = TRUE) 
  
  ### value for: when you already find neighbors for a sample/s just switching plot you will find neighbors for other omics
  # selected_lin_out <- reactive({
  #   input$lin_output
  # })
  
  ### debug
  observeEvent(input$omics_plot, {
    if (length(input$omics_plot) > 1 & (is.null(input$multiomics_method) || input$multiomics_method == "")) {
      return()
    }
  }, ignoreInit = TRUE)
  
  ### when you already find neighbors for a sample/s just switching plot you will find neighbors for other omics
  observeEvent(list(input$reduction_method, input$multiomics_method, input$omics_plot), {
    
    if (is.null(input$omics_plot) || length(input$omics_plot) == 0) {
      return()
    }
    
    if (length(input$omics_plot) > 1 & (is.null(input$multiomics_method) || input$multiomics_method == "")) {
      return()
    }
    
    if (all(input$both_sample %in% selected_samples())) {
      
      if(length(input$both_sample) == 1) {
        
        selected_samples <- reactive({
          input$both_sample
        })
        
        selected_lin_out <- reactive({
          input$lin_output
        })
        
        selected_lineage1 <- reactive({
          input$sel_lineage
        })
        
        x <- find_neighbors(combined_mat = selected_combined_mat(),
                            reduced_mat = selected_reduced_mat(),
                            input_sample = input$both_sample,
                            k = input$num_neighbors,
                            ann = ann_multiomics_v8,
                            type = input$df_selection_output,
                            omics_name = selected_omics_name(),
                            red_method = input$reduction_method,
                            query_lineage = input$lin_output)
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(),
                                 input_sample = input$both_sample,
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v8,
                                 type = input$df_selection_output,
                                 query_lineage = input$lin_output)
        
        piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(),
                                                 input_sample = input$both_sample,
                                                 k = input$num_neighbors,
                                                 ann = ann_multiomics_v8,
                                                 type = input$df_selection_output,
                                                 query_lineage = input$lin_output)
        
        updateSelectizeInput(session, "both_sample", choices = selected_samples(), server = TRUE, selected = selected_samples())
        updateSelectizeInput(session, "lin_output", choices = selected_lin_out(), selected = selected_lin_out() ,server = TRUE)
        updateSelectizeInput(session, "sel_lineage", choices = selected_lineage1(), selected = selected_lineage1())
        
        output$piechart_subtype <- renderPlot({
          piechart_subtype
        })
        
        
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
        
        selected_lin_out <- reactive({
          input$lin_output
        })
        
        selected_lineage1 <- reactive({
          input$sel_lineage
        })
        
        x <- find_neighbors(combined_mat = selected_combined_mat(),
                            reduced_mat = selected_reduced_mat(),
                            selected_samples = selected_samples(),
                            k = input$num_neighbors,
                            ann = ann_multiomics_v8,
                            type = input$df_selection_output,
                            omics_name = selected_omics_name(),
                            red_method = input$reduction_method,
                            query_lineage = input$lin_output)
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(),
                                 selected_samples = selected_samples(),
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v8,
                                 type = input$df_selection_output,
                                 query_lineage = input$lin_output)
        
        piechart_subtype <- get_piechart_subtype(combined_mat = selected_combined_mat(),
                                                 selected_samples = selected_samples(),
                                                 k = input$num_neighbors,
                                                 ann = ann_multiomics_v8,
                                                 type = input$df_selection_output,
                                                 query_lineage = input$lin_output)
        
        updateSelectizeInput(session, "both_sample", choices = selected_samples(), server = TRUE, selected = selected_samples())
        updateSelectizeInput(session, "lin_output", choices = selected_lin_out(), selected = selected_lin_out() ,server = TRUE)
        updateSelectizeInput(session, "sel_lineage", choices = selected_lineage1(), selected = selected_lineage1())
        
        output$piechart_subtype <- renderPlot({
          piechart_subtype
        })
        
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
      
      x_1 <- ann_multiomics_v8 %>%
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
    
    all_names <- ann_multiomics_v8$stripped_cell_line_name
    all_values <- ann_multiomics_v8$sampleID
    
    sub_values <- all_values[all_values %in% lasso_selected_samples()]
    sub_names <- ann_multiomics_v8 %>% filter(sampleID %in% sub_values)
    sub_names <- sub_names$stripped_cell_line_name
    
    choices <- c(choices, setNames(as.list(sub_values), sub_names))
    
    updateSelectizeInput(session, "both_sample",
                         choices = choices,
                         selected = choices)
  })
  
}

shiny::shinyApp(ui, server)



