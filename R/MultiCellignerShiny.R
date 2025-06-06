#' 
#' MultiCelligner Shiny App
#'
#' @import shinycssloaders
#' @import shiny
#' @import magrittr
#' @import crosstalk
#' @import htmltools
#' @import fontawesome
#' @import shinyjs
#' @return Shiny app
#' @export
#'

MultiCellignerShiny <- function() {addResourcePath("static", system.file("www", package = "MultiCelligner"))
  ;shiny::shinyApp(ui, server)}

ui <- fluidPage(              
  
  shinyjs::useShinyjs(),
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      
      tags$style(HTML("hr { margin-top: 4px !important; margin-bottom: 4px !important; padding: 4px !important; }")),
      
      shiny::fluidRow(shiny::column(3,tags$img(src = "static/MultiCelligner_Logo_2.png", height = "100px"),), 
                      shiny::column(6,div(h3("MultiCelligner"), class = "text-center",style="padding:15px;")),
                      #column(3,tags$img(src = "Mcell_logo.png", height = "100px"))
      ),
      hr(),
      
      shiny::fluidRow(
        shiny::column(12, 
                      div(style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                          selectInput('omics_plot', 'Omics alignment:',
                                      choices = c('Methylation',
                                                  "Mutational signature",
                                                  "Expression"),
                                      multiple = TRUE
                          ))
        )),
      
      shiny::fluidRow(
        shiny::column(6, 
                      div(style = "display: flex; flex-direction: shiny::column; align-items: center; justify-content: center; height: 100%;",
                          selectInput('reduction_method', "Reduction:",
                                      choices = c('UMAP', "tSNE"),
                                      selected = 'UMAP',
                                      width = "150px"))
        ),
        
        shiny::column(6,
                      div(style = "display: flex; flex-direction: shiny::column; align-items: center; justify-content: center; height: 100%;",
                          selectInput('multiomics_method', 'Integration method:',
                                      choices = NULL))
        )),
      
      hr(),
      
      tags$p("Find neighbors", style = "text-align: center; font-size: 18px; font-weight: bold;"),
      
      shiny::fluidRow(
        shiny::column(6,
                      selectizeInput("sel_type",
                                     'Select model type',
                                     choices = NULL, 
                                     multiple = TRUE)),
        shiny::column(6,
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
      
      shiny::fluidRow(
        shiny::column(6,
                      tags$strong("or load from map selection:")),
        shiny::column(2,
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
        ),
        shiny::column(4,
                      shiny::div(style = "align-items: flex-start; justify-content: center;",
                                 actionButton("rm", "Clear Selection", style = "text-align: center;")) %>%
                        tagAppendChild(tags$script(HTML('
    $(document).ready(function(){
      $("#rm").tooltip({
        title: "Remove all the samples",
        placement: "right",
        trigger: "hover",
        html: true
      });
    });
  ')))
        )
      ),
      
      
      hr(),
      
      shiny::fluidRow(
        shiny::column(7,
                      checkboxGroupInput("df_selection_output", 
                                         "Search among the closest:", 
                                         choices = c("Cell lines", "Tumors"), 
                                         selected = 'Tumors', inline = TRUE)),
        shiny::column(5,
                      div("Number of neighbors:", style = "text-align: left; font-weight: bold;"),
                      numericInput("num_neighbors", NULL, value = 25, min = 1, width = "70px"),
        )),
      
      shiny::fluidRow(
        shiny::column(12,
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
            title: "Click on Plot Alignment without sample/s in search bar to get the basic plot without neighbors",
            placement: "right",
            trigger: "hover",
            html: true
          });
        });
      ')),
      
      div(style = "text-align: center;",
          actionButton("subset_btn", "Plot alignment", 
                       style = "
                   font-size: 14px;
                   padding: 15px 70px;
                   background-color: #4FC3F7; /* Azzurro chiaro */
                   color: solid black;             /* Colore del testo */
                   font-weight: bold;        /* Testo in grassetto */
                   border: none;             /* Nessun bordo */
                   border-radius: 6px;       /* Angoli leggermente stondati */
                 ")
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
        type = "tabs",
        selected = 'Neighbors lineages',
        
        hr(),
        
        tabPanel("Neighbors lineages", 
                 plotOutput("piechart",height = "430px") , 
                 
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
        
        #hr(),
        
        tabPanel("Neighbors subtypes", 
                 plotOutput("piechart_subtype",height = "430px"),
        ),
        
      )),
    
    mainPanel(
      tags$div(
        tags$p(h3(textOutput("omics_name"),class = "text-center"), style = "text-align: center; font-size: 18px; font-weight: bold;"),
      ),
      
      shinycssloaders::withSpinner(uiOutput("plot"), type = 6, color = "#4FC3F7", color.background = "white", ), 
      width = 9,
      
      hr(),
      
      tabsetPanel(
        type = "tabs",
        selected = 'Lineage distribution',
        
        tabPanel("Lineage distribution",          
                 plotlyOutput("distribution_s", width = '90%')
        ),
        
        tabPanel("Subtype distribution", 
                 plotlyOutput("distribution_l", width = "90%")
        ))
      
    )
  )
)                                


server <- function(input, output, session) { 
  
  
  output$omics_name <- eventReactive(input$subset_btn,{ 
    if(length(input$omics_plot) == 1) 
      paste(input$reduction_method,paste(input$omics_plot, collapse=" "))
    else paste(input$reduction_method,paste(input$omics_plot, collapse=" "), input$multiomics_method)
  })
  
  selected_omics_name <- reactive({
    
    if(length(input$omics_plot) == 1) {
      
      return(switch(input$omics_plot,
                    "Methylation" = 'Methylation',
                    "Mutational signature" = 'Mutational signature',
                    "Expression" = 'Expression'))
    }
    
    if(length(input$omics_plot) > 1) {
      
      return(switch(input$multiomics_method,
                    "MoNETA" = 'MoNETA multiomics',
                    "MOFA" = 'MOFA multiomics',
                    "SNF" = 'SNF multiomics'))
    }
    
  })
  
  
  selected_reduced_mat <- reactive({
    
    if(length(input$omics_plot) == 1 & input$reduction_method == 'UMAP') {
      
      return(switch(input$omics_plot,
                    "Expression" = exp_umap, 
                    "Methylation" = meth_umap,
                    "Mutational signature" = mut_umap))
    }
    
    if(length(input$omics_plot) == 1 & input$reduction_method == 'tSNE') {
      
      return(switch(input$omics_plot,
                    "Expression" = tsne_exp, 
                    "Methylation" = tsne_meth,
                    "Mutational signature" = tsne_mut))
      
    } else if (length(input$omics_plot) > 1) {
      
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(umap_exp_meth_mut)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(umap_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(umap_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(umap_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(tsne_exp_meth_mut)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(tsne_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(tsne_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(tsne_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(mofa_umap_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(mofa_umap_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(mofa_umap_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(mofa_umap_meth_mut)
        }
        
      }
      
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(mofa_tsne_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(mofa_tsne_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(mofa_tsne_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(mofa_tsne_meth_mut)
        }
      }
      
      if(input$multiomics_method == 'SNF' && input$reduction_method == 'UMAP') {
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(snf_umap_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(snf_umap_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(snf_umap_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(snf_umap_meth_mut)
        }
      }
      
      if(input$multiomics_method == 'SNF' && input$reduction_method == 'tSNE') {
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(snf_tsne_all)
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(snf_tsne_exp_meth)
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(snf_tsne_exp_mut)
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(snf_tsne_meth_mut)
        }
      }
      
    }
    
    return(NULL)  
  }) 
  
  
  selected_plot <- eventReactive(input$subset_btn, {
    
    if(length(input$omics_plot) == 1 & input$reduction_method == 'UMAP') {
      
      return(switch(input$omics_plot,
                    "Expression" = get_alignment_plot(reduced_mat = exp_umap, ann = ann_multiomics_v9),
                    "Methylation" = get_alignment_plot(reduced_mat = meth_umap, ann = ann_multiomics_v9),
                    "Mutational signature" = get_alignment_plot(reduced_mat = mut_umap, ann = ann_multiomics_v9)))
      
    }
    
    if (length(input$omics_plot) == 1 & input$reduction_method == 'tSNE'){
      
      return(switch(input$omics_plot,
                    "Expression" = get_alignment_plot(reduced_mat = tsne_exp, ann = ann_multiomics_v9),
                    "Methylation" = get_alignment_plot(reduced_mat = tsne_meth, ann = ann_multiomics_v9),
                    "Mutational signature" = get_alignment_plot(reduced_mat = tsne_mut, ann = ann_multiomics_v9)))
      
    } else if (length(input$omics_plot) > 1) {
      
      # --- MoNETA - UMAP ---
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = umap_exp_meth_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = umap_exp_meth, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = umap_exp_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = umap_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
      # --- MoNETA - tSNE ---
      if(input$multiomics_method == 'MoNETA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = tsne_exp_meth_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = tsne_exp_meth, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = tsne_exp_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = tsne_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
      # --- MOFA - UMAP ---
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'UMAP') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_umap_all, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = mofa_umap_exp_meth, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_umap_exp_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_umap_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
      # --- MOFA - tSNE ---
      if(input$multiomics_method == 'MOFA' && input$reduction_method == 'tSNE') {
        
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_tsne_all, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = mofa_tsne_exp_meth, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_tsne_exp_mut, ann = ann_multiomics_v9))
        }
        else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = mofa_tsne_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
      # --- SNF - UMAP ---
      if(input$multiomics_method == 'SNF' && input$reduction_method == 'UMAP') {
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_umap_all, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = snf_umap_exp_meth, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_umap_exp_mut, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_umap_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
      # --- SNF - tSNE ---
      if(input$multiomics_method == 'SNF' && input$reduction_method == 'tSNE') {
        if(setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_tsne_all, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Methylation","Expression"))) {
          return(get_alignment_plot(reduced_mat = snf_tsne_exp_meth, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_tsne_exp_mut, ann = ann_multiomics_v9))
        } else if(setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          return(get_alignment_plot(reduced_mat = snf_tsne_meth_mut, ann = ann_multiomics_v9))
        }
      }
      
    }
    
    return(NULL)  
  })
  
  
  last_combined_mat <- reactiveVal(NULL)
  
  selected_combined_mat <- reactive({
    
    ### if there is no omics selected will be use the last combined_mat: 
    ### this is for keep the sample in the search bar
    if (length(input$omics_plot) == 0) {
      return(last_combined_mat())
    }
    
    mat <- NULL
    
    if (length(input$omics_plot) == 1) {
      mat <- switch(input$omics_plot,
                    "Expression" = pca_exp, 
                    "Methylation" = pca_meth_1,
                    "Mutational signature" = combined_mat_mut)
    } else if (length(input$omics_plot) > 1) {
      
      if (input$multiomics_method == 'MoNETA') {
        if (setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          mat <- emb_exp_meth_mut_1
        } else if (setequal(input$omics_plot, c("Methylation","Expression"))) {
          mat <- emb_exp_meth_1
        } else if (setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          mat <- emb_exp_mut_1
        } else if (setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          mat <- emb_meth_mut_1
        }
      }
      
      if (input$multiomics_method == 'MOFA') {
        if (setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          mat <- MOFA_mat_all
        } else if (setequal(input$omics_plot, c("Methylation","Expression"))) {
          mat <- MOFA_mat_exp_meth
        } else if (setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          mat <- MOFA_mat_exp_mut
        } else if (setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          mat <- MOFA_mat_meth_mut
        }
      }
      
      if (input$multiomics_method == 'SNF') {
        if (setequal(input$omics_plot, c("Methylation","Expression","Mutational signature"))) {
          mat <- pca_snf_all
        } else if (setequal(input$omics_plot, c("Methylation","Expression"))) {
          mat <- pca_snf_exp_meth
        } else if (setequal(input$omics_plot, c("Expression","Mutational signature"))) {
          mat <- pca_snf_exp_mut
        } else if (setequal(input$omics_plot, c("Methylation", "Mutational signature"))) {
          mat <- pca_snf_meth_mut
        }
      }
    }
    
    if (!is.null(mat)) {
      last_combined_mat(mat)
    }
    
    return(mat)
  })
  
  ### if that sample is present only in any omics in MoNETA, there will be a warning that said in which omic is present
  observeEvent(input$subset_btn,{
    if("MoNETA" %in% input$multiomics_method & length(input$omics_plot) > 1){
      
      if("Expression" %in% input$omics_plot) {
        if(!all(input$both_sample %in% rownames(pca_exp))){
          
          s1 <- input$both_sample[!input$both_sample %in% rownames(pca_exp)]
          s2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% s1]
          
          msg <- paste("This sample/s:", paste(s2, collapse = ", "), "is not present in expression layer")
          showNotification(msg, type = "warning", duration = 30)
          warning(msg)
        }
      }

      if("Methylation" %in% input$omics_plot){
        if(!all(input$both_sample %in% rownames(pca_meth_1))){
          
          s1 <- input$both_sample[!input$both_sample %in% rownames(pca_meth_1)]
          s2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% s1]
          
          msg <- paste("This sample/s:", paste(s2, collapse = ", "), "is not present in methylation layer")
          showNotification(msg, type = "warning", duration = 30)
          warning(msg)
        }
      }
      
      if("Mutational signature" %in% input$omics_plot) {
        if(!all(input$both_sample %in% rownames(combined_mat_mut))){
          
          s1 <- input$both_sample[!input$both_sample %in% rownames(combined_mat_mut)]
          s2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% s1]
          
          msg <- paste("This sample/s:", paste(s2, collapse = ", "), "is not present in mutational signature layer")
          showNotification(msg, type = "warning", duration = 30)
          warning(msg)
        }
      }

    }
  })
  
  
  a <- c('Cell lines', 'Tumors')
  b <- c('Cell lines (CL)', 'Tumors')
  choices_CL <- setNames(a,b)
  updateSelectizeInput(session, "sel_type", choices = choices_CL, selected = c("Cell lines", "Tumors"))
  
  ### get the multiomics data intergration method in the menu
  updateSelectizeInput(session, "multiomics_method", choices = c('MoNETA','MOFA','SNF'), selected = 'MoNETA')
  
  ### select the linage based on the mat that the user choose
  selected_linages <- reactiveVal(NULL)
  l <- reactiveVal(NULL)
  
  observe({
    lin <- unique(ann_multiomics_v9$lineage[ann_multiomics_v9$sampleID %in% rownames(selected_combined_mat())])
    selected_linages(lin)
    
    sel <- l()
    if (!is.null(sel) && sel %in% lin) {
      updateSelectizeInput(session, "sel_lineage", choices = lin, selected = sel, server = TRUE)
    } else {
      updateSelectizeInput(session, "sel_lineage", choices = lin, selected = character(0), server = TRUE)
    }
  })
  
  ### update the value of l with the chooise of the user
  observeEvent(input$subset_btn, {
    l(input$sel_lineage)
  })
  
  
  ### select the Query lineage/s for neighbors search
  q <- reactiveVal(NULL)
  query_linages <- reactiveVal(NULL)
  
  observe({
    lin_out <- unique(ann_multiomics_v9$lineage[ann_multiomics_v9$sampleID %in% rownames(selected_combined_mat())])
    query_linages(lin_out)
    
    sel <- q()
    if (!is.null(sel) && (any(sel %in% lin_out) || isTRUE(sel == 'All'))) {
      updateSelectizeInput(session, "lin_output", choices = c('All', lin_out), selected = sel)
    } else {
      updateSelectizeInput(session, "lin_output", choices = c('All', lin_out), selected = 'All')
    }
  })
  
  ### update the value of q with the chooise of the user
  observeEvent(input$subset_btn, {
    q(input$lin_output)
  })
  
  ### update the sample/s for each selected omics; using depmap code but show nameID
  r_choices <- reactiveVal()
  
  observe({
    choices <- list() 
    subset_1 <- ann_multiomics_v9
    
    if (input$sel_lineage == "") {
      if ("Cell lines" %in% input$sel_type) {
        cl_values <- rownames(selected_combined_mat())
        if (!is.null(cl_values)) {
          cl_values <- cl_values[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', cl_values)]
          cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v9)
          if (!is.null(cl_names)) {
            choices <- c(choices, setNames(as.list(cl_values), cl_names))
            r_choices(choices)
          }
        }
      }
      
      if ("Tumors" %in% input$sel_type) {
        tumor_values <- rownames(selected_combined_mat())
        if (!is.null(tumor_values)) {
          tumor_values <- tumor_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', tumor_values)]
          choices <- c(choices, setNames(as.list(tumor_values), tumor_values))
          r_choices(choices)
        }
      }
      
      if (all(c("Cell lines", "Tumors") %in% input$sel_type)) {
        all_values <- rownames(selected_combined_mat())
        if (!is.null(all_values)) {
          cl_values <- all_values[!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', all_values)]
          tumor_values <- all_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', all_values)]
          cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v9)
          
          if (!is.null(cl_values) && !is.null(cl_names)) {
            choices <- c(choices, setNames(as.list(cl_values), cl_names))
            r_choices(choices)
          }
          
          if (!is.null(tumor_values)) {
            choices <- c(choices, setNames(as.list(tumor_values), tumor_values))
            r_choices(choices)
          }
        }
      }
    } else {
      if (length(input$sel_lineage) >= 1) {
        subset_1 <- subset_1[ann_multiomics_v9$lineage == input$sel_lineage, ]
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
            cl_names <- get_CL_strp_names(selected_combined_mat(), ann_multiomics_v9)
            if (!is.null(cl_names)) {
              cl_filtered <- intersect(cl_values, all_values_x)
              cl_names_filtered <- intersect(cl_names, all_names_x)
              choices <- c(choices, setNames(as.list(cl_filtered), cl_names_filtered))
              r_choices(choices)
            }
          }
        }
        
        if ("Tumors" %in% input$sel_type) {
          tumor_values <- rownames(selected_combined_mat())
          if (!is.null(tumor_values)) {
            tumor_values <- tumor_values[grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', tumor_values)]
            tumor_filtered <- intersect(tumor_values, all_values_x)
            choices <- c(choices, setNames(as.list(tumor_filtered), tumor_filtered))
            r_choices(choices)
          }
        }
        
        if (all(c("Cell lines", "Tumors") %in% input$sel_type)) {
          choices <- c(choices, setNames(as.list(all_values_x), all_names_x))
          r_choices(choices)
        }
      }
    }
    
    updateSelectizeInput(session, "both_sample", choices = choices, server = TRUE)
    
    ### when user already search neighbors and he switch between the combination of input$ and then click on Plot alignment:
    ### the samples will be saved and update in the search bar!
    observeEvent(list(input$subset_btn, input$sel_lineage),{
      if(length(input$both_sample) >= 1) {
        
        selected_samples <- reactive({
          value <- ann_multiomics_v9$sampleID[ann_multiomics_v9$sampleID %in% input$both_sample]
          name <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% input$both_sample]
          choices <- setNames(as.list(value), name)
          return(choices)
        })
        
        if(!is.null(filtered_data())) {
          updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = lasso_selected_samples(), server = TRUE)
        } else {
          updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = selected_samples(), server = TRUE)
        } 
      } 
    })
    
    if(length(input$omics_plot) > 1) {
      if(!is.null(filtered_data())){
        updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = lasso_selected_samples(), server = TRUE)
      } else {
        updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = selected_samples(), server = TRUE)
      }
      
    }
    
  })
  
  ### when samples is present in one omic, then you switch omic and that samples are no more present, will appear a warning
  observeEvent(input$subset_btn,{
    if(!is.null(filtered_data())){
      if(!all(after_c() %in% rownames(selected_combined_mat()))) {
        
        g1 <- lasso_selected_samples()[!lasso_selected_samples() %in% rownames(selected_combined_mat())]
        g2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% g1]
        
        if(length(g2) > 0) {
          
          msg <- paste("There isn't this input sample/s:", paste(g2, collapse = ", "), "in this omics")
          showNotification(msg, type = "warning", duration = 30)
          warning(msg)
        }
        
        updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = lasso_selected_samples(), server = TRUE)
      }
    } else {
      if(!all(selected_samples() %in% rownames(selected_combined_mat()))) {
        
        g1 <- selected_samples()[!selected_samples() %in% rownames(selected_combined_mat())]
        g2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% g1]
        
        msg <- paste("There isn't this input sample/s:", paste(g2, collapse = ", "), "in this omics")
        showNotification(msg, type = "warning", duration = 30)
        warning(msg)
        
        updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = selected_samples(), server = TRUE)
        
      }
    }
  }, ignoreInit = TRUE)
  
  output$plot <- renderUI({
    selected_plot()
  })
  
  ### when you change omics or red_method or integration_method the piecharts will go away
  observeEvent(list(input$multiomics_method, input$omics_plot, input$reduction_method), {
    output$piechart_subtype <- renderPlot({ NULL })
    output$piechart <- renderPlot({ NULL })
    
    output$distribution_l <- renderPlotly({ NULL })
    output$distribution_s <- renderPlotly({ NULL })
  })
  
  selected_samples <- reactive({
    input$both_sample
  })
  
  #### click show in the shiny to find neighbors
  #### click show in the shiny to get both kind of piechart
  observeEvent(input$subset_btn, {
    
    ### when there is no sample/s in search bar, reolad the omics base plot
    if(is.null(input$both_sample) || length(input$both_sample) == 0) {
      
      output$plot <- renderUI({
        selected_plot()
      })
      
      output$distribution_l <- renderPlotly({ NULL })
      output$distribution_s <- renderPlotly({ NULL })
    }
    
    
    
    if(length(input$both_sample) == 1) {
      
      n <- find_neighbors(combined_mat = selected_combined_mat(), 
                          reduced_mat = selected_reduced_mat(),
                          input_sample = input$both_sample,
                          k = input$num_neighbors,
                          ann = ann_multiomics_v9,
                          type = input$df_selection_output,
                          query_lineage = input$lin_output)
      
      if(!is.null(n)) {
        
        x <- get_alignment_plot(reduced_mat = selected_reduced_mat(),
                                ann = ann_multiomics_v9,
                                dist_top_n = n, input_sample = input$both_sample)
        
        output$plot <- renderUI({
          x  
        })
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                 input_sample = input$both_sample,
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v9,
                                 type = input$df_selection_output,
                                 value = 'lineage',
                                 dist_top_n = n)
        
        piechart_subtype <- get_piechart(combined_mat = selected_combined_mat(), 
                                         input_sample = input$both_sample,
                                         k = input$num_neighbors,
                                         ann = ann_multiomics_v9,
                                         type = input$df_selection_output,
                                         value = 'subtype',
                                         dist_top_n = n)
        
        output$piechart_subtype <- renderPlot({
          piechart_subtype
        })
        
        output$piechart <- renderPlot({
          piechart
        })
        
        p <- c_distribution(dist_top_n = n, ann = ann_multiomics_v9)
        
        output$distribution_l <- renderPlotly({p$lineage_distribution})
        output$distribution_s <- renderPlotly({p$subtype_distribution})
        
      }
      
      else {}
      
    } 
    else if (length(input$both_sample) > 1) {
      selected_samples <- reactive({
        input$both_sample
      })
      
      n <- find_neighbors(combined_mat = selected_combined_mat(), 
                          reduced_mat = selected_reduced_mat(),
                          selected_samples = selected_samples(),
                          k = input$num_neighbors,
                          ann = ann_multiomics_v9,
                          type = input$df_selection_output,
                          query_lineage = input$lin_output)
      
      if(!is.null(n)) {
        
        x <- get_alignment_plot(reduced_mat = selected_reduced_mat(),
                                ann = ann_multiomics_v9,
                                dist_top_n = n, selected_samples = selected_samples())
        
        
        output$plot <- renderUI({
          x  
        })
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                 selected_samples = selected_samples(),
                                 k = input$num_neighbors,
                                 ann = ann_multiomics_v9,
                                 type = input$df_selection_output,
                                 value = 'lineage',
                                 dist_top_n = n)
        
        piechart_subtype <- get_piechart(combined_mat = selected_combined_mat(), 
                                         selected_samples = selected_samples(),
                                         k = input$num_neighbors,
                                         ann = ann_multiomics_v9,
                                         type = input$df_selection_output,
                                         value = 'subtype',
                                         dist_top_n = n)
        
        output$piechart_subtype <- renderPlot({
          piechart_subtype
        })
        
        output$piechart <- renderPlot({
          piechart
        })
        
        p <- c_distribution(dist_top_n = n, ann = ann_multiomics_v9)
        
        output$distribution_l <- renderPlotly({p$lineage_distribution})
        output$distribution_s <- renderPlotly({p$subtype_distribution})
        
      }
    }
    
    else {}
    
  })
  
  #### lasso select code
  
  filtered_data <- reactiveVal()
  after_c <- reactiveVal()
  
  observeEvent(event_data("plotly_selected"), {
    
    if(is.null(input$sel_type)) {
      showNotification("Select the model type and then load the samples from the plot")
      warning("Select the model type and then load the samples from the plot")
      return()
    }
    
    if(length(input$sel_lineage) > 1) {
      showNotification("Deselect the section: Select lineage", type = "warning")
      warning("Deselect the section: Select lineage")
    }
    
    selected_data <- event_data("plotly_selected")
    
    if (!is.null(selected_data)) {
      selected_samples <- selected_data$key
      
      x_1 <- ann_multiomics_v9 %>%
        filter(stripped_cell_line_name %in% selected_samples)
      
      c <- ann_multiomics_v9 %>%
        filter(sampleID %in% selected_samples)
      
      after_c(c)
      
      x_2 <- NULL
      
      if(all(c("Cell lines", "Tumors") %in% input$sel_type)) { 
        x_2 <- x_1 
      }
      
      else if('Tumors' %in% input$sel_type) {
        x_2 <- x_1 %>% dplyr::filter(type == 'Tumor')
      }
      
      else if('Cell lines' %in% input$sel_type) {
        x_2 <- x_1 %>% dplyr::filter(type == 'CL')
      }
      
      filtered_data(x_2$sampleID)
    }
  })
  
  lasso_selected_samples <- reactive({
    req(filtered_data()) 
    
    filtered_data()
  })
  
  observeEvent(input$load_selection, {
    if(is.null(input$sel_type)) {
      return()
    }
    
    req(lasso_selected_samples())
    req(input$sel_type)
    
    selected <- list()
    
    all_names <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% rownames(selected_combined_mat())]
    all_values <- colnames(selected_reduced_mat())
    
    sub_values <- all_values[all_values %in% lasso_selected_samples()]
    sub_names <- ann_multiomics_v9 %>% filter(sampleID %in% sub_values)
    sub_names <- sub_names$stripped_cell_line_name
    
    sub_nccle <- sub_names[!grepl("TCGA|TARGET|TH0|TH1|TH2|TH3|THR", sub_names)]
    
    sub_names_tcga <- sub_values[grepl("TCGA|TARGET|TH0|TH1|TH2|TH3|THR", sub_values)]
    
    sub_names <- c(sub_names_tcga, sub_nccle)
    
    selected <- c(selected, setNames(as.list(sub_values), sub_names))
    
    updateSelectizeInput(session, "both_sample",
                         choices = r_choices(),
                         selected = selected,
                         server = TRUE)
  })
  
  ### remove all the sample form the query bar  
  observeEvent(input$rm, {
    filtered_data(NULL)
    updateSelectizeInput(session, "both_sample", choices = r_choices(), server = TRUE)
  })
  
  ### debug
  observeEvent(input$subset_btn, {
    if(is.null(input$omics_plot)) {
      showNotification("Select at least one omic")
      warning("Select at least one omic")
      return()
    }
  })
  
  ### allow multiomics research when is selected almost two omics
  observe({
    if(length(input$omics_plot) <= 1) {
      shinyjs::disable(id = "multiomics_method")
    }
    else {
      shinyjs::enable(id = "multiomics_method")
    }
  })
  
  
  
}

shiny::shinyApp(ui, server)

