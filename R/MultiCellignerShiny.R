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

MultiCellignerShiny <- function() {shiny::addResourcePath("static", system.file("www", package = "MultiCelligner"))
  ;shiny::shinyApp(ui, server)}


################################## UI ###################################
ui <- fluidPage(              
  
  shinyjs::useShinyjs(),
  
  shiny::sidebarLayout(
    ################# side bar panel #################
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
                          shiny::selectInput('omics_plot', '(Multi-)Omics alignment:',
                                      choices = c('Methylation',
                                                  "Mutational signature",
                                                  "Expression"),
                                      multiple = TRUE
                          )),

        )),
      
      shiny::fluidRow(
        shiny::column(6, 
                      div(style = "display: flex; flex-direction: shiny::column; align-items: center; justify-content: center; height: 100%;",
                          shiny::selectInput('reduction_method', "Reduction:",
                                      choices = c('UMAP', "tSNE"),
                                      selected = 'UMAP',
                                      width = "150px"))
        ),
        
        shiny::column(6,
                      div(style = "display: flex; flex-direction: shiny::column; align-items: center; justify-content: center; height: 100%;",
                          shiny::selectInput('multiomics_method', 'Integration method:',
                                      choices = c('MoNETA','MOFA','SNF'),
                                      selected = 'MoNETA'))
        )),
      
      hr(),
      
      tags$p("Find neighbors", style = "text-align: center; font-size: 18px; font-weight: bold;"),
      
      shiny::fluidRow(
        shiny::column(6,
                      shiny::selectizeInput("sel_type",
                                     'Select model type',
                                     choices = NULL, 
                                     multiple = TRUE)),
        shiny::column(6,
                      shiny::selectizeInput('sel_lineage',
                                     'Select lineage',
                                     choices = c("", sort(unique(ann_multiomics_v9$lineage))),
                                     selected = NULL,
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
      
      shiny::selectizeInput("both_sample", 
                     "Select one or more reference tumor/cell line:", 
                     choices = NULL, 
                     multiple = TRUE),
      
      hr(),
      
      shiny::fluidRow(
        shiny::column(6,
                      tags$strong("or load from map selection:")),
        shiny::column(2,
                      shiny::actionButton("load_selection", "Load", style = "text-align: center;") %>%
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
                                 shiny::actionButton("rm", "Clear Selection", style = "text-align: center;")) %>%
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
                      shiny::checkboxGroupInput("df_selection_output", 
                                         "Search among the closest:", 
                                         choices = c("Cell lines", "Tumors"), 
                                         selected = 'Tumors', inline = TRUE)),
        shiny::column(5,
                      div("Number of neighbors:", style = "text-align: left; font-weight: bold;"),
                      shiny::numericInput("num_neighbors", NULL, value = 25, min = 1, width = "70px"),
        )),
      
      shiny::fluidRow(
        shiny::column(12,
                      shiny::selectizeInput("lin_output", 
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
            title: "Click on Plot Alignment without sample/s in the search bar to get the basic plot without neighbors",
            placement: "right",
            trigger: "hover",
            html: true
          });
        });
      ')),
      
      div(style = "text-align: center;",
          shiny::actionButton("subset_btn", "Plot alignment", 
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
      
      shiny::tabsetPanel(
        type = "tabs",
        selected = 'Neighbors lineages',
        
        hr(),
        
        shiny::tabPanel("Neighbors lineages", 
                 shiny::plotOutput("piechart",height = "430px") , 
                 
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
        
        shiny::tabPanel("Neighbors subtypes", 
                 shiny::plotOutput("piechart_subtype",height = "430px"),
        ),
        
      )),
    
    ################# main panel ################# 
    shiny::mainPanel(
      tags$div(
        tags$p(h3(shiny::textOutput("omics_name"),class = "text-center"), style = "text-align: center; font-size: 18px; font-weight: bold;"),
      ),
      
      shinycssloaders::withSpinner(shiny::uiOutput("plot"), type = 6, color = "#4FC3F7", color.background = "white", ), 
      width = 9,
      
      hr(),
      
      shiny::tabsetPanel(
        type = "tabs",
        selected = 'Lineage distribution',
        
        shiny::tabPanel("Lineage distribution",          
                        plotly::plotlyOutput("distribution_l", width = '90%')
        ),
        
        shiny::tabPanel("Subtype distribution", 
                        plotly::plotlyOutput("distribution_s", width = "90%")
        ))
      
    )
  )
)                                

################################## SERVER ###################################
server <- function(input, output, session) { 
  
  
  output$omics_name <- shiny::eventReactive(input$subset_btn,{ 
    if(length(input$omics_plot) == 1) 
      paste(input$reduction_method,paste(input$omics_plot, collapse=" "))
    else paste(input$reduction_method,paste(input$omics_plot, collapse=" "), input$multiomics_method)
  })
  
  ##################### Load low-dimensional combined matrix ###################
  selected_reduced_mat <- shiny::reactive({
    
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
  
  
  
  ##################### Nearest neighbor search: retrieve combined mat ###################
  last_combined_mat <- shiny::reactiveVal(NULL)
  selected_combined_mat <- shiny::reactive({
    
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
  
  ######## Warning: not perfect overlap between omics for MoNETA ###############
  ### if that sample is present only in some omics in MoNETA, there will be a warning that said in which omics is present
  shiny::observeEvent(input$subset_btn,{
    if(input$multiomics_method == "MoNETA" & length(input$omics_plot) > 1){
      
      for (omics in input$omics_plot){
        pca_mat <- switch(omics,
                          "Expression" = pca_exp, 
                          "Methylation" = pca_meth_1,
                          "Mutational signature" = combined_mat_mut)
        
        if(!all(input$both_sample %in% rownames(pca_mat))){
          
          s1 <- input$both_sample[!input$both_sample %in% rownames(pca_mat)]
          s2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% s1]
          
          msg <- paste("This sample/s:", paste(s2, collapse = ", "), "is not present in", omics, "layer")
          shiny::showNotification(msg, type = "warning", duration = 30)
          warning(msg)
        }
      }
    }
  })
  
  ######## Side-bar options ###############
  a <- c('Cell lines', 'Tumors')
  b <- c('Cell lines (CL)', 'Tumors')
  choices_CL <- setNames(a,b)
  shiny::updateSelectizeInput(session, "sel_type", choices = choices_CL, selected = c("Cell lines", "Tumors"))
  
  
  ### update the value of l with the choice of the user
  l <-  shiny::reactiveVal(NULL)
  shiny::observeEvent(input$subset_btn, {
    l(input$sel_lineage)
  })
  
  
  ### select the Query lineage/s for neighbors search
  q <-  shiny::reactiveVal(NULL)
  query_linages <- shiny::reactiveVal(NULL)
  
  shiny::observe({
    lin_out <- unique(ann_multiomics_v9$lineage[ann_multiomics_v9$sampleID %in% rownames(selected_combined_mat())])
    query_linages(lin_out)
    
    sel <- q()
    if (!is.null(sel) && (any(sel %in% lin_out) || isTRUE(sel == 'All'))) {
      shiny::updateSelectizeInput(session, "lin_output", choices = c('All', lin_out), selected = sel)
    } else {
      shiny::updateSelectizeInput(session, "lin_output", choices = c('All', lin_out), selected = 'All')
    }
  })
  
  ### update the value of q with the choice of the user
  shiny::observeEvent(input$subset_btn, {
    q(input$lin_output)
  })
  
  
  ################# Nearest neighbor search: Select samples ####################
  ### update the sample/s for each selected omics; using depmap code but show nameID
  ###### Define the CHOICES of the SelectizeInput "both_sample" using the rownames of the selected combined mat.
  ## these choices are not controlled by a botton, but are automatically updated when changing the selected type and/or lineage.
  ## Here, no selection is defined.
  r_choices <- shiny::reactive({  
    selected_type <- input$sel_type
    selected_lineage <- if_else(input$sel_lineage == "", NA, input$sel_lineage)
    
    point_ids <- ann_multiomics_v9 %>% 
      dplyr::select(sampleID, stripped_cell_line_name, type, lineage) %>% 
      dplyr:: mutate(type2 = case_when(type == 'Tumor' ~ "Tumors", 
                                       type ==  'CL' ~ "Cell lines")  
      ) %>% 
      dplyr::filter(type2 %in% selected_type)
    
    if (!is.na(selected_lineage)) {
      point_ids <- point_ids %>% 
        dplyr::filter(lineage %in% selected_lineage)
    }
    
    choices <- point_ids %>% select(stripped_cell_line_name, sampleID) %>% tibble::deframe() %>% as.list()
    return(choices)
  })
  
  
  shiny::observeEvent(c(input$sel_type, input$sel_lineage),{
    f_selected_samples <- selected_samples()[selected_samples() %in% r_choices()]
    shiny::updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = f_selected_samples, server = TRUE)
  })
  
  
  
  ###### Define the SELECTION of the SelectizeInput "both_sample", using the typed sampleID ("selected_samples") 
  ## or the highlighted points on the plot ("lasso_selected_samples")
  
  selected_samples <- shiny::reactive({
    value <- ann_multiomics_v9$sampleID[ann_multiomics_v9$sampleID %in% input$both_sample]
    name <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% input$both_sample]
    choices <- setNames(as.list(value), name)
    return(choices)
  })
  
  
  ### When changing omics selection, previous selected samples might not be present anymore. Hence, a Warning message will appear.
  shiny::observeEvent(input$subset_btn,{
    if(!is.null(selected_samples())){ 
      shiny::updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = selected_samples(), server = TRUE)
      
      if(!all(selected_samples() %in% rownames(selected_combined_mat()))) {
        
        g1 <- selected_samples()[!selected_samples() %in% rownames(selected_combined_mat())]
        g2 <- ann_multiomics_v9$stripped_cell_line_name[ann_multiomics_v9$sampleID %in% g1]
        
        msg <- paste("There isn't this input sample/s:", paste(g2, collapse = ", "), "in this omics")
        showNotification(msg, type = "warning", duration = 30)
        warning(msg)
        
      }
    }
  }, ignoreInit = TRUE)
  
  
  
  ############################### Plot: visualization ##############################
  
  selected_plot <- shiny::eventReactive(input$subset_btn, {
    if (!is.null(selected_reduced_mat())){
      print('Generate plot ')
      get_alignment_plot(reduced_mat = selected_reduced_mat(),
                         ann = ann_multiomics_v9, 
                         annot_value = 'lineage')
    }else{
      return(NULL)
    }
    
  })
  
  #### This is important to not show the loading symbol whrn user opens the app
  output$plot <-  shiny::renderUI({
    selected_plot()
  })
  
  
  #### click show in the shiny to find neighbors
  #### click show in the shiny to get both kind of piechart
  shiny::observeEvent(input$subset_btn, {
    
    ### when there is no sample/s in search bar, reolad the omics base plot
    if(is.null(input$both_sample) || length(input$both_sample) == 0) {
      
      output$plot <-  shiny::renderUI({
        selected_plot()
      })
      
      output$distribution_l <- plotly::renderPlotly({ NULL })
      output$distribution_s <- plotly::renderPlotly({ NULL })
      
    } else {
      
      n <- find_neighbors(combined_mat = selected_combined_mat(), 
                          reduced_mat = selected_reduced_mat(),
                          input_sample = input$both_sample,
                          k = input$num_neighbors,
                          ann = ann_multiomics_v9,
                          query_type = input$df_selection_output,
                          query_lineage = input$lin_output)
      
      if(!is.null(n)) {
        
        x <- get_alignment_plot(reduced_mat = selected_reduced_mat(),
                                ann = ann_multiomics_v9,
                                dist_top_n = n, 
                                input_sample = input$both_sample)
        
        output$plot <-  shiny::renderUI({
          x  
        })
        
        piechart <- get_piechart(combined_mat = selected_combined_mat(), 
                                 input_sample = input$both_sample,
                                 ann = ann_multiomics_v9,
                                 type = input$df_selection_output,
                                 value = 'lineage',
                                 dist_top_n = n)
        
        piechart_subtype <- get_piechart(combined_mat = selected_combined_mat(), 
                                         input_sample = input$both_sample,
                                         ann = ann_multiomics_v9,
                                         type = input$df_selection_output,
                                         value = 'subtype',
                                         dist_top_n = n)
        
        output$piechart <-  shiny::renderPlot({
          piechart
        })
        
        output$piechart_subtype <-  shiny::renderPlot({
          piechart_subtype
        })
        
        
        p <- c_distribution(dist_top_n = n, ann = ann_multiomics_v9)
        
        output$distribution_l <- plotly::renderPlotly({p$lineage_distribution})
        output$distribution_s <- plotly::renderPlotly({p$subtype_distribution})
        
      }
    } 
  })
  
  
  ### when you change omics or red_method or integration_method the piecharts will go away
  shiny::observeEvent(list(input$multiomics_method, input$omics_plot, input$reduction_method), {
    output$piechart_subtype <- shiny::renderPlot({ NULL })
    output$piechart <- shiny::renderPlot({ NULL })
    
    output$distribution_l <- plotly::renderPlotly({ NULL })
    output$distribution_s <- plotly::renderPlotly({ NULL })
  })
  
  
  
  ############################# Plot: lasso select  ############################# 
  
  filtered_data <- shiny::reactiveVal()
  
  shiny::observeEvent(plotly::event_data("plotly_selected"), {
    
    if(is.null(input$sel_type)) {
      shiny::showNotification("Select the model type and then load the samples from the plot")
      warning("Select the model type and then load the samples from the plot")
      return()
    }
    
    if(length(input$sel_lineage) == 1) {
      shiny::showNotification(paste("", input$sel_lineage), type = "warning")
      warning("")
    }
    
    selected_data <- plotly::event_data("plotly_selected")
    
    if (!is.null(selected_data)) {
      selected_samples <- selected_data$key ### these correspond to the stripped_cell_line_name
      
      x_1 <- ann_multiomics_v9 %>%
        dplyr::filter(stripped_cell_line_name %in% selected_samples) # select cl and tumors
      
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
  
  ##### "lasso_selected_samples" is equal to "filtered_data", but it is not affected by the "rm" button
  lasso_selected_samples <- shiny::reactive({
    req(filtered_data()) 
    filtered_data()
  })
  
  
  ### load selected sample from the plot
  shiny::observeEvent(input$load_selection, {
    if(is.null(input$sel_type)) {
      return()
    }
    
    req(lasso_selected_samples())
    req(input$sel_type)
    
    # commento perchè se cambio selected_reduced_mat e riclicco load mi restituisce l'intersezione
    # noi invece vogliamo continuare a vedere tutte le query come selected e usare il warning per stampare quelle assenti 
    # all_values <- colnames(selected_reduced_mat())
    # sub_values <- all_values[all_values %in% lasso_selected_samples()]
    sub_values <- lasso_selected_samples()
    
    selected <- ann_multiomics_v9 %>% filter(sampleID %in% sub_values) %>%
      select(stripped_cell_line_name, sampleID) %>% tibble::deframe() %>% as.list()
    
    shiny::updateSelectizeInput(
      session, "both_sample",
      choices = r_choices(),
      selected = selected,
      server = TRUE)
  })
  
  
  ### remove all the sample form the query bar  
  shiny::observeEvent(input$rm, {
    # filtered_data(NULL)
    shiny::updateSelectizeInput(session, "both_sample", choices = r_choices(), selected = NULL, server = TRUE)
  })
  
  
  ############################# MISCELLANEOUS   #############################  
  
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
