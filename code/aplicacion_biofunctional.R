# Define UI
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = div(
    style = "font-size: 14px; font-weight: bold; font-family: Georgia;",
    "BIOFunctional"
  )),
  dashboardSidebar(
    sidebarMenu(
      menuItem("HOME", tabName = "inicio", icon = icon("home")),
      menuItem("GANGO Data Explorer", tabName = "gango"),
      menuItem("KEGG", tabName = "kegg", icon = icon("chart-bar"),
               menuSubItem("Filter", tabName = "opcion3_kegg"),
               menuSubItem("Functional Analysis", tabName = "opcion1_kegg"),
               menuSubItem("Network Analysis", tabName = "opcion2_kegg"),
               menuSubItem("HeatMap",tabName = "Heatmap")
      ),
      menuItem("Gene Ontologies", tabName = "gene_ontologies", icon = icon("dna"),
               menuSubItem("Filter", tabName = "opcion3_gene_ontologies"),
               menuSubItem("Functional Analysis", tabName = "opcion1_gene_ontologies"),
               menuSubItem("Network Analysis", tabName = "opcion2_gene_ontologies"),
               menuSubItem("BarPlot",tabName = "Heatmap_go")
      ),
      menuItem("HELP", tabName = "contact", icon = icon("clipboard-user"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "gango",
              # Parámetros de consulta
              fluidRow(
                box(
                  title = "Query Parameters", status = "primary", solidHeader = TRUE, width = 12,
                  textInput("taxon", "Enter Taxon", value = "Homo sapiens, Mus musculus, Saccharomyces cerevisiae", placeholder = "Enter taxons separated by commas"),
                  textInput("gene", "Enter Gene", value = "BRCA1, Tp53, CDC28", placeholder = "Enter genes separated by commas"),
                  textInput("group", "Enter Group", value = "Mammals, Mammals, Fungi", placeholder = "Enter groups separated by commas"),
                  fileInput("file1", "Choose CSV File", accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                  actionButton("submit", "Submit")
                )
              ),
              
              # Resultados divididos en pestañas GO / KEGG
              fluidRow(
                box(
                  title = "Query Results", status = "primary", solidHeader = TRUE, width = 12,
                  shinydashboard::tabBox(
                    id = "results_tabs", width = 12,
                    tabPanel("GO",
                             DT::dataTableOutput("goResults"),
                             br(),
                             downloadButton("downloadGO", "Download GO Results", class = "btn-primary")
                    ),
                    tabPanel("KEGG",
                             DT::dataTableOutput("keggResults"),
                             br(),
                             downloadButton("downloadKEGG", "Download KEGG Results", class = "btn-primary")
                    )
                  )
                )
              ),
              
              # Volcano plots divididos en pestañas GO / KEGG
              fluidRow(
                box(
                  title = "Volcano Plot", status = "primary", solidHeader = TRUE, width = 12,
                  shinydashboard::tabBox(
                    id = "volcano_tabs", width = 12,
                    tabPanel("GO Volcano",
                             plotOutput("volcanoGO", height = "500px"),
                             br(),
                             downloadButton("downloadVolcanoGO", "Download GO Volcano", class = "btn-danger")
                    ),
                    tabPanel("KEGG Volcano",
                             plotOutput("volcanoKEGG", height = "500px"),
                             br(),
                             downloadButton("downloadVolcanoKEGG", "Download KEGG Volcano", class = "btn-danger")
                    )
                  )
                )
              )
      )
      ,
      tabItem(tabName = "inicio",
              fluidRow(
                box(
                  title = "A comprehensive software for interpretation and visualization of functional analysis of Gene Ontologies and KEGG Pathways",
                  width = 12,
                  solidHeader = TRUE,
                  status = "info",
                  style = "font-size: 14px; font-family: Times New Roman;",
                  p("A comprehensive app designed for the interpretation and visualization of the functional analysis related to KEGG pathways and gene ontologies gives researchers and specialists a tool to get detailed functional information about their data, specifically going deep into biological pathways and gene functions information. By using a variety of techniques and libraries, such as Shiny, htrr, dplyr, tibble, and rvest, we have developed an application that provides a well-designed user-oriented interface with all the facilities to assess their data and start analyzing it directly from scratch through a few steps."),
                  br(),
                  p("The software allows an exhaustive exploration of KEGG pathways and Gene Ontologies, facilitating the analysis of complex biological processes.To achieve this, functions described in the scripts integrate data manipulation methods and web scraping techniques to extract the necessary information from online official databases, Kyoto Encyclopedia of Genes and Genomes (KEGG) and QuickGo. Furthermore, those functions are computed by parallel processing, resulting in efficient petitions to the database servers and allowing the user to get quick results from a large dataset."),
                  br(),
                  p("A fundamental feature in the app is the capability, by the techniques explained above, to obtain ancestral information for KEGG pathways and gene ontologies, making it easier to understand their hierarchy and how each of the samples in a dataset is classified through it.This offers users a way to study the dataset at different levels of taxonomy directly from the raw data.Additionally, the ability to create interactive networks is implemented, aiming to represent all the experimental data to see the relationships between the groups and the ontologies, without disregarding the classification created. This is the main tool to understand the meaning of the relations that will be seen around the system displayed."),
                  br(),
                  p("As a result of all these attributes, the software represents a key tool for the analyst involved in the study of biological pathways, providing an intuitive interface with advanced data processing techniques, allowing researchers to puzzle out the intricacy of the biological functions and obtain insights into the relationships between genes or molecular components.")
                ),
                div(
                  img(src = "https://biysc.org/sites/default/files/ub_facultat_biologia.png", height = 300, width = 600)
                )
              )
      ),
      tabItem(tabName = "kegg",
              fluidRow(
                box(title = "KEGG",
                    width = 12,
                    solidHeader = TRUE,
                    status = "primary"
                )
              )
      ),
      tabItem(tabName = "gene_ontologies",
              fluidRow(
                box(title = "Gene Ontologies",
                    width = 12,
                    solidHeader = TRUE,
                    status = "warning"
                )
              )
      ),
      tabItem(tabName = "contact",
              fluidRow(
                box(title = "Contact Information", width = 12, status = "info",
                    "Developed by Alejandro Rodríguez & Antonio Monleón. Section of Statistics. Department of Genetics, Microbiology and Statistics. UB. For any inquiries or support, please contact us at:",
                    br(),
                    br(),
                    "Alex: alejandro.rodriguez@alum.esci.upf.edu",
                    br(),
                    "Toni: amonleong@ub.edu"
                )
              )
      ),
      tabItem(
        tabName = "opcion3_kegg",
        fluidRow(
          # Caja para cargar el archivo y seleccionar filtros
          box(
            title = "Filter Data",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            fileInput("file_filter_kegg", "Upload CSV File"),
            uiOutput("column_filter_selection"),  # Para seleccionar filtros dinámicamente
            br(),
            downloadButton("downloadFile_filter_kegg", "Download Processed File", class = "btn-primary")
          ),
          
          # Caja para mostrar la tabla de datos filtrados
          box(
            title = "Filtered Data",
            width = 12,
            solidHeader = TRUE,
            status = "info",
            withSpinner(dataTableOutput("table_filter_kegg")),
            br(),
            actionButton("clear_filter_kegg", "Clear Data", icon = icon("trash"), class = "btn-danger")
          )
        )
      )
      ,
      tabItem(tabName = "opcion1_kegg",
              fluidRow(
                box(title = "Functional Analysis",
                    width = 12,
                    solidHeader = TRUE,
                    status = "success",
                    fileInput("file_kegg", "Upload File"),
                    downloadButton("downloadFile_kegg", "Download Processed File"),
                    br(),
                    br(),
                    br(),
                    withSpinner(
                      dataTableOutput("table_kegg")
                    ),
                    br(),
                    actionButton("clear_kegg", "Clear Data", icon = icon("trash"))
                )
              )
      ),
      tabItem(tabName = "opcion2_kegg",
              fluidRow(
                box(title = "Network Analysis Filters", width = 3, status = "danger",
                    fileInput("file_kegg_network", "Upload File"),
                    uiOutput("dynamic_filters"),                   
                    selectInput("group_kegg", "Group:", choices = NULL, multiple = FALSE),
                    selectInput("group_by_kegg", "Group by:", choices = NULL, multiple = FALSE),
                    downloadButton("download_network_kegg", "Download Network"),
                    br(),
                    actionButton("clear_kegg_network", "Clear Data", icon = icon("trash")),
                    uiOutput("legend")
                ),
                box(title = "Network", width = 9, status = "danger",
                    withSpinner(visNetworkOutput("network_kegg", width = "100%", height = "800px")),
                    box(title = tagList("AI conclusions prompt", icon("question-circle", id = "helpIcon")), width = 12, status = "info",
                        textOutput("text_kegg"),
                        br(),
                        downloadButton("download_prompt_kegg", "Download Prompt")
                    )
                )
              )
      ),
      tabItem(
        tabName = "opcion3_gene_ontologies",
        fluidRow(
          # Caja para cargar el archivo y seleccionar filtros
          box(
            title = "Filter Data",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            fileInput("file_filter_gene_ontologies", "Upload CSV File"),
            uiOutput("column_selectors_go"),  # Para seleccionar filtros dinámicamente
            br(),
            downloadButton("downloadFile_filter_go", "Download Processed File", class = "btn-primary")
          ),
          
          # Caja para mostrar la tabla de datos filtrados
          box(
            title = "Filtered Data",
            width = 12,
            solidHeader = TRUE,
            status = "info",
            withSpinner(dataTableOutput("table_filter_go")),
            br(),
            actionButton("clear_filter_go", "Clear Data", icon = icon("trash"), class = "btn-danger")
          )
        )
      ),
      tabItem(tabName = "opcion1_gene_ontologies",
              fluidRow(
                box(title = "Functional Analysis",
                    width = 12,
                    solidHeader = TRUE,
                    status = "success",
                    fileInput("file_gene_ontologies", "Upload File"),
                    downloadButton("downloadFile_gene_ontologies", "Download Processed File"),
                    br(),
                    br(),
                    br(),
                    withSpinner(
                      dataTableOutput("table_gene_ontologies")
                    ),
                    br(),
                    actionButton("clear_go", "Clear Data", icon = icon("trash"))
                )
              )
      ),
      tabItem(tabName = "Heatmap",
              fluidRow(
                box(title = "Heatmap Analysis",
                    width = 12,
                    solidHeader = TRUE,
                    status = "success",
                    fileInput("file_heatmap", "Upload KEGG Data"),
                    uiOutput("dynamic_filters_heatmap"),
                    selectInput("group_by_heatmap", "Group By:", choices = c("Domain", "Subdomain", "Relation")),
                    selectInput("group_heatmap", "Group", choices = NULL),
                    downloadButton("download_heatmap", "Download Heatmap")
                ),
                box(title = "Heatmap Visualization",
                    width = 12,
                    solidHeader = TRUE,
                    status = "info",
                    withSpinner(plotlyOutput("heatmap_plot", height = "600px"))
                )
              )
      ),
      tabItem(
        tabName = "Heatmap_go",
        fluidRow(
          box(
            title = "GO Data Analysis",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            fileInput("file_heatmap_go", "Upload GO Data"),
            
            # Filtros obligatorios
            selectInput("ont_description_heatmap_go", "Ontology Description", choices = NULL),  # Filtro obligatorio
            selectInput("group_heatmap_go", "Group", choices = NULL),  # Filtro obligatorio
            
            # Filtros dinámicos generados
            uiOutput("dynamic_filters_heatmap_go"),  # Filtros dinámicos
            
            downloadButton("download_barplot_up_go", "Download Up Regulation Plot"),
            downloadButton("download_barplot_down_go", "Download Down Regulation Plot")
          ),
          box(
            title = "Up Regulation Visualization",
            width = 12,
            solidHeader = TRUE,
            status = "info",
            withSpinner(plotlyOutput("barplot_up_go", height = "600px"))
          ),
          box(
            title = "Down Regulation Visualization",
            width = 12,
            solidHeader = TRUE,
            status = "warning",
            withSpinner(plotlyOutput("barplot_down_go", height = "600px"))
          )
        )
      )
      ,
      tabItem(tabName = "Volcano",
              fluidRow(
                box(title = "Volcano Analysis",
                    width = 12,
                    solidHeader = TRUE,
                    status = "success",
                    fileInput("file_volcano", "Upload KEGG Data"),
                    selectInput("sample_volcano", "Sample", choices = NULL),
                    selectInput("disease_volcano", "Disease", choices = NULL),
                    selectInput("group_by_volcano", "Group By:", choices = c("Domain", "Subdomain", "Relation")),
                    selectInput("ontology_volcano", "Ontology", choices = NULL),
                    downloadButton("download_volcano", "Download Volcano")
                ),
                box(title = "Volcano Visualization",
                    width = 12,
                    solidHeader = TRUE,
                    status = "info",
                    withSpinner(plotlyOutput("volcano_plot", height = "600px"))
                )
              )
      ),
      tabItem(tabName = "Volcano_go",
              fluidRow(
                box(title = "Volcano Analysis",
                    width = 12,
                    solidHeader = TRUE,
                    status = "success",
                    fileInput("file_volcano_go", "Upload KEGG Data"),
                    selectInput("sample_volcano_go", "Sample", choices = NULL),
                    selectInput("disease_volcano_go", "Disease", choices = NULL),
                    selectInput("ont_description_volcano_go", "Ontology Description", choices = NULL), 
                    downloadButton("download_volcano_go", "Download Volcano")
                ),
                box(title = "Volcano Visualization",
                    width = 12,
                    solidHeader = TRUE,
                    status = "info",
                    withSpinner(plotlyOutput("volcano_plot_go", height = "600px"))
                )
              )
      ),
      tabItem(tabName = "opcion2_gene_ontologies",
              fluidRow(
                box(title = "Network Analysis Filters", width = 3, status = "danger",
                    fileInput("file_go_network", "Upload File"),
                    selectInput("ont_description", "Ontology Description", choices = NULL, multiple = FALSE),
                    selectInput("group_go", "Group:", choices = NULL, multiple = FALSE),
                    uiOutput("dynamic_filters_go"),  # Filtros dinámicos
                    downloadButton("download_network_go", "Download Network"),
                    br(),
                    actionButton("clear_go_network", "Clear Data", icon = icon("trash")),
                    uiOutput("legend_go")
                ),
                box(title = "Network", width = 9, status = "danger",
                    withSpinner(visNetworkOutput("network_go", width = "100%", height = "800px")),
                    box(title = tagList("AI conclusions prompt", icon("question-circle", id = "helpIcon2")), width = 12, status = "info",
                        textOutput("text_go"),
                        br(),
                        downloadButton("download_prompt_go", "Download Prompt")
                    )
                )
              )
      )
      
    )
  )
)




#Define server logic
server <- function(input, output, session) {
  
  observeEvent(input$submit, {
    # --- 1. Leer inputs ---
    if (!is.null(input$file1)) {
      file_data <- read.csv(input$file1$datapath)
      taxon <- unique(trimws(file_data$taxon))
      gene  <- unique(trimws(file_data$gene))
      group <- unique(trimws(file_data$group))
    } else {
      req(input$taxon, input$gene, input$group)
      taxon <- trimws(strsplit(input$taxon, ",")[[1]])
      gene  <- trimws(strsplit(input$gene,  ",")[[1]])
      group <- trimws(strsplit(input$group, ",")[[1]])
    }
    if (length(taxon)!=length(gene) || length(taxon)!=length(group)) {
      showModal(modalDialog(
        title = "Input Error",
        "The number of taxons, genes, and groups must match.",
        easyClose = TRUE, footer = NULL
      ))
      return()
    }
    
    # --- 2. Ejecutar GANGO ---
    result_data <- GANGO(taxon, gene, group)
    if (is.null(result_data)) {
      showModal(modalDialog(
        title = "No Data Found",
        "GANGO returned no results. Please check your inputs.",
        easyClose = TRUE, footer = NULL
      ))
      return()
    }
    
    # --- 3. Extraer significant_results para GO y KEGG ---
    go_sig   <- result_data$GO_results$significant_results
    kegg_sig <- result_data$KEGG_results$significant_results
    
    # === GO ===
    # Depuración
    cat("=== GO Significant Results ===\n")
    cat("Rows:", nrow(go_sig), "\n")
    print(head(go_sig))
    
    # Tabla GO
    output$goResults <- DT::renderDataTable({
      clean_go <- go_sig %>%
        rename_with(~ make.names(.x, TRUE)) %>%
        filter(!is.na(ONTOLOGY) & !is.na(EA_VALUE) & !is.na(pvalue) & !is.na(FDR))
      DT::datatable(clean_go, options = list(pageLength = 10))
    })
    # Download GO table
    output$downloadGO <- downloadHandler(
      filename = function() paste0("GANGO_GO_significant_", Sys.Date(), ".csv"),
      content = function(file) write.csv(go_sig, file, row.names = FALSE)
    )
    # Volcano GO
    output$volcanoGO <- renderPlot({
      data <- go_sig
      ggplot(data, aes(x = EA_VALUE, y = -log10(FDR), color = UP_DOWN, label = ONTOLOGY)) +
        geom_point(aes(size = -log10(FDR)), alpha = .8, shape = 21, fill = "white") +
        scale_color_manual(values = c("UP"="red","DOWN"="blue","NEUTRAL"="grey")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
        geom_vline(xintercept = 0) +
        geom_text(data = subset(data, -log10(FDR)>3), size = 3) +
        labs(title="Volcano Plot - GO (Significant)", x="Enrichment Score", y="-log10(FDR)") +
        theme_minimal(base_size=14)
    })
    # Download Volcano GO
    output$downloadVolcanoGO <- downloadHandler(
      filename = function() "Volcano_GO_significant.png",
      content = function(file) {
        png(file, width=1200, height=900, res=150)
        print(ggplot(go_sig, aes(x = EA_VALUE, y = -log10(FDR), color = UP_DOWN, label = ONTOLOGY)) +
                geom_point(aes(size = -log10(FDR)), alpha = .8, shape = 21, fill = "white") +
                scale_color_manual(values = c("UP"="red","DOWN"="blue","NEUTRAL"="grey")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
                geom_vline(xintercept = 0) +
                geom_text(data = subset(go_sig, -log10(FDR)>3), size = 3) +
                labs(title="Volcano Plot - GO (Significant)", x="Enrichment Score", y="-log10(FDR)") +
                theme_minimal(base_size=14))
        dev.off()
      }
    )
    
    # === KEGG ===
    # Depuración
    cat("=== KEGG Significant Results ===\n")
    cat("Rows:", nrow(kegg_sig), "\n")
    print(head(kegg_sig))
    
    # Tabla KEGG
    output$keggResults <- DT::renderDataTable({
      req(kegg_sig)
      
      kegg_clean <- kegg_sig %>%
        rename_with(~make.names(.x, TRUE)) %>%
        filter(
          !is.na(KEGG_Pathway) & 
            !is.na(EA_VALUE) & 
            !is.na(pvalue) & 
            !is.na(FDR)
        ) %>%
        mutate(Row = row_number()) %>%
        select(Row, everything()) %>%
        rename(
          ONTOLOGY = KEGG_Pathway,
          ONT_DESCRIPTION = KEGG_Name
        )
      
      DT::datatable(
        kegg_clean,
        rownames = FALSE,
        options = list(pageLength = 10),
        callback = DT::JS("table.columns(0).header().to$().text('');")
      )
    })
    
    
    # Download KEGG table
    output$downloadKEGG <- downloadHandler(
      filename = function() paste0("GANGO_KEGG_significant_", Sys.Date(), ".csv"),
      content = function(file) {
        kegg_clean <- kegg_sig %>%
          rename_with(~make.names(.x, TRUE)) %>%
          filter(
            !is.na(KEGG_Pathway) & 
              !is.na(EA_VALUE) & 
              !is.na(pvalue) & 
              !is.na(FDR)
          ) %>%
          rename(
            ONTOLOGY = KEGG_Pathway,
            ONT_DESCRIPTION = KEGG_Name
          )
        write.csv(kegg_clean, file, row.names = FALSE)
      }
    )
    
    
    # Volcano KEGG
    output$volcanoKEGG <- renderPlot({
      data <- kegg_sig
      ggplot(data, aes(x = EA_VALUE, y = -log10(FDR), color = UP_DOWN, label = KEGG_Pathway)) +
        geom_point(aes(size = -log10(FDR)), alpha = .8, shape = 21, fill = "white") +
        scale_color_manual(values = c("UP"="red","DOWN"="blue","NEUTRAL"="grey")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
        geom_vline(xintercept = 0) +
        geom_text(data = subset(data, -log10(FDR)>3), size = 3) +
        labs(title="Volcano Plot - KEGG (Significant)", x="Enrichment Score (EA_VALUE)", y="-log10(FDR)") +
        theme_minimal(base_size=14)
    })
    # Download Volcano KEGG
    output$downloadVolcanoKEGG <- downloadHandler(
      filename = function() "Volcano_KEGG_significant.png",
      content = function(file) {
        png(file, width=1200, height=900, res=150)
        print(ggplot(kegg_sig, aes(x = EA_VALUE, y = -log10(FDR), color = UP_DOWN, label = KEGG_Pathway)) +
                geom_point(aes(size = -log10(FDR)), alpha = .8, shape = 21, fill = "white") +
                scale_color_manual(values = c("UP"="red","DOWN"="blue","NEUTRAL"="grey")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
                geom_vline(xintercept = 0) +
                geom_text(data = subset(kegg_sig, -log10(FDR)>3), size = 3) +
                labs(title="Volcano Plot - KEGG (Significant)", x="Enrichment Score (EA_VALUE)", y="-log10(FDR)") +
                theme_minimal(base_size=14))
        dev.off()
      }
    )
  })
  
  
  # Render popovers
  output$helpIcon <- renderUI({
    bsPopover(id = "popover_helpIcon", title = "AI Conclusions Prompt",
              content = "This section provides insights and conclusions drawn by the AI based on the analysis performed. It helps in understanding the significant findings and their implications.",
              placement = "right", trigger = "click")
  })
  
  output$helpIcon2 <- renderUI({
    bsPopover(id = "popover_helpIcon2", title = "AI Conclusions Prompt",
              content = "This section provides insights and conclusions drawn by the AI based on the analysis performed. It helps in understanding the significant findings and their implications.",
              placement = "right", trigger = "click")
  })
  
  # Activate popovers
  observe({
    shinyjs::runjs("
      $('#helpIcon').click(function() {
        $('#popover_helpIcon').popover('toggle');
      });
      
      $('#helpIcon2').click(function() {
        $('#popover_helpIcon2').popover('toggle');
      });
    ")
  })
  
  observeEvent(input$file_filter_kegg, {
    loadedData <- reactive({
      req(input$file_filter_kegg)
      df <- read.csv(input$file_filter_kegg$datapath, stringsAsFactors = FALSE)
      
      # Columnas obligatorias
      required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "metabolic_domain", "metabolic_subdomain")
      
      # Columnas disponibles para seleccionar (excluyendo las obligatorias)
      selectable_columns <- setdiff(names(df), required_columns)
      
      # UI para seleccionar las columnas de filtrado
      output$column_filter_selection <- renderUI({
        checkboxGroupInput("selected_filter_columns", "Select additional filter columns:", 
                           choices = selectable_columns, selected = NULL)
      })
      
      # Retornar el dataset original
      return(df)
    })
    
    # Renderizar la tabla filtrada con columnas obligatorias + seleccionadas
    output$table_filter_kegg <- renderDataTable({
      req(loadedData())  # Asegurar que los datos están cargados
      df <- loadedData()
      
      # Columnas obligatorias
      required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "metabolic_domain", "metabolic_subdomain")
      
      # Columnas seleccionadas por el usuario
      selected_columns <- input$selected_filter_columns
      
      # Si el usuario no selecciona ninguna, solo mostramos las obligatorias
      final_columns <- unique(c(required_columns, selected_columns))
      
      # Filtrar solo las columnas requeridas y seleccionadas
      df_filtered <- df[, intersect(names(df), final_columns), drop = FALSE]
      
      return(df_filtered)
    })
    
    # Descargar datos procesados con columnas obligatorias + seleccionadas
    output$downloadFile_filter_kegg <- downloadHandler(
      filename = function() {
        paste("filtered_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(loadedData())
        df <- loadedData()
        
        # Columnas obligatorias + seleccionadas
        required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "metabolic_domain", "metabolic_subdomain")
        selected_columns <- input$selected_filter_columns
        final_columns <- unique(c(required_columns, selected_columns))
        
        df_filtered <- df[, intersect(names(df), final_columns), drop = FALSE]
        write.csv(df_filtered, file, row.names = FALSE)
      }
    )
  })
  
  
  observeEvent(input$file_filter_gene_ontologies, {
    
    loadedData <- reactive({
      req(input$file_filter_gene_ontologies)
      df <- read.csv(input$file_filter_gene_ontologies$datapath, stringsAsFactors = FALSE)
      
      # Verificar si la columna ONT_NAME existe, si no, obtener los nombres desde QuickGO
      if (!"ONT_NAME" %in% names(df)) {
        unique_ontologies <- unique(df$ONTOLOGY)
        ontology_map <- setNames(future_sapply(unique_ontologies, get_ontology_name), unique_ontologies)
        df$ONT_NAME <- ontology_map[df$ONTOLOGY]  # Asigna nombres
      }
      
      # Columnas obligatorias
      required_columns <- c("ONTOLOGY", "ONT_NAME", "ONT_DESCRIPTION", "EA_VALUE", 
                            "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN")
      
      # Columnas disponibles para filtrado (excluyendo las obligatorias)
      selectable_columns <- setdiff(names(df), required_columns)
      
      # UI para seleccionar columnas adicionales
      output$column_selectors_go <- renderUI({
        checkboxGroupInput("selected_filter_columns_go", "Select additional filter columns:", 
                           choices = selectable_columns, selected = NULL)
      })
      
      return(df)
    })
    
    # Renderizar tabla filtrada con columnas obligatorias + seleccionadas
    output$table_filter_go <- renderDataTable({
      req(loadedData())  
      df <- loadedData()
      
      required_columns <- c("ONTOLOGY", "ONT_NAME", "ONT_DESCRIPTION", "EA_VALUE", 
                            "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN")
      
      selected_columns <- input$selected_filter_columns_go
      final_columns <- unique(c(required_columns, selected_columns))
      
      df_filtered <- df[, intersect(names(df), final_columns), drop = FALSE]
      return(df_filtered)
    })
    
    # Descargar datos procesados con columnas obligatorias + seleccionadas
    output$downloadFile_filter_go <- downloadHandler(
      filename = function() {
        paste("filtered_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(loadedData())
        df <- loadedData()
        
        required_columns <- c("ONTOLOGY", "ONT_NAME", "ONT_DESCRIPTION", "EA_VALUE", 
                              "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN")
        selected_columns <- input$selected_filter_columns_go
        final_columns <- unique(c(required_columns, selected_columns))
        
        df_filtered <- df[, intersect(names(df), final_columns), drop = FALSE]
        write.csv(df_filtered, file, row.names = FALSE)
      }
    )
  })
  
  # Lógica para procesar el archivo de KEGG
  observeEvent(input$file_kegg, {
    req(input$file_kegg)
    
    filename <- input$file_kegg$name
    data_kegg <- read.csv(input$file_kegg$datapath)
    
    # Llama a la función para procesar los datos de KEGG
    processed_kegg_data <- ancestors_kegg(data_kegg)
    
    # Generar la respuesta para descargar el archivo procesado
    output$downloadFile_kegg <- downloadHandler(
      filename = function() {
        "processed_kegg_data.csv"  # Especifica el nombre del archivo correctamente
      },
      content = function(file) {
        write.csv(processed_kegg_data, file, row.names = FALSE)
      }
    )
    
    # Mostrar la tabla procesada de KEGG
    output$table_kegg <- renderDataTable({  
      processed_kegg_data[, c("ONTOLOGY","ONT_DESCRIPTION", "GROUP_1","GROUP_2","UP_DOWN","metabolic_domain","metabolic_subdomain")]
    })
  })
  # Lógica para limpiar los datos de análisis funcional de KEGG
  observeEvent(input$clear_kegg, {
    # Eliminar la tabla y el archivo cargado
    output$table_kegg <- renderDataTable({ NULL })
    unlink(input$file_kegg$datapath)  # Eliminar el archivo cargado
  })
  
  # Lógica para procesar el archivo de Gene Ontologies
  observeEvent(input$file_gene_ontologies, {
    req(input$file_gene_ontologies)  # Asegurar que se haya cargado un archivo
    
    original_data <- read.csv(input$file_gene_ontologies$datapath, stringsAsFactors = FALSE)
    print(colnames(original_data))
    # Eliminar filas con NA en cualquier columna antes de procesar
    original_data <- na.omit(original_data)  # Alternativamente: original_data <- drop_na(original_data)
    print(colnames(original_data))
    # Llamar a la función ancestors para Gene Ontologies
    processed_data <- ancestors_gene_ontologies(original_data$ONTOLOGY, original_data$GROUP)
    colnames(processed_data) <- c("GROUP", "ONTOLOGY", "ANCESTORS")
    
    print(head(processed_data))
    
    
    colnames(processed_data) <- c("GROUP", "ONTOLOGY", "ANCESTORS")
    # Modificar FILE.GOEA.ANCESTORS
    processed_data$ANCESTORS <- lapply(processed_data$ANCESTORS, function(x) {
      x <- unlist(strsplit(x, ", "))
      x <- x[!x %in% c("GO:0003674", "GO:0005575", "GO:0008150")]
      x <- trimws(x)  
      paste(x, collapse = ", ")
    })
    
    # Realizar el merge por la columna "ONTOLOGY"
    original_data <- merge(original_data, processed_data[, c("ONTOLOGY", "ANCESTORS")], by = "ONTOLOGY", all.x = TRUE)
    print(head(original_data))
    # Summarise
    summary_data <- original_data %>%
      group_by(GROUP, ONTOLOGY, UP_DOWN) %>%
      summarise(
        ONT_DESCRIPTION = unique(ONT_DESCRIPTION),
        ONT_NAME = unique(ONT_NAME),
        EA_VALUE = mean(EA_VALUE),
        GROUP_1 = toString(unique(GROUP_1)),
        GROUP_2 = toString(unique(GROUP_2)),
        UP_DOWN = unique(UP_DOWN),
        ANCESTORS = toString(unique(ANCESTORS)),
        .groups = 'drop'
      ) %>%
      filter(ANCESTORS != "NA" & ANCESTORS != "")
    print(head(summary_data))
    # Asegurar que ANCESTORS en summary_data es character
    summary_data <- summary_data %>%
      mutate(ANCESTORS = as.character(ANCESTORS))
    
    # Asegurar que ANCESTORS en original_data es character antes del merge
    original_data <- original_data %>%
      mutate(ANCESTORS = as.character(ANCESTORS))
    
    # Unir nuevamente con original_data para recuperar columnas faltantes
    original_data <- left_join(summary_data, original_data, by = c("GROUP", "ONTOLOGY", "UP_DOWN", 
                                                                   "ONT_DESCRIPTION", "ONT_NAME", 
                                                                   "EA_VALUE", "GROUP_1", "GROUP_2", "ANCESTORS"))
    
    print(head(original_data))
    # Eliminar duplicados si es necesario
    original_data <- original_data %>%
      distinct()
    
    # Continuar con el procesamiento
    original_data <- original_data %>%
      mutate(first_ancestor = sapply(strsplit(ANCESTORS, ", "), find_first_matching_ancestor))
    
    original_data <- original_data %>%
      mutate(first_ancestor_name = sapply(first_ancestor, get_first_ancestor_name))
    
    # Eliminar la columna ANCESTORS si no la necesitas
    data_without_ancestors <- original_data[, !names(original_data) %in% "ANCESTORS"]
    
    # Definir la función downloadHandler para descargar el archivo original_data
    output$downloadFile_gene_ontologies <- downloadHandler(
      filename = function() {
        "go.csv"
      },
      content = function(file) {
        write.csv(data_without_ancestors, file, row.names = FALSE)
      }
    )
    
    # Mostrar la tabla procesada en la interfaz de usuario
    output$table_gene_ontologies <- renderDataTable({  # Usar renderDataTable en lugar de renderTable
      data_without_ancestors
    })
  })
  
  # Lógica para limpiar los datos de análisis funcional de KEGG
  observeEvent(input$clear_go, {
    # Eliminar la tabla y el archivo cargado
    output$table_gene_ontologies <- renderDataTable({ NULL })
    unlink(input$file_gene_ontologies$datapath)  # Eliminar el archivo cargado
  })
  
  
  
  observeEvent(input$file_kegg_network, {
    req(input$file_kegg_network)
    dataset <- read.csv(input$file_kegg_network$datapath, stringsAsFactors = FALSE)
    dataset <- analyze_regulation(dataset) 
    
    # Define required columns
    required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "metabolic_domain", "metabolic_subdomain")
    
    # Get filterable columns
    filter_columns <- setdiff(names(dataset), required_columns)
    
    # Actualiza el selectInput de group_by_kegg con las opciones de agrupación
    updateSelectInput(session, "group_by_kegg", choices = c("Domain", "Subdomain", "Relation"))
    
    output$dynamic_filters <- renderUI({
      req(dataset)  
      filter_columns <- setdiff(names(dataset), required_columns)
      
      if (length(filter_columns) == 0) {
        return(tagList("No valid columns for filtering"))
      }
      
      tagList(
        lapply(filter_columns, function(col) {
          selectInput(col, col, choices = unique(dataset[[col]]), multiple = FALSE)
        })
      )
    })
    
    # Reactive expression para filtrar dinámicamente los datos
    filtered_data_reactive <- reactive({
      req(dataset)
      filtered_data <- dataset
      
      for (col in filter_columns) {
        if (!is.null(input[[col]]) && input[[col]] != "") {
          filtered_data <- filtered_data %>% filter(.data[[col]] == input[[col]])
        }
      }
      
      return(filtered_data)
    })
    
    # Actualiza las opciones de group_kegg basadas en los datos filtrados
    observe({
      filtered_data <- filtered_data_reactive()
      available_groups <- unique(filtered_data$GROUP)
      updateSelectInput(session, "group_kegg", choices = c("All groups", available_groups))
    })
    
    # Observa cambios en filtros, group_kegg y group_by_kegg para actualizar la red
    observe({
      req(input$group_kegg, input$group_by_kegg)
      
      filtered_data <- filtered_data_reactive()
      
      if (input$group_kegg != "All groups") {
        filtered_data <- filtered_data %>% filter(GROUP == input$group_kegg)
      }
      
      # Obtener la columna de agrupación correctamente
      group_by_col <- switch(input$group_by_kegg,
                             "Domain" = "metabolic_domain",
                             "Subdomain" = "metabolic_subdomain",
                             "Relation" = "ONT_DESCRIPTION")
      
      if (is.null(group_by_col) || !(group_by_col %in% names(filtered_data))) {
        return(NULL)
      }
      
      # Crear listas de nodos y aristas
      all_nodes <- list()
      all_edges <- list()
      
      group_colors <- rainbow(length(unique(filtered_data[[group_by_col]])))
      names(group_colors) <- unique(filtered_data[[group_by_col]])
      
      unique_ontologies <- unique(filtered_data$ONTOLOGY)
      
      for (ontology in unique_ontologies) {
        ontology_data <- filtered_data %>% filter(ONTOLOGY == ontology)
        all_groups <- unique(c(ontology_data$GROUP_1, ontology_data$GROUP_2))
        
        nodes <- data.frame(
          id = paste(ontology, all_groups, sep = "_"), 
          label = ontology, 
          group = all_groups,
          Domain = unique(ontology_data$metabolic_domain),
          Subdomain = unique(ontology_data$metabolic_subdomain),
          Relation = unique(ontology_data$ONT_DESCRIPTION),
          color.background = group_colors[all_groups], 
          title = unique(ontology_data$ONT_DESCRIPTION)
        )
        
        arrow_types <- ifelse(ontology_data$UP_DOWN == "UP", "to", 
                              ifelse(ontology_data$UP_DOWN == "DOWN", "to", "none"))
        
        edges <- data.frame(
          from = paste(ontology, ontology_data$GROUP_1, sep = "_"),
          to = paste(ontology, ontology_data$GROUP_2, sep = "_"),
          arrows = arrow_types,
          color = ifelse(ontology_data$UP_DOWN == "UP", "blue", 
                         ifelse(ontology_data$UP_DOWN == "DOWN", "red", "black")),
          title = paste("EA_value:", ontology_data$EA_VALUE)
        )
        
        all_nodes[[ontology]] <- nodes
        all_edges[[ontology]] <- edges
      }
      
      all_nodes_combined <- do.call(rbind, all_nodes)
      all_edges_combined <- do.call(rbind, all_edges)
      
      # Renderizar la red
      output$network_kegg <- renderVisNetwork({
        visNetwork(all_nodes_combined, edges = all_edges_combined, main = "KEGG Network", width = "100%", height = "100%") %>%
          visNodes(color = list(border = "black"), shadow = TRUE) %>%
          visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 2)), width = 3) %>%
          visOptions(highlightNearest = TRUE, selectedBy = input$group_by_kegg) %>%
          visLegend(main = "Groups", useGroups = TRUE, zoom = FALSE)
      })
    })
    
    output$download_network_kegg <- downloadHandler(
      filename = function() { paste('network-', Sys.Date(), '.html', sep = '') },
      content = function(con) {
        visNetwork(all_nodes_combined, edges = all_edges_combined, main = "KEGG Network", width = "100%", height = "100%") %>%
          visNodes(color = list(border = "black"), shadow = TRUE) %>%
          visEdges(arrows = "to") %>%
          visOptions(highlightNearest = TRUE, selectedBy = input$group_by_kegg) %>%
          visLegend(main = "Groups", useGroups = TRUE, zoom = FALSE) %>%
          visSave(con)
      }
    )
    
    arrow_colors <- c("Up-regulated" = "blue", "Down-regulated" = "red", "Neutral" = "black")
    
    output$legend <- renderUI({
      tagList(
        tags$h3("Legend"),
        tags$h4("Arrows"),
        lapply(names(arrow_colors), function(arrow) {
          tags$div(style = "margin-bottom: 10px;",
                   tags$div(style = paste("display: inline-block; width: 20px; height: 2px; background-color:", arrow_colors[arrow], ";")), 
                   arrow)
        })
      )
    })
    
    output$text_kegg <- renderText({
      "Download a comprehensive prompt designed to assist in the analysis of Gene Ontology (GO) enrichment across your study. 
      This prompt includes detailed information on the most enriched ontologies, samples analyzed, descriptions of GO terms, and more."
    })
    
    output$download_prompt_kegg <- downloadHandler(
      filename = function() { paste("prompt-", Sys.Date(), ".txt", sep = "") },
      content = function(file) {
        filtered_data <- filtered_data_reactive()
        
        if (nrow(filtered_data) == 0) {
          writeLines("No data available for the selected filters.", file)
          return()
        }
        
        top_ontologies <- filtered_data %>%
          arrange(desc(EA_VALUE)) %>%
          head(20)
        
        intro_text <- "Write an exhaustive analysis focusing on biological experimental conclusions to learn how the diseases are doing on a Gene Ontology (GO) enrichment dataset for CU or EC, taking the main ideas for all the dataset without specifying in each of them. The dataset includes enrichment information for the 15 most enriched ontologies in the dataset. Each entry in the dataset has the following attributes:

                          Ontology: The KEGG identifier.
                          Sample: Sample used for enrichment analysis.
                          Description: Description of the KEGG term.
                          Disease: Disease studied.
                          Group_1: First group for comparison.
                          Group_2: Second group for comparison.
                          ea_value: Enrichment value.
                          metabolic_domain: The first level hierarchy pathway in the KEGG hierarchy.
                          metabolic_subdomain: Second level of hierarchy pathway in the KEGG hierarchy.
                          
                          So here are the gene ontologies to analyze:\\n\\n"
        
        ontology_list <- apply(top_ontologies, 1, function(row) {
          paste("- Ontology:", row["ONTOLOGY"],
                "| Sample:", row["Sample"], 
                "| Description:", row["ONT_DESCRIPTION"],
                "| Disease:", row["Disease"],
                "| Group_1:", row["GROUP_1"],
                "| Group_2:", row["GROUP_2"],
                "| ea_value:", row["EA_VALUE"],
                "| metabolic_domain:", row["metabolic_domain"],
                "| metabolic_subdomain:", row["metabolic_subdomain"])
        })
        
        final_text <- paste(intro_text, paste(ontology_list, collapse = "\n"), sep = "\n")
        
        writeLines(final_text, file)
      }
    )
  })
  
  
  # Lógica para limpiar los datos de análisis funcional de KEGG
  observeEvent(input$clear_kegg_network, {
    # Eliminar la tabla y el archivo cargado
    output$network_kegg <- renderVisNetwork({ NULL })
    unlink(input$file_kegg_network$datapath)  # Eliminar el archivo cargado
  })
  
  
  
  observeEvent(input$file_heatmap, {
    req(input$file_heatmap)
    dataset <- read.csv(input$file_heatmap$datapath, stringsAsFactors = FALSE)
    dataset <- analyze_regulation(dataset)
    
    required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "metabolic_domain", "metabolic_subdomain")
    filter_columns <- setdiff(names(dataset), required_columns)
    
    updateSelectInput(session, "group_by_heatmap", choices = c("Domain", "Subdomain", "Relation"))
    
    output$dynamic_filters_heatmap <- renderUI({
      req(dataset)
      if (length(filter_columns) == 0) {
        return(tagList("No valid columns for filtering"))
      }
      tagList(
        lapply(filter_columns, function(col) {
          selectInput(col, col, choices = unique(dataset[[col]]), multiple = FALSE)
        })
      )
    })
    
    # Reactive expression para obtener datos filtrados dinámicamente
    filtered_data_reactive <- reactive({
      req(dataset)
      filtered_data <- dataset
      for (col in filter_columns) {
        if (!is.null(input[[col]]) && input[[col]] != "") {
          filtered_data <- filtered_data %>% filter(.data[[col]] == input[[col]])
        }
      }
      return(filtered_data)
    })
    
    observe({
      filtered_data <- filtered_data_reactive()
      updateSelectInput(session, "group_heatmap", choices = c("All groups", unique(filtered_data$GROUP)))
    })
    
    observe({
      req(input$group_heatmap, input$group_by_heatmap)
      
      filtered_data <- filtered_data_reactive()
      
      if (input$group_heatmap != "All groups") {
        filtered_data <- filtered_data %>% filter(GROUP == input$group_heatmap)
      }
      
      # Obtener la columna de agrupación correctamente
      group_by_col <- switch(input$group_by_heatmap,
                             "Domain" = "metabolic_domain",
                             "Subdomain" = "metabolic_subdomain",
                             "Relation" = "ONT_DESCRIPTION")
      
      # Verificar si la columna de agrupación es válida
      if (is.null(group_by_col) || !(group_by_col %in% names(filtered_data))) {
        return(NULL)
      }
      
      # Aplicar el filtrado basado en la agrupación
      filtered_data <- filtered_data %>%
        filter(!is.na(.data[[group_by_col]])) %>%
        group_by(GROUP, .data[[group_by_col]]) %>%
        summarise(EA_VALUE = mean(EA_VALUE, na.rm = TRUE), .groups = 'drop')
      
      output$heatmap_plot <- renderPlotly({
        req(nrow(filtered_data) > 0)  # Asegura que haya datos antes de graficar
        
        p <- ggplot(filtered_data, aes(x = .data[[group_by_col]], y = GROUP, fill = EA_VALUE)) +
          geom_tile(color = "white") +
          scale_fill_gradient(low = "blue", high = "red") +
          labs(
            title = paste("Heatmap - Agrupado por:", input$group_by_heatmap),
            x = input$group_by_heatmap, 
            y = "Grupo", 
            fill = "EA_VALUE"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggplotly(p)
      })
    })
    
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste('heatmap-', Sys.Date(), '.png', sep = '')
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png", width = 8, height = 6)
      }
    )
  })
  
  
  
  
  observeEvent(input$file_volcano, {
    req(input$file_volcano)
    dataset <- read.csv(input$file_volcano$datapath, stringsAsFactors = FALSE)
    dataset <- analyze_regulation(dataset)
    
    # Actualización de selectInput para los valores únicos
    updateSelectInput(session, "sample_volcano", choices = unique(dataset$sample))
    updateSelectInput(session, "disease_volcano", choices = unique(dataset$Disease))
    updateSelectInput(session, "group_by_volcano", choices = c("Domain", "Subdomain", "Relation"))
    
    # Observador para actualizar las ontologías en base al nivel de agrupamiento
    observeEvent(input$group_by_volcano, {
      filtered_data <- dataset %>%
        filter(sample == input$sample_volcano, Disease == input$disease_volcano)
      
      if (input$group_by_volcano == "Domain") {
        unique_ontologies <- unique(filtered_data$metabolic_domain)
        updateSelectInput(session, "ontology_volcano", choices = unique_ontologies)
      } else if (input$group_by_volcano == "Subdomain") {
        unique_ontologies <- unique(filtered_data$metabolic_subdomain)
        updateSelectInput(session, "ontology_volcano", choices = unique_ontologies)
      } else {
        unique_ontologies <- unique(filtered_data$ONT_DESCRIPTION)
        updateSelectInput(session, "ontology_volcano", choices = unique_ontologies)
      }
    })
    
    observeEvent(input$ontology_volcano, {
      req(input$ontology_volcano)
      
      # Filtrar datos basados en las selecciones del usuario
      filtered_data <- dataset %>%
        filter(sample == input$sample_volcano, Disease == input$disease_volcano)
      
      group_by_col <- switch(input$group_by_volcano,
                             "Domain" = "metabolic_domain",
                             "Subdomain" = "metabolic_subdomain",
                             "Relation" = "ONT_DESCRIPTION")
      
      # Filtrar para obtener las ontologías bajo el grupo seleccionado
      volcano_data <- filtered_data %>%
        filter(!!sym(group_by_col) == input$ontology_volcano) %>%
        select(log2FoldChange, pvalue, gene)
      
      # Transformar el p-value a -log10(p-value) para el volcano plot
      volcano_data$negLogPval <- -log10(volcano_data$pvalue)
      
      # Crear el volcano plot
      output$volcano_plot <- renderPlotly({
        # Crear el ggplot para el volcano plot
        p <- ggplot(volcano_data, aes(x = log2FoldChange, y = negLogPval, text = gene)) +
          geom_point(aes(color = negLogPval > -log10(0.05) & abs(log2FoldChange) > 1), 
                     size = 2, alpha = 0.6) +
          scale_color_manual(values = c("black", "red")) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
          labs(x = "log2(Fold Change)", y = "-log10(p-value)", 
               title = paste("Volcano Plot: Sample:", input$sample_volcano, 
                             "Disease:", input$disease_volcano)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        # Hacer el plot interactivo con plotly
        ggplotly(p, tooltip = "text")
      })
    })
    
    # Opción para descargar el volcano plot
    output$download_volcano <- downloadHandler(
      filename = function() {
        paste('volcano-', Sys.Date(), '.png', sep = '')
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png", width = 8, height = 6)
      }
    )
  })
  
  observeEvent(input$file_go_network, {
    req(input$file_go_network)
    dataset <- read.csv(input$file_go_network$datapath, stringsAsFactors = FALSE)
    dataset <- analyze_regulation(dataset) 
    
    # Definir columnas requeridas
    required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "GROUP_1", "GROUP_2", "UP_DOWN", "first_ancestor", "first_ancestor_name", "ONT_NAME")
    
    # Obtener las columnas filtrables (excluyendo las obligatorias)
    filter_columns <- setdiff(names(dataset), required_columns)
    
    # Actualiza el selectInput de group_go con las opciones de agrupación
    updateSelectInput(session, "group_go", choices =  unique(dataset$GROUP))
    updateSelectInput(session, "ont_description", choices = unique(dataset$ONT_DESCRIPTION))
    
    # Renderiza los filtros dinámicos
    output$dynamic_filters_go <- renderUI({
      req(dataset)  
      filter_columns <- setdiff(names(dataset), required_columns)
      
      tagList(
        lapply(filter_columns, function(col) {
          selectInput(col, col, choices = unique(dataset[[col]]), multiple = FALSE)
        })
      )
    })
    
    # Expresión reactiva para filtrar los datos dinámicamente
    filtered_data_reactive <- reactive({
      req(dataset)
      filtered_data <- dataset
      
      for (col in filter_columns) {
        if (!is.null(input[[col]]) && input[[col]] != "") {
          filtered_data <- filtered_data %>% filter(.data[[col]] == input[[col]])
        }
      }
      
      return(filtered_data)
    })
    
    # Actualiza las opciones de group_go basadas en los datos filtrados
    observe({
      filtered_data <- filtered_data_reactive()
      available_groups <- unique(filtered_data$GROUP)
      updateSelectInput(session, "group_go", choices = c("All groups", available_groups))
    })
    
    # Observa los cambios en los filtros, group_go y otros filtros dinámicos para actualizar la red
    observe({
      req(input$group_go, input$ont_description)
      
      filtered_data <- filtered_data_reactive()
      
      # Filtrar por el ont_description
      filtered_data <- filtered_data %>% filter(ONT_DESCRIPTION == input$ont_description)
      
      if (input$group_go!= "All groups") {
        filtered_data <- filtered_data %>% filter(GROUP == input$group_go)
      }
      
      # Crear listas de nodos y aristas
      all_nodes <- list()
      all_edges <- list()
      
      group_colors <- rainbow(length(unique(filtered_data$GROUP)))
      names(group_colors) <- unique(filtered_data$GROUP)
      
      unique_ontologies <- unique(filtered_data$ONTOLOGY)
      
      for (ontology in unique_ontologies) {
        ontology_data <- filtered_data %>% filter(ONTOLOGY == ontology)
        all_groups <- unique(c(ontology_data$GROUP_1, ontology_data$GROUP_2))
        
        # Crear nodos para todos los grupos únicos
        nodes <- data.frame(
          id = paste(ontology, all_groups, sep = "_"), 
          label = ontology,
          group = all_groups,
          first_anc = unique(ontology_data$first_ancestor_name),
          color.background = group_colors[all_groups], 
          title = unique(ontology_data$ONT_NAME)
        )
        
        # Crear tipos de flechas en función de los valores de UP_DOWN
        arrow_types <- ifelse(ontology_data$UP_DOWN == "UP", "to", 
                              ifelse(ontology_data$UP_DOWN == "DOWN", "to", "none"))
        
        # Crear aristas basadas en los datos filtrados
        edges <- data.frame(
          from = paste(ontology, ontology_data$GROUP_1, sep = "_"),
          to = paste(ontology, ontology_data$GROUP_2, sep = "_"),
          arrows = arrow_types,
          color = ifelse(ontology_data$UP_DOWN == "UP", "blue", 
                         ifelse(ontology_data$UP_DOWN == "DOWN", "red", "black")),
          title = paste("EA_value:", ontology_data$EA_VALUE)
        )
        
        all_nodes[[ontology]] <- nodes
        all_edges[[ontology]] <- edges
      }
      
      # Unir todos los nodos y aristas
      all_nodes_combined <- do.call(rbind, all_nodes)
      all_edges_combined <- do.call(rbind, all_edges)
      
      # Verifica que los datos de nodos y aristas estén correctamente formateados
      if (nrow(all_nodes_combined) == 0 || nrow(all_edges_combined) == 0) {
        return(NULL)  # No hacer nada si no hay nodos o aristas
      }
      
      # Renderizar la red
      output$network_go <- renderVisNetwork({
        visNetwork(all_nodes_combined, edges = all_edges_combined, main = "GO Network", width = "100%", height = "100%") %>%
          visNodes(color = list(border = "black"), shadow = TRUE) %>%
          visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 2)), width = 3) %>%
          visOptions(highlightNearest = TRUE, selectedBy = "first_anc") %>%
          visLegend(main = "Groups", useGroups = TRUE, zoom = FALSE)
      })
    })
    
    arrow_colors <- c("Up-regulated" = "blue", "Down-regulated" = "red", "Neutral" = "black")
    
    
    # Renderizado dinámico de la leyenda
    output$legend_go <- renderUI({
      tagList(
        tags$h3("Legend"),
        tags$h4("Arrows"),
        lapply(names(arrow_colors), function(arrow) {
          tags$div(style = "margin-bottom: 10px;",
                   tags$div(style = paste("display: inline-block; width: 20px; height: 2px; background-color:", arrow_colors[arrow], ";")), 
                   arrow)
        })
      )
    })
    
    output$text_go <-renderText({
      "Download a comprehensive prompt designed to assist in the analysis of Gene Ontology (GO) enrichment across your study. 
        This prompt includes detailed information on the most enriched ontologies, samples analyzed, descriptions of GO terms, and more. It is a valuable 
        tool for scientists and researchers conducting in-depth biological data analysis."
    })
    
    
    
    output$download_prompt_go <- downloadHandler(
      filename = function() {
        paste("prompt-", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        top_ontologies <- filtered_data %>%
          arrange(desc(EA_VALUE)) %>%
          head(20)
        
        intro_text <- "Write an exhaustive analysis focusing on biological experimental conclusions to learn how the diseases are doing on a Gene Ontology (GO) enrichment dataset for CU or EC, taking the main ideas for all the dataset without specifying in each of them. The dataset includes enrichment information for the 15 most enriched ontologies in the dataset. Each entry in the dataset has the following attributes:

Ontology: The GO identifier.
Sample: Sample used for enrichment analysis.
Description: Description of the GO term.
Disease: Disease studied.
Group_1: First group for comparison.
Group_2: Second group for comparison.
ea_value: Enrichment value.
first_ancestor: The most general ancestor node in the GO hierarchy.

So here are the gene ontologies to analyze:\\n\\n"
        
        ontology_list <- character(nrow(top_ontologies))
        
        for (i in seq_len(nrow(top_ontologies))) {
          ontology_info <- paste(
            "- Ontology:", top_ontologies$ONTOLOGY[i],
            "| Sample:", top_ontologies$sample[i],
            "| Description:", top_ontologies$ONT_DESCRIPTION[i],
            "| Disease:", top_ontologies$Disease[i],
            "| Group_1:", top_ontologies$GROUP_1[i],
            "| Group_2:", top_ontologies$GROUP_2[i],
            "| ea_value:", top_ontologies$EA_VALUE[i],
            "| first_ancestor:", top_ontologies$first_ancestor_name[i],
            sep = " "
          )
          
          # Append each ontology_info to ontology_list
          ontology_list[i] <- ontology_info
        }
        
        # Combine all ontology_list elements into a single string separated by newline characters
        ontology_list_text <- paste(ontology_list, collapse = "\\n")
        
        # Combine the introductory text and the ontology list text
        final_text <- paste(intro_text, ontology_list_text, sep = "")
        
        # Write the final text to the file
        writeLines(final_text, file)
      })
  })
  
  output$download_network_go <- downloadHandler(
    filename = function() {
      paste('network-', Sys.Date(), '.html', sep='')
    },
    content = function(con) 
      
      visNetwork(nodes, edges = filtered_edges, main = "Sample", width = "100%") %>%
      visNodes(color = list(background = "white", border = "black"), size = "value") %>%
      visGroups(groupname = "group", color = list(border = "black", background = rainbow(length(group_names), start = 0, end = 1)), legend = TRUE) %>%
      visPhysics(
        enabled = TRUE,
        repulsion = list(nodeDistance = 1)
      ) %>%
      visOptions(highlightNearest = TRUE, selectedBy = "first_anc") %>%
      visLegend(main = "Group", useGroups = TRUE, zoom = FALSE) %>% 
      visSave(con)
  )
  
  
  # Lógica para limpiar los datos de análisis funcional de KEGG
  observeEvent(input$clear_go_network, {
    # Eliminar la tabla y el archivo cargado
    output$network_go <- renderDataTable({ NULL })
    unlink(input$file_go_network$datapath)  # Eliminar el archivo cargado
    
    
  })
  
  # Lógica para procesar el archivo y generar los gráficos de barras con estandarización Z-score
  observeEvent(input$file_heatmap_go, {
    req(input$file_heatmap_go)
    dataset <- read.csv(input$file_heatmap_go$datapath, stringsAsFactors = FALSE)
    dataset <- analyze_regulation(dataset)
    
    # Definir columnas requeridas (obligatorias)
    required_columns <- c("ONTOLOGY", "ONT_DESCRIPTION", "EA_VALUE", "GROUP", "UP_DOWN", "first_ancestor_name", "first_ancestor", "ONT_NAME","GROUP_1","GROUP_2")
    
    # Obtener las columnas filtrables (excluyendo las obligatorias)
    filter_columns <- setdiff(names(dataset), required_columns)
    
    # Actualiza el selectInput de "ont_description" y "group"
    updateSelectInput(session, "ont_description_heatmap_go", choices = unique(dataset$ONT_DESCRIPTION))
    updateSelectInput(session, "group_heatmap_go", choices = unique(dataset$GROUP))
    
    # Renderiza los filtros dinámicos
    output$dynamic_filters_heatmap_go <- renderUI({
      req(dataset)  
      filter_columns <- setdiff(names(dataset), required_columns)
      
      tagList(
        lapply(filter_columns, function(col) {
          selectInput(col, col, choices = unique(dataset[[col]]), multiple = FALSE)
        })
      )
    })
    
    # Expresión reactiva para filtrar los datos dinámicamente
    filtered_data_reactive <- reactive({
      req(dataset)
      filtered_data <- dataset
      
      # Filtrar por los filtros dinámicos
      for (col in filter_columns) {
        if (!is.null(input[[col]]) && input[[col]] != "") {
          filtered_data <- filtered_data %>% filter(.data[[col]] == input[[col]])
        }
      }
      
      return(filtered_data)
    })
    
    # Actualiza las opciones de "group" basadas en los datos filtrados
    observe({
      filtered_data <- filtered_data_reactive()
      available_groups <- unique(filtered_data$GROUP)
      updateSelectInput(session, "group_heatmap_go", choices = c("All groups", available_groups))
    })
    
    # Gráficos de barras: UP y DOWN regulación
    observe({
      req(input$group_heatmap_go, input$ont_description_heatmap_go)
      
      # Filtrar los datos según los filtros aplicados
      filtered_data <- filtered_data_reactive()
      
      # Filtrar por el ont_description
      filtered_data <- filtered_data %>% filter(ONT_DESCRIPTION == input$ont_description_heatmap_go)
      
      # Filtrar por group si no es "All groups"
      if (input$group_heatmap_go != "All groups") {
        filtered_data <- filtered_data %>% filter(GROUP == input$group_heatmap_go)
      }
      
      # Asegurarse de que los valores de EA sean > 0
      filtered_data <- filtered_data %>% filter(EA_VALUE > 0)
      
      # Estandarización Z-score
      z_score <- function(x) {
        (x - mean(x)) / sd(x)
      }
      
      up_data <- filtered_data %>%
        filter(UP_DOWN == "UP") %>%
        group_by(ONTOLOGY, first_ancestor_name) %>%  # Agrupar también por first_ancestor
        summarize(EA_VALUE = mean(EA_VALUE, na.rm = TRUE), .groups = 'drop') %>%
        mutate(EA_VALUE = z_score(EA_VALUE))
      
      down_data <- filtered_data %>%
        filter(UP_DOWN == "DOWN") %>%
        group_by(ONTOLOGY, first_ancestor_name) %>%  # Agrupar también por first_ancestor
        summarize(EA_VALUE = mean(EA_VALUE, na.rm = TRUE), .groups = 'drop') %>%
        mutate(EA_VALUE = z_score(EA_VALUE))
      
      # Crear el gráfico de barras para "up"
      output$barplot_up_go <- renderPlotly({
        p_up <- ggplot(up_data, aes(x = ONTOLOGY, y = EA_VALUE, fill = first_ancestor_name)) +  # Cambiar fill a first_ancestor
          geom_bar(stat = "identity", position = "dodge") +
          labs(
            title = paste("Regulación Up - Muestra:", input$sample_heatmap_go, 
                          "y Enfermedad:", input$disease_heatmap_go),
            x = "Ontologías",
            y = "Valor de EA"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
          )
        
        ggplotly(p_up, tooltip = c("x", "y", "fill"))
      })
      
      # Crear el gráfico de barras para "down"
      output$barplot_down_go <- renderPlotly({
        p_down <- ggplot(down_data, aes(x = ONTOLOGY, y = EA_VALUE, fill = first_ancestor_name)) +  # Cambiar fill a first_ancestor
          geom_bar(stat = "identity", position = "dodge") +
          labs(
            title = paste("Regulación Down - Muestra:", input$sample_heatmap_go, 
                          "y Enfermedad:", input$disease_heatmap_go),
            x = "Ontologías",
            y = "Valor de EA"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
          )
        
        ggplotly(p_down, tooltip = c("x", "y", "fill"))
      })
    })
    
    # Descargar los gráficos de barras
    output$download_barplot_up_go <- downloadHandler(
      filename = function() {
        epaste('barplot_up_go-', Sys.Date(), '.png', sep = '')
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png", width = 8, height = 6)
      }
    )
    
    output$download_barplot_down_go <- downloadHandler(
      filename = function() {
        paste('barplot_down_go-', Sys.Date(), '.png', sep = '')
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png", width = 8, height = 6)
      }
    )
  })
}

shinyApp(ui,server)