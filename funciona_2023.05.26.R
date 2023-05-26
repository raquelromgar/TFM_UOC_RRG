###############################################################################
###############################################################################
##                                                                           ##
## Esta es una aplicaci?n web Shiny. Puede ejecutar la aplicaci?n haciendo   ##
##           clic en el bot?n 'Ejecutar aplicaci?n' de arriba.               ##
##                                                                           ##  
##                                                                           ##
###############################################################################
###############################################################################

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(dplyr)
library(readr)
library(knitr)
library(rmarkdown)


############################################################
##    Definir la interfaz de usuario para la aplicación   ##
############################################################

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "cosmo"),
  
  titlePanel("OncoReporter"),
  
  tabsetPanel(id = "tablas",
    tabPanel(
      title = "Archivos y datos",
      fluidRow(
        h3("Sube tus archivos a la categoria correspondiente")
      ),
      fluidRow(
        column(5,
               fileInput("SNV",label = "Fichero para SNVs, MNVs e INDELs", multiple = F, accept = c(".tsv")),
               fileInput("CNV",label = "Fichero para CNVs", multiple = F, accept = c(".tsv")),
               br(),
        ),
        
        column(6,
               h4("Datos"),
               selectInput("sexo", label = "Selecciona el sexo",
                           choices = c("Masculino" = "M", "Femenino" = "F")),
               textInput("edad", label = "Edad"),
               textInput("tumor", label = "Tipo de tumor"),
               textInput("localizacion", label = "Localización"),
               textInput("céls_tum", label = "Porcentaje de células tumorales"),
        ),
      ),
      fluidRow(
        #Ejecutar priorización
        actionButton("submit", label= "Priorización de variantes")
      )
    ),
    
    tabPanel(
      title = "Resultados SNVs, MNVs e INDELs",
      br(),
      #Tablas Resultados
      DT::dataTableOutput("filtered_table_SNV"),
      br(),
      #Descarga resultados
      downloadButton("desc_res_SNV", label = "Guardar resultado (SNVs, MNVs e INDELs)"),
    ), 
    tabPanel(
      title = "Resultados CNVs",
      br(),
      #Tablas Resultados
      DT::dataTableOutput("filtered_table_CNV"),
      br(),
      #Descarga resultados
      downloadButton("desc_res_CNV", label = "Guardar resultado (CNVs)"),
    ), 
    tabPanel(
      title = "Clasificación de las variantes",
      br(),
      #Tablas Clasificación
      DT::dataTableOutput("classification"),
      actionButton("generar_informe", "Generar informe"),
      downloadButton("report", "Descargar Informe"),
      br(),
    ),
  )
)


############################################################
##    Definir la lógica del servidor para la aplicación   ##
############################################################     

# Definir la lógica del servidor
server <- function(input, output, session) {
  
    # Función para leer los archivos TSV
  read_snv_file <- function() {
    req(input$SNV) # Verificar si se ha subido un archivo
    
    # Leer el archivo TSV y retornar el dataframe
    snv_data <- read_tsv(input$SNV$datapath, comment = "#")
    return(snv_data)   
  }
  
  read_cnv_file <- function() {
    req(input$CNV)   # Verificar si se ha subido un archivo
    
    # Leer el archivo TSV y retornar el dataframe
    cnv_data <- read_tsv(input$CNV$datapath, comment = "#")
    return(cnv_data)
  }
  
  # Función para aplicar filtros al dataframe de SNV
  filter_snv_data <- function(data) {
    # Aplicar los filtros que desees en esta función
    filtered_data_snv <- data %>%
      filter(Filter == "PASS") %>%
      filter(!(Type %in% c("CNV", "EXPR_CONTROL", "REF", "ASSAYS_5P_3P"))) %>%
      filter(`Variant Effect` != "synonymous" | is.na(`Variant Effect`)) %>%
      filter(`UCSC Common SNPs`!="YES" | is.na(`UCSC Common SNPs`)) %>%
      filter(!grepl("benign", ClinVar, ignore.case=TRUE) | is.na(ClinVar)) %>%
      filter(`P-Value`<=0.05) %>%
      filter(grepl("exonic", Location, ignore.case=TRUE)) %>%
      filter(`Allele Frequency %`>=3) %>%
      filter(Coverage>=300)
    return(filtered_data_snv)
  }
  
  # Función para aplicar filtros al dataframe de CNV
  filter_cnv_data <- function(data) {
    # Aplicar los filtros que desees en esta función
    filtered_data_cnv <- data %>%
      filter(`CNV P-Value`<=0.05) %>%
      filter(`Copy Number`>4 | `Copy Number`<=0.01)
    return(filtered_data_cnv)
  }
  
  
  # Observador para leer el archivo y aplicar los filtros cuando se pulsa el botón "submit"
  observeEvent(input$submit, {
    snv_data <- read_snv_file() # Leer el archivo TSV
    filtered_data_snv <- filter_snv_data(snv_data) # Aplicar filtros de SNV
    cnv_data <- read_cnv_file() # Leer el archivo TSV
    filtered_data_cnv <- filter_cnv_data(cnv_data) # Aplicar filtros de CNV
    
    # Cambiar a la página de resultados cuando se puse el botón submit
    updateTabsetPanel(session, "tablas", "Resultados SNVs, MNVs e INDELs")
    
    # Mostrar las tablas con los valores filtrados
    output$filtered_table_SNV <- DT::renderDataTable({
      datatable(
        filtered_data_snv,
        filter = "top",
        extensions = 'FixedHeader',
        caption = htmltools::tags$caption(
          style = 'caption-side: bottom; text-align: center;',
          'Tabla 1: ', htmltools::em('Variantes de tipo SNV, MNV o INDEL')
        ),
        # Opciones adicionales de la tabla
        options = list(
          dom = 'lBfrtip',
          fixedHeader = TRUE,
          columnDefs = list(list(
          targets = c(42,58,59),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 60 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 60) + '...</span>' : data;",
            "}")
        )),
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#5dadbd', 'color': '#fff'});",
          "}")
        ), 
          callback = JS('table.page().draw(false);')
   ) 
})

    output$filtered_table_CNV <- DT::renderDataTable({
      datatable(
        filtered_data_cnv,
        filter = "top",
        extensions = 'FixedHeader',
        caption = htmltools::tags$caption(
          style = 'caption-side: bottom; text-align: center;',
          'Tabla 2: ', htmltools::em('Variantes de tipo CNV')
        ),
        # Opciones adicionales de la tabla
        options = list(
          dom = 'lBfrtip',
          fixedHeader = TRUE,
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#5dadbd', 'color': '#fff'});",
            "}")
        ), 
        callback = JS('table.page().draw(false);')
      ) 
    })
    
    #Función de clasificación
    Varsome<-data.frame(Gen = character(),
                        cDNA = character(),
                        dbSNP = character(),
                        url = character(),
                        stringsAsFactors = FALSE)
    
    
    classif<-function(Varsome, df){
      for (i in 1:nrow(df)){
        Gen<-df$Genes[i]
        cDNA<- paste0(df$Transcript[i], ":", df$Coding[i])
        dbSNP<- df$dbSNP[i]
        url <- paste0("https://varsome.com/variant/hg38/", df$dbSNP[i])
        
        nuevaFila <- data.frame(Gen = Gen, cDNA = cDNA, dbSNP = dbSNP, url = url)
        
        Varsome<-rbind(Varsome, nuevaFila)
      }
      
      return(Varsome)
    }
    
    varsome_vector <- classif(Varsome, filtered_data_snv)
    
    
    # Generar la tabla de clasificación
    output$classification <- DT::renderDataTable({
      datatable(varsome_vector, escape = FALSE)
    })
    
  })
    
 
  
    #BOTONES DESCARGA 
  observeEvent(input$generar_informe, {
    output$report <- downloadHandler(
      filename = "Informe de variantes.pdf",
      content =  function(file){
        # Copia el archivo del informe en un directorio temporal 
        # antes de procesarlo, en caso de que no tengamos permisos 
        # de escritura en el directorio de trabajo actual 
        # (que puede suceder cuando se implementa).
        tempReport <- file.path(tempdir(), "OncoReporter.Rmd")
        file.copy("OncoReporter.Rmd", tempReport, overwrite = TRUE)
        
        # Configuramos los parámetros a pasar al documento Rmd
        params <- list(Sexo: input$sexo,
                       Edad: input$edad,
                       Tipo_de_tumor: input$tumor,
                       Localización_del_tumor: input$localización,
                       Porcentaje_de_células_tumorales: input$céls_tum,
                       SNV: 'r filtered_data_snv',
                       CNV: 'r filtered_data_cnv')
        
        # Genera el documento, pasa la lista `params` y lo evalúa 
        # en un elemento secundario del entorno global (esto aísla 
        # el código del documento del código de esta aplicación).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)