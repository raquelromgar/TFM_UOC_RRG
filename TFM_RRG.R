###############################################################################
###############################################################################
##                                                                           ##
## Esta es una aplicación web Shiny. Puede ejecutar la aplicación haciendo   ##
##           clic en el botón 'Ejecutar aplicación' de arriba.               ##
##                                                                           ##  
##                                                                           ##
###############################################################################
###############################################################################


library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(dplyr)



############################################################
##    Definir la interfaz de usuario para la aplicación   ##
############################################################

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "cosmo"),
  
  # add logout button UI
  div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
  # add login panel UI function
  shinyauthr::loginUI(id = "login"),
  
  titlePanel("OncoReporter"),
  
  tabsetPanel(id = "tablas",
    tabPanel(
      title = "Archivos y datos",
      fluidRow(
        h3("Sube tus archivos a la categoria correspondiente"),
      ),
      fluidRow(
        column(5,
               fileInput("SNV",label = "Fichero para SNVs, MNVs e INDELs", multiple = F, accept = c(".tsv")),
               fileInput("CNV",label = "Fichero para CNVs", multiple = F, accept = c(".tsv")),
               br(),
        ),
        
        column(6,
               h4("Datos"),
               selectInput("Sexo", label = "Selecciona el sexo",
                           choices = c("Masculino" = "M", "Femenino" = "F")),
               textInput("edad", label = "Edad"),
               textInput("tumor", label = "Tipo de tumor"),
               textInput("localizacion", label = "Localización"),
               textInput("%tum", label = "Porcentaje de células tumorales"),
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
      title = "Tratamientos",
      br(),
    )
  )
)
     

# Definar la lógica del servidor
server <- function(input, output, session) {
  
  
  # Función para leer el archivo TSV
  read_tsv_file <- function() {
    req(input$SNV, input$CNV) # Verificar si se ha subido un archivo
    
    # Leer el archivo TSV y retornar el dataframe
    snv_data <- read_tsv(input$SNV$datapath, comment = "#")
    return(snv_data)
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
    snv_data <- read_tsv_file() # Leer el archivo TSV
    filtered_data_snv <- filter_snv_data(snv_data) # Aplicar filtros de SNV
    cnv_data <- read_tsv_file() # Leer el archivo TSV
    filtered_data_cnv <- filter_cnv_data(cnv_data) # Aplicar filtros de CNV
    
    # Cambiar a la página de resultados cuando se puse el botón submit
    updateTabsetPanel(session, "tablas", "Resultados SNVs, MNVs e INDELs")
    
    # Mostrar la tabla con los valores filtrados
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
    
    #BOTONES DESCARGA 
    output$desc_res_SNV <- downloadHandler(
      filename = "Resultado_priorizacion_SNVs.tsv",
      content =  function(file){
        write_tsv(filtered_data_snv, file)
      },
      contentType = "text/plain"
    )
    
    output$desc_res_CNV <- downloadHandler(
      filename = "Resultado_priorizacion_CNVs.tsv",
      content =  function(file){
        write_tsv(filtered_data_cnv, file)
      },
      contentType = "text/plain"
    )
    
  })

  output$logoUOC <- renderImage({
    list(src = "Images/logoUOC2.png",contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
    
}

# Run the application
shinyApp(ui = ui, server = server)