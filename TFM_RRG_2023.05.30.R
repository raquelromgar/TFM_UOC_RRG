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

# Agregamos la definición de las variables filtered_data_snv y filtered_data_cnv 
# en el entorno global (antes de la función server) para que puedan ser 
# accedidas por otras partes del código.
filtered_data_snv <- NULL
filtered_data_cnv <- NULL
clasif_variantes <- NULL


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
               fileInput("SNV",label = "Fichero para SNVs, MNVs e INDELs", 
                         multiple = F, accept = c(".tsv")),
               fileInput("CNV",label = "Fichero para CNVs", 
                         multiple = F, accept = c(".tsv")),
               br(),
        ),
        
        column(6,
               h4("Datos"),
               selectInput("sexo", label = "Selecciona el sexo",
                           choices = c("Masculino" = "M", "Femenino" = "F")),
               textInput("edad", label = "Edad"),
               textInput("tumor", label = "Tipo de tumor"),
               textInput("localizacion", label = "Localización"),
               textInput("cels_tum", label = "Porcentaje de células tumorales"),
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
    ), 
    tabPanel(
      title = "Resultados CNVs",
      br(),
      #Tablas Resultados
      DT::dataTableOutput("filtered_table_CNV"),
      br(),
    ), 
    tabPanel(
      title = "Clasificación de las variantes",
      br(),
      #Tablas Clasificación
      DT::dataTableOutput("classification"),
      br(),
      textOutput("no_classification"),
      actionButton("generar_informe", "Generar informe"),
      downloadButton("report", "Descargar Informe"),
      br(),
    ),
  ),
  # Agregar el footer
  tags$footer(
    style = "text-align: center; margin-top: 20px; color: cornflowerblue;
    font-family: Calibri;",
    "Aplicación web desarrollada por Raquel Romero García para el TFM titulado 
    'Generación de informes genéticos para la toma de decisiones oncológicas' 
    del Máster de Bioinformática y Bioestadística de la UOC."
  )
)


############################################################
##    Definir la lógica del servidor para la aplicación   ##
############################################################     

# Definimos de las funciones filter_cnv_data y filter_snv_data 
# fuera de la función server, para que estén disponibles en el entorno global.

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
  
  
  # Observador para leer el archivo y aplicar los filtros cuando se pulsa 
  # el botón "submit"
  observeEvent(input$submit, {
    snv_data <- read_snv_file() # Leer el archivo TSV
    filtered_data_snv <- filter_snv_data(snv_data) # Aplicar filtros de SNV
    cnv_data <- read_cnv_file() # Leer el archivo TSV
    filtered_data_cnv <- filter_cnv_data(cnv_data) # Aplicar filtros de CNV
    
    # Cambiar a la página de resultados cuando se puse el botón submit
    updateTabsetPanel(session, "tablas", "Resultados SNVs, MNVs e INDELs")
    
    # Mostrar las tablas con los valores filtrados
    output$filtered_table_SNV <- DT::renderDataTable(server = FALSE, {
      datatable(
        filtered_data_snv,
        extensions = c('Buttons', 'FixedHeader'),
        caption = htmltools::tags$caption(
          style = 'caption-side: bottom; text-align: center;',
          'Tabla 1: Variantes de tipo SNV, MNV o INDEL', 
        ),
        # Opciones adicionales de la tabla
        options = list(
          dom = 'lBfrtip',
          fixedHeader = TRUE,
          buttons = list('print', 'copy',
                         list(extend = 'collection', 
                              buttons = c('excel'), 
                              text = 'Download')),
          columnDefs = list(list(
          targets = c(42,58,59),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 60 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 60) + 
            '...</span>' : data;",
            "}")
        )),
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#5dadbd', 
          'color': '#fff'});",
          "}")
        ), 
          callback = JS('table.page().draw(false);')
   ) 
})

    output$filtered_table_CNV <- DT::renderDataTable({
      datatable(
        filtered_data_cnv,
        extensions = c('Buttons', 'FixedHeader'),
        caption = htmltools::tags$caption(
          style = 'caption-side: bottom; text-align: center;',
          'Tabla 2: Variantes de tipo CNV', 
        ),
        # Opciones adicionales de la tabla
        options = list(
          dom = 'lBfrtip',
          fixedHeader = TRUE,
          buttons = list('print', 'copy',
                         list(extend = 'collection', 
                              buttons = c('excel'), 
                              text = 'Download')),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#5dadbd',
            'color': '#fff'});",
            "}")
        ), 
        callback = JS('table.page().draw(false);')
      ) 
    })
    
    #Función de clasificación
    tabla<-data.frame(Gen = character(),
                        cDNA = character(),
                        dbSNP = character(),
                        url = character(),
                        Clasificación = character(),
                        stringsAsFactors = FALSE)
    
    
    classif <- function(tabla, df) {
      for (i in 1:nrow(df)) {
        Gen <- df$Genes[i]
        cDNA <- paste0(df$Transcript[i], ":", df$Coding[i])
        dbSNP <- df$dbSNP[i]
        url <- ""
        if (!is.na(dbSNP[i])) {
          url <- sprintf('<a href="https://varsome.com/variant/hg38/%s" 
                     target="_blank">https://varsome.com/variant/hg38/%s</a>',
                         dbSNP[i], dbSNP[i])
        }
        Clasificación <- ""
        
        nuevaFila <- data.frame(Gen = Gen, cDNA = cDNA, dbSNP = dbSNP, 
                                url = url, Clasificación = Clasificación,
                                stringsAsFactors = FALSE)
        
        tabla <- rbind(tabla, nuevaFila)
      }
      
      return(tabla)
    }
    
    
    if (nrow(filtered_data_snv) != 0) {
      clasif_variantes <- classif(tabla, filtered_data_snv)
      
    # Agregar una columna con desplegable de opciones a la tabla
      output$classification <- DT::renderDataTable({
        datatable(
          clasif_variantes,
          escape = FALSE,
          caption = htmltools::tags$caption(
            style = 'caption-side: bottom; text-align: center;',
            'Tabla 3: Tabla de clasificación de variantes'
          ),
          # Opciones adicionales de la tabla
          options = list(
            dom = 'lBfrtip',
            fixedHeader = TRUE,
            columnDefs = list(
              list(
                targets = c(5),
                render = JS(
                  "function(data, type, row, meta) {",
                  "if (type === 'display' && data === '') {",
                  "return '<select><option value=\"Sin clasificar\">Sin clasificar</option>' +",
                  "'<option value=\"Benigna\">Benigna</option>' +",
                  "'<option value=\"Probablemente benigna\">Probablemente benigna</option>' +",
                  "'<option value=\"De significado incierto\">De significado incierto</option>' +",
                  "'<option value=\"Probablemente patogénica\">Probablemente patogénica</option>' +",
                  "'<option value=\"Patogénica\">Patogénica</option></select>';",
                  "} else {",
                  "return data;",
                  "}",
                  "}"
                )
              )
            ),
            initComplete = JS(
              "function(settings, json) {",
              "$(this.api().table().header()).css({'background-color': 
              '#5dadbd', 'color': '#fff'});",
              "}"
            )
          )
        )
      })
    }else{
      output$no_classification <- renderText("Ninguna variante que clasificar")
  }
 
  })
    
  
  # Enlace de descarga para el informe generado
  output$report <- downloadHandler(
    filename = "Informe de variantes.pdf",
    content = function(file) {
      params <- list(sexo = input$sexo,
                     edad = input$edad,
                     tumor = input$tumor,
                     localizacion = input$localizacion,
                     cels_tum = input$cels_tum,
                     filtered_data_snv = filtered_data_snv,
                     filtered_data_cnv = filtered_data_cnv,
                     clasif_variantes = clasif_variantes
                     )
      id <- showNotification(
        "Rendering report...", 
        duration = NULL, 
        closeButton = FALSE
      )
      on.exit(removeNotification(id), add = TRUE)
      
      rmarkdown::render("OncoReporter.Rmd", 
                        output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)