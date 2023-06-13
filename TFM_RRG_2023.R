###############################################################################
###############################################################################
##                                                                           ##
## Esta es una aplicación web Shiny. Puede ejecutar la aplicación haciendo   ##
##           clic en el botón 'Ejecutar aplicación' de arriba.               ##
##                                                                           ##  
###############################################################################
###############################################################################


# Cargamos las librerias necesarias

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(knitr)
library(rmarkdown)



# Copia el informe en un directorio temporal. Esto es muy importante al implementar
# la aplicación, ya que a menudo no se podrá escribir en el directorio de trabajo.

report_path <- tempfile(fileext = ".Rmd")
file.copy("OncoReporter.Rmd", report_path, overwrite = TRUE)

render_report <- function(input, output, params) {
  rmarkdown::render(input,
                    output_file = output,
                    params = params,
                    envir = new.env(parent = globalenv())
  )
}


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
               br(),
               fileInput("SNV",label = "Fichero para SNVs, MNVs e INDELs", 
                         multiple = F, accept = c(".tsv")),
               br(),
               fileInput("CNV",label = "Fichero para CNVs", 
                         multiple = F, accept = c(".tsv")),
               br(),
        ),
        
        column(6,
               h4("Datos"),
               selectInput("sexo", label = "Selecciona el sexo",
                           choices = c("Masculino", "Femenino")),
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
      br(),
      br(),
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

# Definimos de las funciones "filter_cnv_data", "filter_snv_data" y "classif"
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

# Función de clasificación
classif <- function(tabla, df) {
  for (i in 1:nrow(df)) {
    Gen <- df$Genes[i]
    Tipo<- df$Type[i]
    Frec.Alélica <- df$`Allele Frequency %`[i]
    Nºcopias <- df$`Copy Number`[i]
    cDNA <- paste0(df$Transcript[i], ":", df$Coding[i])
    AA <- df$`Amino Acid Change`[i]
    Clasificacion <- ""
    
    nuevaFila <- data.frame(Gen = Gen, Tipo = Tipo, Frec.Alélica = Frec.Alélica, 
                            Nºcopias = Nºcopias, cDNA = cDNA, AA = AA,
                            Clasificacion = Clasificacion, 
                            stringsAsFactors = FALSE)
    
    tabla <- rbind(tabla, nuevaFila)
  }
  
  return(tabla)
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
  
  datos_clasificados <- reactiveVal(NULL)
  
  # Observador para leer el archivo y aplicar los filtros cuando se pulsa el botón "submit"
  observeEvent(input$submit, {
    snv_data <- read_snv_file() # Leer el archivo TSV
    filtered_data_snv <- filter_snv_data(snv_data) # Aplicar filtros de SNV
    cnv_data <- read_cnv_file() # Leer el archivo TSV
    filtered_data_cnv <- filter_cnv_data(cnv_data) # Aplicar filtros de CNV
    
    # Cambiar a la página de resultados cuando se pulse el botón submit
    updateTabsetPanel(session, "tablas", "Resultados SNVs, MNVs e INDELs")
    
    # Mostrar las tablas con los valores filtrados
    output$filtered_table_SNV <- DT::renderDataTable(server = FALSE, {
      datatable(
        filtered_data_snv,
        extensions = c('Buttons', 'FixedHeader'),
        caption = htmltools::tags$caption(
          style = 'caption-side: bottom; text-align: left;',
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
          style = 'caption-side: bottom; text-align: left;',
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
    
    # Crear tabla de variantes para la clasificación
   crear_tabla <- function(filtered_data_snv, filtered_data_cnv) {
     tabla <- data.frame(Gen = character(),
                         Tipo = character(),
                         Frec.Alélica = double(),
                         Nºcopias = double(),
                         cDNA = character(),
                         AA = character(),
                         Clasificacion = character(),
                         stringsAsFactors = FALSE)
     
     if (nrow(filtered_data_snv) != 0) {
       clasif_variantes <- classif(tabla, filtered_data_snv)
     }else{
       clasif_variantes<-data.frame()
     }
     
     if (nrow(filtered_data_cnv) != 0) {
       clasif_variantes <- rbind(clasif_variantes, classif(tabla, filtered_data_cnv))
       }
     
     if (nrow(filtered_data_snv) == 0 && nrow(filtered_data_cnv) == 0) {
       output$no_classification <- renderText({
         texto_estilizado <- paste(toupper("Ninguna variante que clasificar"))
         HTML(texto_estilizado)
       })
     }
     
     return(clasif_variantes)
   }
   
   datos_clasificados(crear_tabla(filtered_data_snv, filtered_data_cnv))
   
   # Agregar una columna con desplegable de opciones a la tabla
   output$classification <- DT::renderDataTable({
     datos <- datos_clasificados()  # Obtener los datos originales
     
     # Agregar una columna adicional con HTML personalizado para el desplegable
     datos$Clasificacion <- lapply(datos$Clasificacion, function(option) {
       select <- paste0("<select>",
                        "<option value='Sin clasificar'", ifelse(option == "Sin clasificar", " selected", ""), ">Sin clasificar</option>",
                        "<option value='Benigna'", ifelse(option == "Benigna", " selected", ""), ">Benigna</option>",
                        "<option value='Probablemente benigna'", ifelse(option == "Probablemente benigna", " selected", ""), ">Probablemente benigna</option>",
                        "<option value='De significado incierto'", ifelse(option == "De significado incierto", " selected", ""), ">De significado incierto</option>",
                        "<option value='Probablemente patogénica'", ifelse(option == "Probablemente patogénica", " selected", ""), ">Probablemente patogénica</option>",
                        "<option value='Patogénica'", ifelse(option == "Patogénica", " selected", ""), ">Patogénica</option>",
                        "</select>")
       paste0(select)
     })
     
     datatable(
       datos,
       escape = FALSE,
       caption = htmltools::tags$caption(
         style = 'caption-side: bottom; text-align: center;',
         'Tabla 3: Tabla de clasificación de variantes'
       ),
       options = list(
         dom = 'lBfrtip',
         fixedHeader = TRUE,
         initComplete = JS(
           "function(settings, json) {",
           "$(this.api().table().header()).css({'background-color': '#5dadbd', 'color': '#fff'});",
           "}"
         ),
         drawCallback = JS(
           "function(settings) {",
           "$(this.api().table().header()).css({'background-color': '#5dadbd', 'color': '#fff'});",
           "}"
         )
       )
     )
   })
   
   
   
   
   
   
   
   
  })
    
  
  # Descargar el informe de variantes generado
  
  output$report <- downloadHandler(
    filename = "Informe de variantes.html",
    content = function(file) {
      params <- list(sexo = input$sexo,
                     edad = input$edad,
                     tumor = input$tumor,
                     localizacion = input$localizacion,
                     cels_tum = input$cels_tum,
                     datos_clasificados = datos_clasificados())
      callr::r(
        render_report,
        list(input = report_path, output = file, params = params)
      )
    }
  )
  

}

# Corremos la aplicación
shinyApp(ui = ui, server = server)