###############################################################################
###############################################################################
##                                                                           ##
## This is a Shiny web application. You can run the application by clicking  ##
## the 'Run App' button above.                                               ##  
##                                                                           ##
###############################################################################
###############################################################################


library(shiny)
library(shinydashboard)
library(shinyauthr)
library(tidyverse)
library(DT)
library(dplyr)


####################################
##    Define UI for application   ##
####################################

ui <- fluidPage(
  
  # add logout button UI
  div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
  # add login panel UI function
  shinyauthr::loginUI(id = "login"),
  
  titlePanel("Onco_Reporter"),
  
  tabsetPanel(
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
actionButton("run", label= "Priorización de variantes")
)
),

    tabPanel(
      title = "Resultados",
      br(),
      #Descarga resultados
      downloadButton("desc_res_SNV", label = "Guardar resultado (SNVs, MNVs e INDELs)"),
      downloadButton("desc_res_CNV", label = "Guardar resultado (CNVs)"),
      br(), br(),
      #Tablas Resultados
      textOutput("res_SNV"),
      br(),
      dataTableOutput("resultado.SNV"),
      br(),
      textOutput("res_CNV"),
      br(),
      dataTableOutput("resultado.CNV"),
      br(),
    ),

    tabPanel(
      title = "Tratamientos",
      br(),
      )
)
)

      
# Define server logic for application ----

server <- function(input, output, session) {
    
  observeEvent(input$lectura, {
    {
    #SNVs, MNVs e INDELS
    if (is.null(input$SNV) == FALSE &
        is.null(input$CNV)==FALSE){
        withProgress(message = "Procesando los archivos...", value = 0, {
        #SNV
        cat("Inicio del procesamiento del archivo de SNVs, MNVs e INDELs", "\n")
        archivo_SNV<-input$SNV
        incProgress(0.01, detail = "Iniciando lectura y procesamiento del archivo de SNVs, MNVs e INDELs. Leyendo archivo TSV")
        tsvData <- read_tsv(archivo_SNV$datapath, comment = "#")
        cat("Archivo tsv correspondiente con SNVs, MNVs e INDELs leido","\n")
        incProgress(0.09, detail = "Archivo de SNVs, MNVs e INDELs leido.")
        
        #CNV
        cat("Inicio del procesamiento del archivo de CNVs", "\n")
        archivo_CNV<-input$CNV
        incProgress(0.01, detail = "Iniciando lectura y procesamiento del archivo de CNVs. Leyendo archivo TSV")
        tsvData <- read_tsv(archivo_CNV$datapath, comment = "#")
        cat("Archivo tsv correspondiente con CNVs leido","\n")
        incProgress(0.09, detail = "Archivo de CNVs leido.")
         })
       } 
    } 
    
##Filtros
observeEvent(input$run,{
  {
    output$res_SNV <- renderText("Resultado de la priorizaci?n para SNVs, MNVs e INDELS")
    res_SNV<- reactive({
      SNV_variants <- tsvData %>%
      filter(Filter == "PASS") %>%
      filter(!(Type %in% c("CNV", "EXPR_CONTROL", "REF", "ASSAYS_5P_3P"))) %>%
      filter(`Variant Effect` != "synonymous" | is.na(`Variant Effect`)) %>%
      filter(`UCSC Common SNPs`!="YES" | is.na(`UCSC Common SNPs`)) %>%
      filter(!grepl("benign", ClinVar, ignore.case=TRUE) | is.na(ClinVar)) %>%
      filter(`P-Value`<=0.05) %>%
      filter(grepl("exonic", Location, ignore.case=TRUE)) %>%
      filter(`Allele Frequency %`>=3) %>%
      filter(Coverage>=300)
      return(SNV_variants)})
    output$resultado.SNV<-DT::renderDataTable(DT::datatable({res_SNV()}))
    
    output$res_CNV <- renderText("Resultado de la priorizaci?n para CNVs")
    res_CNV<-reactive({
      CNV_variants <- tsvData %>%
      filter(`CNV P-Value`<=0.05) %>%
      filter(`Copy Number`>4 | `Copy Number`<=0.01)
      return(CNV_variants)})
    output$resultado.CNV<-DT::renderDataTable(DT::datatable({res_CNV()}))
  
}

#BOTONES DESCARGA 
output$desc_res_SNV <- downloadHandler(
  filename = "Resultado_priorizaci?n_SNVs.csv",
  content =  function(file){
  write.csv2(res_SNV() , file)
    },
  contentType = "text/plain"
  )
  
output$desc_res_CNV <- downloadHandler(
  filename = "Resultado_priorizaci?n_CNVs.csv",
  content =  function(file){
  write.csv2(res_CNV() , file)
    },
  contentType = "text/plain"
  )
  
  })
  })
      
     
output$logoUOC <- renderImage({
  list(src = "media/logoUOC2.png",contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
} 
  
  # Run the app ----
shinyApp(ui = ui, server = server)
