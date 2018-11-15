library(shiny)
library(readxl)
library(stringr)
library(stringi)
library(DT)

source("cytoscript_vinput.R")

# Define server logic required to plot various variables against mpg

  function(input, output) {
    
    ##check file loaded
    getData <- reactive({
      if(is.null(input$file)){return(NULL)}else{return(TRUE)}
    })
    output$fileUploaded <- reactive({
      return(!is.null(getData()))
    })
    outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
    
    ## output$read=renderText({datasetInput()})
    
    # You can access the value of the widget with input$file, e.g.
    CytotableString <- reactive({
      cyto=""
      
      if(!is.null(input$text))
      {
        cyto<-input$text
      }
      ##if(!is.null(input$file))
      ##{
      ##  req(input$file)
      ##  cyto<-sapply(as.data.frame(read.delim(input$file$datapath,header=FALSE,sep="\t")),as.character)
      ##}
      CytoConverter(cyto)
      })
    
    CytotableFile<- reactive({
      cyto=""

      if(!is.null(input$file))
      {
        req(input$file)
        ##cyto<-readLines(input$file$datapath,encoding = 'UTF-8')
        cyto<-sapply(as.data.frame(read.delim(input$file$datapath,header=FALSE,sep="\t",fileEncoding = "UTF-8")),as.character)
      }
      CytoConverter(cyto)
    })
    
    ##select data to download
    datasetInput <- reactive({
      switch(input$dataset,"Both"=list(rbind(CytotableString()[[1]],CytotableFile()[[1]]),rbind(CytotableString()[[2]],CytotableFile()[[2]])),
             "Text"= CytotableString(),
             "File"= CytotableFile())
    })
    
    
    output$tableresult<-reactive({table_Result()})
    ##download data   table
      output$downloadData <- downloadHandler(
          filename = function() {
             paste('CytoConverter_Results', '.txt', sep='')
          },
           content = function(file) {
             write.table(datasetInput()[[1]], file,quote=F,sep='\t',row.names = F)
          }
        
      )
       
      output$downloadErrorLog <- downloadHandler(
        filename = function() {
          paste('CytoConverter_Error_Log' ,'.txt', sep='')
        },
        content = function(file) {
          write.table(datasetInput()[[2]], file,quote=F,sep='\t',row.names = F)
        }
        
      )
      
      ##download script button
      output$downloadScript <- downloadHandler(
        filename =
          paste('CytoConverter', '.zip', sep='')
        ,
        content = function(fname) {
          files2zip <- dir('Cytogenetic_software', full.names = TRUE)
          zip(zipfile = fname, files = files2zip)}, contentType = "application/zip"
        
      )
      
      output$table<-renderTable({CytotableString()[[1]]})
      output$filetable<-renderTable({CytotableFile()[[1]]}) 
     
      output$tablelog<-renderTable({CytotableString()[[2]]})
      output$filetablelog<-renderTable({CytotableFile()[[2]]})
      
      ##make download button work
     ## output$downloadData <- downloadHandler(
      ##filename = function() {
      ##  paste("sample", ".txt", sep = "")
      ##},
      ##content = function(file) {
        ##write.table(output$table, file, row.names = FALSE,sep="\t")
      ##}
      ##)
  }

##reactive(
##  if(!is.null(input$text))
##  {
##    cyto<-input$text 
##  }
##  if(!is.null(input$file))
##  {
##    cyto<-input$file
##  }
##)