library(shiny)
library(stringr)
library(stringi)
##library(DT)

source("cytoscript_vinput.R")
source("cyto_graph.R")
source("plot_cyto_graph.R")
example_table<-read.delim(file="cyto_example.txt",sep='\t',header=F)
colnames(example_table)<-NULL
example_result_table<-read.delim(file="cyto_result.txt",sep='\t',header=T)
# Define server logic required to plot various variables against mpg

  function(input, output) {
    
    ##variable for hg type
  
    
    
    
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
      CytoConverter(cyto,input$radio)
      })
    
    CytotableFile<- reactive({
      cyto=""

      if(!is.null(input$file))
      {
        req(input$file)
        cyto<-readLines(input$file$datapath,encoding = 'UTF-8')
        cyto<-sapply(as.data.frame(read.delim(input$file$datapath,header=FALSE,sep="\t")),as.character)
        cyto<-({as.matrix(cyto)})
      }
      
      CytoConverter(cyto,input$radio)
    })
    
    ##return booleaen on length of error log
    output$nonemptyerrorlogstring<-reactive({if(nrow(CytotableString()[[2]])>0){return(TRUE)}else{return(FALSE)}})
    output$nonemptyerrorlogfile<-reactive({if(nrow(CytotableFile()[[2]])>0){return(TRUE)}else{return(FALSE)}})
    outputOptions(output, 'nonemptyerrorlogstring', suspendWhenHidden=FALSE)
    outputOptions(output, 'nonemptyerrorlogfile', suspendWhenHidden=FALSE)
    
    
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
      
      output$downloadExample <- downloadHandler(
        filename = function() {
          paste('CytoConverter_Example_File' ,'.txt', sep='')
        },
        content = function(file) {
          write.table(example_table,file=file,quote=F,sep='\t',row.names = F,fileEncoding = "UTF-8")
        }
        
      )
      
      output$downloadResultExample <- downloadHandler(
        filename = function() {
          paste('CytoConverter_Result_Example_File' ,'.txt', sep='')
        },
        content = function(file) {
          write.table(example_result_table,file=file,quote=F,sep='\t',row.names = F,fileEncoding = "UTF-8")
        }
        
        
      )
      
      ##function to plot stuff-strings
      make_plot_string<-reactive({ 
        cyto_graph(CytotableString()[[1]],input$radio)
        
      })
      
      ##function to plot stuff-table
      make_plot_table<-reactive({ 
        cyto_graph(CytotableFile()[[1]],input$radio)
        
      })        
  
      output$func<-reactive({str(make_plot_string())})
      output$table<-renderTable({CytotableString()[[1]]})
      output$filetable<-renderTable({CytotableFile()[[1]]}) 
      
      
      ##change this to recive stuff from make plot string and make plot table
      output$singleplot<-renderPlot({
          
      
        
        if(!is.null(CytotableString()[[1]]))
        {
          list_from_cyto<-cyto_graph(CytotableString()[[1]],input$radio)
          rect_maker<-list_from_cyto[[1]]
          xbegin<-list_from_cyto[[2]]
          xcoord_master<-list_from_cyto[[3]]
          y_above<-list_from_cyto[[4]]
          y_below<-list_from_cyto[[5]]
          sorted_reflist<-list_from_cyto[[6]]
          cum_length_coords<-list_from_cyto[[7]]
          start_cum_length<-list_from_cyto[[8]]
          uniq_coord_name<-list_from_cyto[[9]]  
        }
        
        
        
            if(length(uniq_coord_name)<50){
              ylabel=T
            }else{
              ylabel=F
            }
          

        
        
        plot.new()
        
        plot.window(c(0,1),c(0,1),mar=rep(0,4))
        
        rect(xleft=xcoord_master, xright=1,ybottom=y_below,ytop=y_above,col="gray90")
        
        if(nrow(rect_maker)>0)
        {
          apply(rect_maker,1,
                function(x){ 
                  if(x[5]=="Gain")
                  {
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col="red",border=NA)
                    
                  }else if(x[5]=="Loss"){
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col=rgb(0,0,1,alpha=0.5),border=NA)
                    
                  }else if(x[5]=="Double"){
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col="orange",border=NA)
                    
                  }
                  
                })
          
          
          ##take lines away if number of samples is over 20
          if(length(uniq_coord_name)<20)
          {
            sapply(unique((rect_maker[,3])),function(x){lines(x=c(xcoord_master,1),y=c(x[1],x[1]))})
          }
          
          if(ylabel)
          {
            text(x=xcoord_master,y=c(unique((rect_maker[,3]+rect_maker[,4])/2)),labels=uniq_coord_name,cex=0.8,pos=2)
          }
        }
        
        lines(x=c(xcoord_master,xcoord_master),y=c(y_above+0.02,y_above+0.08))
        sapply(cum_length_coords*(1-xcoord_master)+xcoord_master,function(x){lines(x=c(x,x),y=c(y_above+0.02,y_above+0.08));lines(x=c(x,x),y=c(y_above,y_below),col="white")})
        lines(x=c(xcoord_master,1),y=c(y_above+0.02,y_above+0.02))
        text(x=c(start_cum_length*(1-xcoord_master)+xcoord_master+as.numeric(sorted_reflist[,2])/sum(as.numeric(sorted_reflist[,2]))/2*(1-xcoord_master)),y=(y_above*2+0.1)/2,labels=gsub("chr","",sorted_reflist[,1]),cex=1,offset=0)
        
        
        rect(xleft=xcoord_master, xright=1,ybottom=y_below,ytop=y_above,col=NA)
        
        ##legend(x=0.20,y= 0.2,uniq_coord_name)
        legend(x=1,y= y_below-0.03,legend=c("Gain","Hemizygous Loss","Homozygous Loss"),fill=c("red",rgb(0,0,1,alpha=0.5),"orange"),xjust=1,yjust=1,cex=0.7)
        ##xlab("Chromosome")
        ##ylab("Sample")
        
          },height=500,width=1200)
     
      
      
      output$tableplot<-renderPlot({
        
        
        
        if(!is.null(CytotableFile()[[1]]))
        {
          list_from_cyto<-cyto_graph(CytotableFile()[[1]],input$radio)
          rect_maker<-list_from_cyto[[1]]
          xbegin<-list_from_cyto[[2]]
          xcoord_master<-list_from_cyto[[3]]
          y_above<-list_from_cyto[[4]]
          y_below<-list_from_cyto[[5]]
          sorted_reflist<-list_from_cyto[[6]]
          cum_length_coords<-list_from_cyto[[7]]
          start_cum_length<-list_from_cyto[[8]]
          uniq_coord_name<-list_from_cyto[[9]]  
        }
        
        
        
        if(length(uniq_coord_name)<50){
          ylabel=T
        }else{
          ylabel=F
        }
        
        
        
        
        plot.new()
        
        plot.window(c(0,1),c(0,1),mar=rep(0,4))
        
        rect(xleft=xcoord_master, xright=1,ybottom=y_below,ytop=y_above,col="gray90")
        
        if(nrow(rect_maker)>0)
        {
          apply(rect_maker,1,
                function(x){ 
                  if(x[5]=="Gain")
                  {
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col="red",border=NA)
                    
                  }else if(x[5]=="Loss"){
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col=rgb(0,0,1,alpha=0.5),border=NA)
                    
                  }else if(x[5]=="Double"){
                    rect(xleft=x[1], xright=x[2],ybottom=x[3],ytop=x[4],col="orange",border=NA)
                    
                  }
                  
                })
          
          
          ##take lines away if number of samples is over 20
          if(length(uniq_coord_name)<20)
          {
            sapply(unique((rect_maker[,3])),function(x){lines(x=c(xcoord_master,1),y=c(x[1],x[1]))})
          }
          
          if(ylabel)
          {
            text(x=xcoord_master,y=c(unique((rect_maker[,3]+rect_maker[,4])/2)),labels=uniq_coord_name,cex=0.8,pos=2)
          }
        }
        
        lines(x=c(xcoord_master,xcoord_master),y=c(y_above+0.02,y_above+0.08))
        sapply(cum_length_coords*(1-xcoord_master)+xcoord_master,function(x){lines(x=c(x,x),y=c(y_above+0.02,y_above+0.08));lines(x=c(x,x),y=c(y_above,y_below),col="white")})
        lines(x=c(xcoord_master,1),y=c(y_above+0.02,y_above+0.02))
        text(x=c(start_cum_length*(1-xcoord_master)+xcoord_master+as.numeric(sorted_reflist[,2])/sum(as.numeric(sorted_reflist[,2]))/2*(1-xcoord_master)),y=(y_above*2+0.1)/2,labels=gsub("chr","",sorted_reflist[,1]),cex=1,offset=0)
        
        
        rect(xleft=xcoord_master, xright=1,ybottom=y_below,ytop=y_above,col=NA)
        
        ##legend(x=0.20,y= 0.2,uniq_coord_name)
        legend(x=1,y= y_below-0.03,legend=c("Gain","Hemizygous Loss","Homozygous Loss"),fill=c("red",rgb(0,0,1,alpha=0.5),"orange"),xjust=1,yjust=1,cex=0.7)
        ##xlab("Chromosome")
        ##ylab("Sample")
      },height=800,width=1200)
      
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