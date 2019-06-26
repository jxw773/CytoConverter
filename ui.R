navbarPage( "CytoConverter : Cytogenetic Nomenclature to Genomic Coordinate Translator",
tabPanel("Main Page",
  fluidPage(
  tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    ##titlePanel("Cytogenetics Nomenclature to Genomic Coordinate Translator"),
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        ##select file to download
        selectInput("dataset", "Display:",
                    choices = c("Both","Text","File")),
  
        textInput("text", label = h3("Text input"), value = ""),
  
        fileInput("file", label = h3("File input"),accept=""),
       
        
        ##buttons for reference type
        radioButtons("radio", label = h3("Reference Build"),
                     choices = list("GRCh38" = "GRCh38", "hg19" = "hg19", "hg18" = "hg18"), 
                     selected = "GRCh38"),
        
        hr(),
        
        fluidRow(column(3, verbatimTextOutput("value"))),
  
        ##download data button
       
        downloadButton("downloadData", "Download Results"),
       br(),
       br(),
        downloadButton("downloadErrorLog", "Download Error Log"),
       br(),
       br(),
       hr(),
       h3("Instructions"),
       "1. Program accepts a single karyotype in the 'Text input' box or a tab-delimited UTF-8 encoded text file with sample name in column 1 and karyotype in column 2 in the 'File input' box. An example of a file input and output, in UTF-8 encoding, is shown below:",
       
       br(),
       br(),
       "ABC\t46,XY,der(1;19)(q10;p10)",
       br(),
       "DEF\t47,XX,+der(10)t(10;21)(p13;q21)",
       br(),
       "P10\t47,X,+X[30]/48,XX,+7,+9[50]",
       
       br(),
       br(),
       downloadButton("downloadExample", "Download File Input Example"),
       downloadButton("downloadResultExample", "Download File Output Example"),
       br(),
       br(),
       "2. Program returns a text tab-delimited file that can be downloaded using the 'Download Results' button, and an error log that can be viewed under the tab 'Error Log' and downloaded using the 'Download Error Log' button. The results are presented in a table with the sample name (adding '_n' to the sample name for gains/losses present in the nth clone in the karyotype), beginning coordinate, ending coordinate, whether the region was lost or gained, and number of cells present of a clone out of the total number of cells in a sample with that gain or loss, if that information is provided in the input karyotype. The sample name is given as 'sample' if the sample name is not provided. The error log contains the sample name_n, karyotype, and the issue or warning encountered while parsing the karyotype." ,
       
       br(),
       br(),
       "3. Please be patient with results generated from large input files-it may take while after file is uploaded for table to be generated.",
       br(),
       br(),
       "4. Check error log after table has been generated.",
       br(),
       br(),
       "5. Application will disconnect if inactive for a period of time.",
       br(),
       br(),
       "6. If the user wishes to download the script for use on their own machine, the download link for the program is provided",
       a("here.",href="https://github.com/jxw773/CytoConverter")
       

       ),
      
      mainPanel (

        ##text output
        conditionalPanel(
        condition = "(input.dataset == 'Text'||input.dataset == 'Both') && input.text != ''",
        headerPanel("Text Result"),
        br(),
        ##indicates error in text output
        conditionalPanel(
          condition = "output.nonemptyerrorlogstring",
          br(),
          p("Please Check Error Log",style="color:red"),
          br()
        ),
        
        tableOutput("table"),
       
        plotOutput("singleplot", height="500px",width="1200px",
                               click = "plot_click",  # Equiv, to click=clickOpts(id="plot_click")
                               hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                               brush = brushOpts(id = "plot_brush"))
        
        
        ),
    
        

        ##table output
        conditionalPanel(
        condition = "(input.dataset == 'File'||input.dataset == 'Both') && output.fileUploaded",
        
        headerPanel("File Result"),
        
        ##indicates error in table output
       conditionalPanel(
          condition = "output.nonemptyerrorlogfile" ,
          br(),
          p("Please Check Error Log",style="color:red"),
          br()
        ),
        
        tableOutput("filetable"),
        plotOutput("tableplot", height="800px",width="1200px",
                   click = "plot_click",  # Equiv, to click=clickOpts(id="plot_click")
                   hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                   brush = brushOpts(id = "plot_brush")),
        value="File Result"
        )
        
      )
      )
    
    
      
    )
  ),
##something about here is messing code up 
##error log
tabPanel("Error Log",

           mainPanel(
             conditionalPanel(
               condition = "(input.dataset == 'Text'||input.dataset == 'Both') && input.text != ''",
               headerPanel("Text Error Log"),
               tableOutput("tablelog")
             ),
             
             conditionalPanel(
               condition = "(input.dataset == 'File'||input.dataset == 'Both') && output.fileUploaded",
               headerPanel("File Error Log"),
               tableOutput("filetablelog")
             )##,
             
             ## conditionalPanel(
              ## condition ="true",
              ##  textOutput("fileUploaded"),
              ##  textOutput("nonemptyerrorlogfile"),
              ##  textOutput("nonemptyerrorlogstring")
              ##)
             
           )

         
),
##help page
tabPanel("About",
         mainPanel(
           h3("Description"),
           "Program extracts net gain and loss of material relative to a normal karyotype from cytogenetic nomenclature into genomic coordinates. Nomenclature must conform to the ISCN 2016. Genomic coordinates output is according to build Gchr38. ",
           h3("Instructions"), 
           "Paste karyotype into text input",
           br(),
           "or",
           br(),
           "Upload a tab delimited table with 2 columns: one with sample name and one with the associated karyotype into file input. If file returns an error, try changing encoding of file to UTF-8",
           br(),
           h3("Known Issues"),
           "Program does not parse chromosome abberations containing question marks. The counting of sex chromosomes is sometimes buggy. Program does not parse fish readings.	Derivative chromosomes involving translocations of the same chromosome number (or sex chromosome) will generate error. In the ISCN 2016, the relevant chromosome is denoted by underlining; however, underlining is not supported in plain text files or in the R. 
          The counting of chromosomes does not work well when question marks are present in the karyotypes or the karyotypes do not follow the ISCN. Haploidy or polyploidy is to be interpreted with caution under such circumstances.
           CytoConverter will not processes karyotypes such as 48,XX, del(8)(p11),dup(8)(q20)x2 properly because an overall count of aberrations per chromosomes is not kept as of now. Please input the karyotype as 48,XX, +del(8)(p11),dup(8)(q20)x2 instead 
           Program allows for some deviations from the ISCN 2016 to be processed correctly but please don't push it (e.g. 46XX, +7, using - in place of ~, is not case sensitive, errors in spacing are allowed). 
           ",
           h3("Contact Information"),
           "cytoconverter@gmail.com",
           h3("Download Script"),
          a("Script Repository",href="https://github.com/jxw773/CytoConverter")
           )
         )

)
