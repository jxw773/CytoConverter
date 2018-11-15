navbarPage("Cytogenetic Nomenclature to Genomic Coordinate Translator",
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
       
        
  
        ##download data button
       
        downloadButton("downloadData", "Download Results"),
       br(),
       br(),
        downloadButton("downloadErrorLog", "Download Error Log"),
       br(),
       br(),
       hr(),
       h3("Instructions"),
       "1. Program accepts tab-delimited text file with sample name in column 1 and karyotype in column 2 in file input or a karyotype in text input",
       br(),
       br(),
       "2. Program returns tab-delimited text file and error log that can be viewed under tab error log",
       br(),
       br(),
       "3. Please be patient with results generated from large input files-it may take while after file is uploaded for table to be generated",
       br(),
       br(),
       "4. Check error log",
       br(),
       br(),
       "5. Application will disconnect if inactive for an hour"

       ),
      
      mainPanel (
        conditionalPanel(
        condition = "(input.dataset == 'Text'||input.dataset == 'Both') && input.text != ''",
        headerPanel("Text Result"),
        tableOutput("table")
        ),
    
        conditionalPanel(
        condition = "(input.dataset == 'File'||input.dataset == 'Both') && output.fileUploaded",
        headerPanel("File Result"),
        tableOutput("filetable"), 
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
             )
             
             
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
           a("cytoconverter@gmail.com",href="mailto:cytoconverter@gmail.com"),
           h3("Download Script"),
           downloadButton("downloadScript", "Download Script")
           )
         )

)
