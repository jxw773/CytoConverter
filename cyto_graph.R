
#' Cytogenetic Graph Data Preparation
#'
#' @description
#' This function prepares cytogenetic data for visualization by processing CytoConverter 
#' results and organizing them for plotting. It handles coordinate transformation,
#' chromosome ordering, and data structure preparation for the plotting functions.
#' 
#' This function is primarily used internally by plot_cyto_graph() and should typically
#' not be called directly by end users.
#'
#' @param cyto_list Data frame containing CytoConverter results with chromosomal aberrations
#' @param ref_list Reference genome build specification. Options:
#'   \itemize{
#'     \item "GRCh38" - GRCh38/hg38 human genome build (default)
#'     \item "hg19" - hg19 human genome build  
#'     \item "hg18" - hg18 human genome build
#'     \item "hg17" - hg17 human genome build
#'     \item Custom cytoband matrix - User-provided cytoband data
#'   }
#' @param include_normals_graph Boolean flag to include normal samples for comparison (default: FALSE)
#' @param list_of_samples Vector of sample names to include when include_normals_graph is TRUE
#'
#' @return List containing processed data structures for plotting:
#'   \item{rect_maker}{Data frame with rectangle coordinates for plotting aberrations}
#'   \item{coordlist}{Processed chromosome coordinates}
#'   \item{y_coordlist}{Y-axis positioning data for samples}
#'   \item{uniq_coord_name}{Unique sample identifiers}
#'   
#' @details
#' The function processes cytogenetic data through several steps:
#' \itemize{
#'   \item Loads appropriate cytoband reference data based on ref_list parameter
#'   \item Converts chromosomal coordinates to plotting coordinates
#'   \item Handles X and Y chromosome naming (converted to 23 and 24 respectively)
#'   \item Organizes data for rectangular plotting regions representing aberrations
#'   \item Calculates cumulative chromosome lengths for continuous plotting
#' }
#'
#' @note This function expects cyto_list to contain columns for chromosome, start position,
#' end position, and aberration type. Column names should match CytoConverter output format.
#'
#' @examples
#' \dontrun{
#' # Typically called internally by plot_cyto_graph()
#' result <- CytoConverter(karyotype_data)
#' graph_data <- cyto_graph(result$Results, ref_list = "GRCh38")
#' }
#'
#' @seealso 
#' \code{\link{plot_cyto_graph}} for the main plotting function that uses this data
#' \code{\link{cyto_graph_fusion}} for fusion-specific data preparation
#'
##setting up blank plot
cyto_graph<-function(cyto_list,ref_list="GRCh38",include_normals_graph=F,list_of_samples=NULL){
  
  ##if it is string, set it equal to one of the stuff  
  if(length(ref_list)<=1)
  {
    if(ref_list =="GRCh38")
    {
      ref_list <-
        sapply(as.data.frame(
          read.delim("Builds/cytoBand_GRCh38.txt", header = FALSE)
        ), as.character)
    }else if(ref_list =="hg19"){
      ref_list <-
        sapply(as.data.frame(
          read.delim("Builds/cytoBand_hg19.txt", header = FALSE)
        ), as.character)
    }else if(ref_list =="hg18"){
      ref_list <-
        sapply(as.data.frame(
          read.delim("Builds/cytoBand_hg18.txt", header = FALSE)
        ), as.character)
    }else if(ref_list=="hg17"){
      ref_list <-
        sapply(as.data.frame(
          read.delim("Builds/cytoBand_hg17.txt", header = FALSE)
        ), as.character)
      
    }else if(is.null(ref_list))
    {
      ##default is grch38
      ref_list <-
        sapply(as.data.frame(
          read.delim("Builds/cytoBand_GRCh38.txt", header = FALSE)
        ), as.character) 
    }else{
      return("ref_list incorrectly specified")
    }
    
    
    ref_list <-as.data.frame(ref_list[sapply(unique(ref_list[,1]),function(x){grep(x,ref_list[,1])[length(grep(paste(x,"$",sep=""),ref_list[,1]))]}),][,c(1,3)])
    ref_list<-apply(ref_list,2,as.character)  
  }else{
    ref_list<-apply(ref_list,2,as.character)  
  }
  
  # Include normal samples in graph if requested (for comparison visualization)
  if(include_normals_graph){
    # Create dummy entries for normal samples with neutral aberrations
    temp<-cbind(list_of_samples,"chr1",0,0,"Gain",NA)
    
    # Handle single row case (convert to vector for proper processing)
    if(nrow(cyto_list)==1)
    {
      cyto_list<-as.vector(cyto_list) 
    }
    
    # Merge normal samples with aberration data
    if(nrow(cyto_list)>0)
    {
      colnames(temp)<-colnames(cyto_list)
      cyto_list<-rbind(temp,cyto_list)
    }else{
      cyto_list<-temp
      colnames(temp) <- c(
        "Sample ID", "Chr", "Start", "End", "Type", "Percent Present"
      )
    }
  }
  
  # Process aberration data for overlapping loss detection and visualization setup
  if( nrow(cyto_list) >= 1)
  {
    # Detect overlapping losses for "Double" classification
    # This algorithm identifies cases where the same sample has multiple 
    # overlapping loss regions, which should be marked as "Double" events
    double_loss<-cyto_list[which(cyto_list[,5]=="Loss"),]
    double_loss[,1]<-as.character(double_loss[,1])
  
    if(nrow(double_loss)>1)
    {
      loss_overlap<-data.frame()
      chr_list<-unique(double_loss[,2])  # Get unique chromosomes with losses
      name_list<-as.character(unique(double_loss[,1]))  # Get unique sample names
      
    # Iterate through each chromosome and sample combination
    for(i in 1:length(chr_list))
    {
          for(k in 1:length(name_list))
          {
            # Get all loss events for this chromosome-sample combination
            chr_table<-double_loss[intersect(which(double_loss[,2]==chr_list[i]) , which(double_loss[,1]==name_list[k] )),]
            
            # Process only if there are multiple loss events in the same sample/chromosome
            if(nrow(chr_table) >=2)
            {
              # Sort by chromosome position for overlap detection
              chr_table <- chr_table[order(chr_table[,2],chr_table[,3]),]
              
              # Compare each loss event with subsequent events for overlaps
              for(d in 1:(nrow(chr_table)-1))
              {
                overlap=F
                if(!is.na(chr_table[d,1])){
                  for(j in (d+1):(nrow(chr_table))){
                      if(!is.na(chr_table[j,1])){
                      
                        # Extract coordinate ranges for overlap testing
                        first <- as.numeric(chr_table[d,3:4])  # Start-End of first loss
                        sec <- as.numeric(chr_table[j,3:4])    # Start-End of second loss
                        
                        # Check for genomic overlap using %overlaps% function
                        if(first %overlaps% sec)
                        {
                          overlap=T
                          # Mark overlapping regions for "Double" classification
                         ## }
                         
                          
                          ##if(first[1] >  sec[2]){
                          ##  chr_table[d,2] <- sec[2]
                         ## }
                          if(first[2] >  sec[2]){
                            chr_table[d,4] <- sec[2]
                          }
                          
                          if(first[1] <  sec[1]){
                            chr_table[d,3] <- sec[1]
                          }
                          chr_table[j,] <- rep(NA,6)
                          
                        }
                    }
                  }
  
                }
                
                if(overlap)
                {
                  loss_overlap<-rbind(loss_overlap,chr_table[d,])
                }
  
              }
            }
            
      
          }
  
    }
      if(is.vector(loss_overlap)&& !is.na(loss_overlap[1]))
      {
        
        loss_overlap[5]<-"Double"
        cyto_list<-rbind(cyto_list,loss_overlap)
        
      }else if(nrow(loss_overlap)>0)
      {
        
        loss_overlap <- loss_overlap[which(!is.na(loss_overlap[,1])),]
        loss_overlap[,5]<-"Double"
        cyto_list<-rbind(cyto_list,loss_overlap)
      }  
    }
    ##if not assume ref list was inputted 
  }
  sorted_reflist<-ref_list[c(order(as.numeric(gsub("chr","",ref_list[1:22,1]))),23:24),]
  coords<-as.numeric(sorted_reflist[,2])
  ##length_coords<-coordss[2:length(coordss)]-coordss[1:(length(coordss)-1)]
  cum_length_coords=cumsum(coords/(sum(coords)))
  start_cum_length=c(0,cum_length_coords[1:length(cum_length_coords)-1])
  ##setting up lables
  
  ##x coord starting point
  xbegin=0
  ##default parameters 
  xcoord_master=0
  
  #determine how big graph should be
  ##if(length(uniq_coord_name)<=2)
  ##{
  ##  y_above=0.90
  ##  y_below=0.70
  ##}
  ##else if(length(uniq_coord_name)>10)
  ##{
  y_above=0.90
  y_below=0.15
  ##}else{
  ##  y_above=0.90
  ##  y_below=0.40
  ##}
  
  uniq_coord_name<-vector()
  rect_maker<-data.frame()
  matched_coord_names<-vector()
  y_coordlist<-matrix(ncol=2,nrow=0)
  
  if(!is.null(cyto_list) && nrow(cyto_list) > 0)
  {
    ##plotting data values
    coord_name<-cyto_list[,1]
    uniq_coord_name<-unique(coord_name)

    
    
    ##y coordinates (by name)

    if(length(uniq_coord_name)>1)
    {
      ##if(nrow(cyto_list)>1)
      ## {
      ##  cyto_list<-cyto_list[match(cyto_list[,1],uniq_coord_name),]
      ##}
      
      
     ##coord_name<-coord_name[match(coord_name,uniq_coord_name)]
      
      
      ##matched_coord_names<-cyto_list
      matched_coord_names<-lapply(uniq_coord_name,function(x){y=which(as.character(x)==as.character(coord_name));cbind(as.vector(y),rep(which(x==uniq_coord_name),length(y)))})
    
      }else{
      y=which(uniq_coord_name==coord_name)
      matched_coord_names<-cbind(as.vector(y),rep(1,length(y)))
    }
   
   ## if(is.list(matched_coord_names))
    ##{
   ##   temp_coord_matrix<-data.frame()
   ##   for(i in 1:length(matched_coord_names))
   ##   {
   ##     rbind(temp_coord_matrix,matched_coord_names[[i]])
    ##  }
  ##  }
  ##  
   
  # Build y-coordinate list from matched coordinate names
  # Handle both list and matrix formats for coordinate data
  if(is.list(matched_coord_names))
  {
    for(i in 1:length(matched_coord_names))
    {
      y_coordlist<-rbind(y_coordlist,matched_coord_names[[i]])
    }
  }else{
    y_coordlist<-matched_coord_names
  }
}
  # Calculate x-coordinate master position based on maximum sample name length
  # This determines where the chromosome plot area begins (after sample labels)
  xcoord_master=max(nchar(as.character(uniq_coord_name)))*0.005+xbegin
  
  
  # Sort cyto_list according to y_coordinate positions for proper vertical stacking
  if(length(y_coordlist)>0 && nrow(cyto_list)>1 )
  {
    cyto_list<-cyto_list[order(cyto_list[,1],y_coordlist[,1]),]

  }
  
  # Handle single aberration case (vector format)
  if(is.vector(cyto_list)){
    coords_listed<-cyto_list[2:5]  # Extract chromosome, start, end, type
    temp_coords<-gsub("chr","",coords_listed[1])  # Remove "chr" prefix for numeric processing
    
    # Convert sex chromosomes to numeric values for coordinate calculations
    temp_coords[grep("Y",temp_coords)]<-24  # Y chromosome = 24
    temp_coords[grep("X",temp_coords)]<-23  # X chromosome = 23
    coordlist<-as.numeric(temp_coords);
    
    # Apply same conversion to coordinate list for consistent indexing
    coords_listed[grep("X",coords_listed[,1]),1]<-23
    coords_listed[grep("Y",coords_listed[,1]),1]<-24      
    
    # Calculate y-axis coordinates for sample positioning
    # Distributes samples evenly across vertical plot space
    y_area_coord=cbind((y_coordlist[,2]-1)*-(y_above-y_below)/length(uniq_coord_name)+y_above,y_above-y_coordlist[,2]*(y_above-y_below)/length(uniq_coord_name))
    
    
    # Calculate rectangle x-coordinates based on genomic positions
    # Maps genomic coordinates to proportional positions on plot x-axis
    
    # X-axis start position: cumulative chromosome position + relative position within chromosome
    xstart<-(start_cum_length[coordlist]+(as.numeric(coords_listed[2]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
    # X-axis end position: similar calculation for end coordinate
    xend<-(start_cum_length[coordlist]+(as.numeric(coords_listed[3]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
    
    
    # Assemble all rectangle coordinates and aberration type information
    rect_maker<-as.data.frame(cbind(xstart,xend,y_area_coord,coords_listed[4]))
    rect_maker[1:4]<-apply(rect_maker[,1:4],2,function(x){as.numeric(as.character(x))})
    rect_maker[5]<-as.character(rect_maker[,5])
    
  }else{
    # Handle multiple aberrations case (matrix format)
    coords_listed<-cyto_list[,2:5]
    temp_coords<-gsub("chr","",coords_listed[,1])
      temp_coords[grep("Y",temp_coords)]<-24
      temp_coords[grep("X",temp_coords)]<-23
      coordlist<-as.numeric(temp_coords);
      
      ##adjust chrom name for x and y for input data
      coords_listed[grep("X",coords_listed[,1]),1]<-23
      coords_listed[grep("Y",coords_listed[,1]),1]<-24
      
      
      y_area_coord=cbind((y_coordlist[,2]-1)*-(y_above-y_below)/length(uniq_coord_name)+y_above,y_above-y_coordlist[,2]*(y_above-y_below)/length(uniq_coord_name))
      
      
      ##calculating where rectangle should start from how long the sample is
      
      
      
      ##x coordinates
      xstart<-(start_cum_length[coordlist]+(as.numeric(coords_listed[,2]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
      xend<-(start_cum_length[coordlist]+(as.numeric(coords_listed[,3]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
      
      
      ##all info for coordinates
      rect_maker<-as.data.frame(cbind(xstart,xend,y_area_coord,coords_listed[,4]))
      rect_maker[,1:4]<-apply(rect_maker[,1:4],2,function(x){as.numeric(as.character(x))})
      rect_maker[,5]<-as.character(rect_maker[,5])
      
    }
    
  
  sorted_reflist<-as.data.frame(sorted_reflist)
  sorted_reflist[,1]<-as.character(sorted_reflist[,1])
  sorted_reflist[,2]<-as.numeric(as.character(sorted_reflist[,2]))
  return(list(rect_maker,xbegin,xcoord_master,y_above,y_below,sorted_reflist,cum_length_coords,start_cum_length,uniq_coord_name))
}
