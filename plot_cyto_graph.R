#' Standard Cytogenetic Plotting Function
#' 
#' @description
#' This function creates standard visualizations for chromosomal gains and losses from 
#' CytoConverter analysis. It takes the output of cyto_graph and plots the information 
#' as a heatmap-style graph showing genomic aberrations across samples.
#' 
#' For structural rearrangements and fusion events, use cyto_graph_fusion() instead,
#' which provides specialized visualization capabilities for complex karyotypes.
#' 
#' @param cyto_list Data frame containing CytoConverter results with gains and losses
#' @param list_from_cyto Pre-computed output from cyto_graph (optional, unnecessary if cyto_list provided)
#' @param ref_list Reference genome build ("GRCh38", "hg19", "hg18", "hg17") 
#' @param ylabel Boolean flag to enable/disable sample name labels on y-axis (auto-determined if NULL)
#' @param include_normals_graph Boolean flag to include normal samples for comparison (default: FALSE)
#' @param list_of_samples Vector of normal sample names when include_normals_graph is TRUE
#' 
#' @return Creates a plot showing chromosomal aberrations (no return value)
#' 
#' @details
#' This function visualizes:
#' \itemize{
#'   \item **Gains**: Red rectangles indicating chromosomal gains
#'   \item **Losses**: Semi-transparent blue rectangles indicating losses  
#'   \item **Double aberrations**: Orange rectangles for overlapping gain/loss regions
#'   \item **Chromosome structure**: Gray background showing chromosome boundaries
#'   \item **Sample labels**: Optional y-axis labels (auto-disabled for >50 samples)
#' }
#' 
#' **Visualization features**:
#' \itemize{
#'   \item Proportional chromosome sizing based on actual genomic lengths
#'   \item Automatic sample spacing and positioning
#'   \item Chromosome boundary lines and tick marks
#'   \item Customizable sample labeling
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic plotting with CytoConverter results
#' data <- read.table("cyto_example.txt", sep="\t", header=FALSE)
#' result <- CytoConverter(data)
#' plot_cyto_graph(result$Results)
#' 
#' # Plot with specific reference genome
#' plot_cyto_graph(result$Results, ref_list="hg19")
#' 
#' # Include sample labels for small datasets
#' plot_cyto_graph(result$Results, ylabel=TRUE)
#' 
#' # Use pre-computed cyto_graph output for efficiency
#' graph_data <- cyto_graph(result$Results, "GRCh38")
#' plot_cyto_graph(list_from_cyto=graph_data)
#' 
#' # Include normal samples for comparison
#' plot_cyto_graph(result$Results, include_normals_graph=TRUE, 
#'                 list_of_samples=c("Normal1", "Normal2"))
#' }
#' 
#' @seealso 
#' \code{\link{cyto_graph}} for data preparation
#' \code{\link{cyto_graph_fusion}} for fusion-specific plotting
#' \code{\link{CytoConverter}} for generating input data
#' \itemize{
#'   \item Gains (red regions)
#'   \item Hemizygous losses (blue regions)  
#'   \item Homozygous losses (orange regions)
#' }
#' 
#' Note: This function is optimized for standard gains and losses. For fusion events
#' and structural rearrangements detected with count_fusions=TRUE, use the specialized
#' cyto_graph_fusion() function which provides enhanced visualization for complex
#' chromosomal aberrations.
#' 
#' @examples
#' \dontrun{
#' # Standard plotting for gains/losses
#' result <- CytoConverter(data)
#' plot_cyto_graph(result$Results)
#' 
#' # For fusion data, use cyto_graph_fusion instead:
#' fusion_result <- CytoConverter(data, count_fusions = TRUE)
#' plot_cyto_graph_fusion(fusion_result$Results)
#' }
#' 
#' @seealso 
#' \code{\link{cyto_graph_fusion}} for fusion-specific visualization
#' \code{\link{cyto_graph}} for data preparation
#' 
#' @export
plot_cyto_graph<-function(cyto_list=NULL,list_from_cyto=NULL,ref_list="GRCh38",ylabel=NULL,include_normals_graph=F,list_of_samples=NULL){
  
  if(is.null(ylabel)){
    if(length(uniq_coord_name) < 50){
      ylabel=T
    }else{
      ylabel=F
    }
  }
   if(!is.null(cyto_list))
  {
    list_from_cyto<-cyto_graph(cyto_list,ref_list,include_normals_graph,list_of_samples)
    rect_maker<-list_from_cyto[[1]]
    xbegin<-list_from_cyto[[2]]
    xcoord_master<-list_from_cyto[[3]]
    y_above<-list_from_cyto[[4]]
    y_below<-list_from_cyto[[5]]
    sorted_reflist<-list_from_cyto[[6]]
    cum_length_coords<-list_from_cyto[[7]]
    start_cum_length<-list_from_cyto[[8]]
    uniq_coord_name<-list_from_cyto[[9]]  
  }else if(!is.null(list_from_cyto)){
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
    
    if(ylabel == T)
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
}