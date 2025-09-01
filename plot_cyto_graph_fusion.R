#' Fusion-Specific Cytogenetic Plotting Function with Automatic Color Assignment
#' 
#' @description
#' This function creates visualizations for chromosomal fusions and structural
#' aberrations from CytoConverter analysis with automatic color assignment for
#' different fusion types. It takes the output of cyto_graph_fusion and plots
#' the information using distinct colors for each fusion type.
#' 
#' @param cyto_list Data frame containing CytoConverter fusion results
#' @param list_from_cyto Pre-computed output from cyto_graph_fusion (optional)
#' @param ref_list Reference genome build ("GRCh38", "hg19", "hg18", "hg17") 
#' @param ylabel Boolean flag to enable/disable sample name labels on y-axis (auto-determined if NULL)
#' @param include_normals_graph Boolean flag to include normal samples for comparison (default: FALSE)
#' @param list_of_samples Vector of normal sample names when include_normals_graph is TRUE
#' 
#' @return Creates a plot showing chromosomal fusions with automatic color assignment (no return value)
#' 
#' @details
#' This function visualizes:
#' \itemize{
#'   \item **Fusion types**: Each unique fusion type gets a distinct color
#'   \item **Automatic coloring**: Colors assigned based on fusion type (excluding |chrom portions)
#'   \item **Colorblind-friendly palette**: Uses scientifically validated color schemes
#'   \item **Consistent mapping**: Same fusion type always gets same color
#'   \item **Informative legend**: Shows fusion types and their corresponding colors
#'   \item **Chromosome structure**: Gray background showing chromosome boundaries
#'   \item **Sample labels**: Optional y-axis labels (auto-disabled for >50 samples)
#' }
#' 
#' **Fusion Type Color Assignment**:
#' The function automatically extracts fusion types from tags like:
#' \itemize{
#'   \item "#translocation_balanced|chrom_1" → "translocation_balanced" (blue)
#'   \item "#derivative_chrom::translocation" → "derivative_chrom::translocation" (red)
#'   \item "#insertion_chrom::inserted_piece|chrom_1" → "insertion_chrom::inserted_piece" (green)
#'   \item "#ring_chrom" → "ring_chrom" (yellow)
#' }
#' 
#' **Visualization features**:
#' \itemize{
#'   \item Proportional chromosome sizing based on actual genomic lengths
#'   \item Automatic sample spacing and positioning
#'   \item Chromosome boundary lines and tick marks
#'   \item Customizable sample labeling
#'   \item Color-coded legend for fusion types
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic plotting with CytoConverter fusion results
#' data <- read.table("cyto_fusion_examples.txt", sep="\t", header=FALSE)
#' result <- CytoConverter(data, count_fusions = TRUE)
#' plot_cyto_graph_fusion(result$Results)
#' 
#' # Plot with specific reference genome
#' plot_cyto_graph_fusion(result$Results, ref_list="hg19")
#' 
#' # Include sample labels for small datasets
#' plot_cyto_graph_fusion(result$Results, ylabel=TRUE)
#' 
#' # Use pre-computed cyto_graph_fusion output for efficiency
#' graph_data <- cyto_graph_fusion(result$Results, "GRCh38")
#' plot_cyto_graph_fusion(list_from_cyto=graph_data)
#' 
#' # Include normal samples for comparison
#' plot_cyto_graph_fusion(result$Results, include_normals_graph=TRUE, 
#'                        list_of_samples=c("Normal1", "Normal2"))
#' }
#' 
#' @seealso 
#' \code{\link{cyto_graph_fusion}} for data preparation
#' \code{\link{plot_cyto_graph}} for standard gain/loss plotting
#' \code{\link{CytoConverter}} for generating input data
#' 
#' @note This function is optimized for fusion events and structural rearrangements.
#' For standard gains and losses, use plot_cyto_graph() which provides
#' specialized visualization for chromosomal aberrations without fusion tags.
#' 
#' @export
plot_cyto_graph_fusion<-function(cyto_list=NULL,list_from_cyto=NULL,ref_list="GRCh38",ylabel=NULL,include_normals_graph=F,list_of_samples=NULL){
  
  # Source the fusion colors module
  source("modules/fusion_colors.R")
  
  if(!is.null(cyto_list))
  {
    list_from_cyto<-cyto_graph_fusion(cyto_list,ref_list,include_normals_graph,list_of_samples)
    rect_maker<-list_from_cyto[[1]]
    xbegin<-list_from_cyto[[2]]
    xcoord_master<-list_from_cyto[[3]]
    y_above<-list_from_cyto[[4]]
    y_below<-list_from_cyto[[5]]
    sorted_reflist<-list_from_cyto[[6]]
    cum_length_coords<-list_from_cyto[[7]]
    start_cum_length<-list_from_cyto[[8]]
    uniq_coord_name<-list_from_cyto[[9]]
    fusion_color_mapping<-list_from_cyto[[10]]
    fusion_legend_data<-list_from_cyto[[11]]
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
    fusion_color_mapping<-list_from_cyto[[10]]
    fusion_legend_data<-list_from_cyto[[11]]
  }
  
  # Auto-determine ylabel setting based on number of samples
  if(is.null(ylabel)){
    if(length(uniq_coord_name) < 50){
      ylabel=T
    }else{
      ylabel=F
    }
  }
  
  plot.new()
  
  plot.window(c(0,1),c(0,1),mar=rep(0,4))
  
  rect(xleft=xcoord_master, xright=1,ybottom=y_below,ytop=y_above,col="gray90")
  
  if(nrow(rect_maker)>0)
  {
    apply(rect_maker,1,
          function(x){ 
            # Use the assigned color if available, otherwise default to light gray
            plot_color <- if("Color" %in% names(x) && !is.na(x["Color"])) {
              x["Color"]
            } else {
              "#CCCCCC"  # Light gray for non-fusion or unassigned
            }
            
            rect(xleft=as.numeric(x[1]), xright=as.numeric(x[2]),
                 ybottom=as.numeric(x[3]),ytop=as.numeric(x[4]),
                 col=plot_color,border=NA)
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
  
  # Create legend for fusion types if there are any
  if(length(fusion_color_mapping) > 0) {
    legend(x=1,y= y_below-0.03,legend=fusion_legend_data$display_name,
           fill=fusion_legend_data$color,xjust=1,yjust=1,cex=0.7,title="Fusion Types")
  } else {
    # Fallback legend for non-fusion data
    legend(x=1,y= y_below-0.03,legend=c("No Fusion Data"),fill=c("#CCCCCC"),xjust=1,yjust=1,cex=0.7)
  }
}