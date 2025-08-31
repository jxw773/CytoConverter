#' Plot Chromosomal Gains and Losses
#' 
#' @description
#' Creates a graphical visualization of chromosomal gains and losses across samples.
#' The function generates a heatmap-style plot showing gains (red), losses (blue), 
#' and overlapping regions (orange) mapped to chromosomal coordinates.
#' 
#' @param cyto_list Data frame. Results table from CytoConverter containing gains/losses.
#'   If provided, the function will process this data through cyto_graph internally.
#' @param list_from_cyto List. Pre-processed output from cyto_graph function. Use this
#'   if you have already processed the data to avoid redundant computation.
#' @param ref_list Character. Reference genome build for coordinate mapping. Options:
#'   "GRCh38" (default), "hg19", "hg18", "hg17". Can also be a custom reference matrix.
#' @param ylabel Logical or NULL. Whether to display sample names on the y-axis. 
#'   If NULL, automatically determined based on number of samples (<50 shows labels).
#'
#' @return Generates a plot displayed in the current graphics device. No return value.
#'
#' @examples
#' # Basic plotting with CytoConverter results
#' result <- CytoConverter(sample_data)
#' plot_cyto_graph(cyto_list = result$Result)
#' 
#' # Advanced plotting with custom parameters
#' plot_cyto_graph(cyto_list = result$Result, 
#'                ref_list = "hg19",
#'                ylabel = TRUE)
#' 
#' # Plot using pre-processed data
#' graph_data <- cyto_graph(result$Result, "GRCh38")
#' plot_cyto_graph(list_from_cyto = graph_data, ylabel = FALSE)
#'
#' @details
#' The visualization displays:
#' - Red regions: Chromosomal gains
#' - Blue regions: Chromosomal losses  
#' - Orange regions: Overlapping gains and losses (double aberrations)
#' - Gray background: Normal chromosomal regions
#' 
#' Chromosomes are arranged horizontally with samples displayed vertically.
#' The plot automatically scales based on the reference genome build and
#' number of samples provided.
#'
#' @export
#' 
#' @seealso 
#' \code{\link{cyto_graph}} for data preprocessing
#' \code{\link{CytoConverter}} for generating input data
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