
##setting up blank plot
cyto_graph<-function(cyto_list,ref_list="GRCh38"){
  
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
  }
  
  ##cyto_list<-cyto_list[order(cyto_list[,1],order(cyto_list[,2],cyto_list[,3])),]
  
    
  double_loss<-cyto_list[which(cyto_list[,5]=="Loss"),]
  double_loss[,1]<-as.character(double_loss[,1])

  if(nrow(double_loss)>1)
  {
    loss_overlap<-data.frame()
    chr_list<-unique(double_loss[,2])
    name_list<-as.character(unique(double_loss[,1]))
    
  for(i in 1:length(chr_list))
  {
        for(k in 1:length(name_list))
        {
          chr_table<-double_loss[intersect(which(double_loss[,2]==chr_list[i]) , which(double_loss[,1]==name_list[k] )),]
          
          ##don't do this, check for complete overlap first, if so , skip any that are in completely 
          
          ##chr_table_ctr<-cbind("test",unique.data.frame(chr_table[,2:4]),"Loss")
          if(nrow(chr_table) >=2)
          {
            chr_table <- chr_table[order(chr_table[,2],chr_table[,3]),]
            
            for(d in 1:(nrow(chr_table)-1))
            {
              overlap=F
              if(!is.na(chr_table[d,1])){
                for(j in (d+1):(nrow(chr_table))){
                    if(!is.na(chr_table[j,1])){
                    
                      first <- as.numeric(chr_table[d,3:4])
                      sec <- as.numeric(chr_table[j,3:4])
                      if(first %overlaps% sec)
                      {
                        overlap=T
                        ##if(first[2] >  sec[1]){
                         ## chr_table[d,3] <- sec[1]
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
   }
    y_coordlist<-matrix(ncol=2,nrow=0)
    if(is.list(matched_coord_names))
    {
      for(i in 1:length(matched_coord_names))
      {
        y_coordlist<-rbind(y_coordlist,matched_coord_names[[i]])
      }
    }else{
      y_coordlist<-matched_coord_names
    }
  
    xcoord_master=max(nchar(as.character(uniq_coord_name)))*0.005+xbegin
    
    
    ##calculate x coords and y coords according to y_coord_list by resorting cytolist
    if(length(y_coordlist)>0 && nrow(cyto_list)>1)
    {
      cyto_list<-cyto_list[order(cyto_list[,1],y_coordlist[,1]),]
    }
    
    if(is.vector(cyto_list)){
      coords_listed<-cyto_list[2:5]
      temp_coords<-gsub("chr","",coords_listed[1])
      temp_coords[grep("Y",temp_coords)]<-24
      temp_coords[grep("X",temp_coords)]<-23
      coordlist<-as.numeric(temp_coords);
      
      ##adjust chrom name for x and y for input data
      coords_listed[grep("X",coords_listed[,1])]<-23
      coords_listed[grep("Y",coords_listed[,1])]<-24      
      
      y_area_coord=cbind((y_coordlist[,2]-1)*-(y_above-y_below)/length(uniq_coord_name)+y_above,y_above-y_coordlist[,2]*(y_above-y_below)/length(uniq_coord_name))
      
      
      ##calculating where rectangle should start from how long the sample is
      
      
      
      ##x coordinates
      xstart<-(start_cum_length[coordlist]+(as.numeric(coords_listed[2]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
      xend<-(start_cum_length[coordlist]+(as.numeric(coords_listed[3]))/(sum(coords)))*(1-xcoord_master)+xcoord_master
      
      
      ##all info for coordinates
      rect_maker<-as.data.frame(cbind(xstart,xend,y_area_coord,coords_listed[4]))
      rect_maker[1:4]<-apply(rect_maker[,1:4],2,function(x){as.numeric(as.character(x))})
      rect_maker[5]<-as.character(rect_maker[,5])
      
    }else{
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
