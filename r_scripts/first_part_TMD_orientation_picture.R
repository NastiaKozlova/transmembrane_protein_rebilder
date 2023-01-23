#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(magick)
#install.packages("ggraph")
library(ggraph)

library(igraph)
#library(network)
setwd(part_start)
df_start<-read.csv("results/df_first_part.csv",stringsAsFactors = F)
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
df_seq<-data.frame(matrix(nrow = max(c(df_start$first_part_finish,df_start$second_part_finish)),ncol = 2))
colnames(df_seq)<-c("amino","TMD")
df_seq$amino<-c(1:nrow(df_seq))
df_topology<-df_topology%>%filter(type=="Transmembrane")
df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))
for (i in 1:nrow(df_topology)) {
  df_seq$TMD[df_seq$amino>=df_topology$seq_beg[i]&df_seq$amino<=df_topology$seq_end[i]]<-df_topology$type[i]
}
df_seq<-df_seq%>%filter(!is.na(TMD))
part<-paste0(part_start,"results/first_part/")
setwd(part)
i<-1
if (!dir.exists(paste0("make_picture_tcl_center"))){dir.create(paste0("make_picture_tcl_center"))}
i<-1
#make_picture_tcl_center
for (i in 1:nrow(df_start)) {
  
#  df_ring<-read.csv(paste0("interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),stringsAsFactors = F)
#  df_ring<-left_join(df_ring,df_seq,by=c("NodeId1"="amino"))
#  df_ring<-left_join(df_ring,df_seq,by=c("NodeId2"="amino"))
#  colnames(df_ring)<-c("NodeId1", "NodeId2","Interaction", "resid_1", "resid_2", "TMD_1", "TMD_2" )
#  df_ring<-df_ring%>%filter(!is.na(TMD_1))
#  df_ring<-df_ring%>%filter(!is.na(TMD_2))
#  df_ring<-df_ring%>%filter(TMD_1!=TMD_2)
  
  
  #  df_ring<-df_ring%>%mutate(TMD_1=paste0("TM",TMD_1))
  #  df_ring<-df_ring%>%mutate(TMD_2=paste0("TM",TMD_2))
#  df_ring<-df_ring%>%select(TMD_1,TMD_2,bonds_quantity)
  pdb_name<-paste0(part,"structure/",df_start$name[i],"/",df_start$name[i],"_",df_start$group_number[i],".pdb")
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0(#'cd ', part_name,"complex_structure_center/\n\n",
                        'mol new {',pdb_name,'} type {pdb}')
  #  b<-paste('(resid ',df_ring$NodeId1,' and name CA and resname ',df_ring$NodeId2," and name CA)")
  #  a<-paste0(b,collapse = " or ")
    df_tcl[1,2]<-paste0(#'set all [atomselect ',(i-1),' "',a,'"]\n',
    #                    'set i ',(i-1),"\n",
    #                    'foreach atom [$all list] {\n',
    #                    '  label add Atoms ',(i-1),'/$atom\n',
    #                    '  incr i\n}\n',
    #                    '$all delete\n\n',
                        'color Display Background white\n',
                        'color Labels Atoms black\n',
                        'color Labels Bonds black\n\n')
    a<-c(df_start$first_part_start[i] :  df_start$first_part_finish[i])
    a<-a[a%in%df_seq$amino]
    b<-c(df_start$second_part_start[i] :  df_start$second_part_finish[i])
    b<-b[b%in%df_seq$amino]
    df_tcl[1,3]<-paste('mol modselect 0 ',i-1,' all \n',
                        'mol modmaterial 0 ',(i-1),' Transparent\n',
                        'mol modstyle 0 ' ,i-1, ' NewCartoon\n',
                        'mol modcolor 0 ' ,i-1, ' Type \n\n')#,
    
    df_tcl[1,4]<-paste('mol selection  protein and',
                       ' resid ',paste0(a,collapse = " "),
                       ' \n',
                       'mol addrep ',(i-1),'\n',
                       'mol modmaterial 1 ',(i-1),' Opaque\n',
                       'mol modstyle 1 ' ,i-1, ' NewCartoon\n',#,
                       'mol modcolor 1 ' ,i-1, ' ColorID 3 \n')
    
    df_tcl[1,5]<-paste('mol selection  protein and',
                        ' resid ',paste0(b,collapse = " "),
                        ' \n',
                       'mol addrep ',(i-1),'\n',
                        'mol modmaterial 2 ',(i-1),' Opaque\n',
                        'mol modstyle 2 ' ,i-1, ' NewCartoon\n',#,
                          'mol modcolor 2 ' ,i-1, ' ColorID 21 \n')
#                          )#,
#  df_tcl[1,5]<-paste0('mol selection all\n',
#3                      'mol material Transparent\n',
#                      'mol addrep ',(i-1),'\n',
#                      'mol modstyle 2 ',(i-1),' NewCartoon\n',
#                      'mol modcolor 2 ',(i-1),' Type \n')
    
    #  df_tcl[1,4]<-paste0('mol selection resid ',paste0(c(158:161,431:434),collapse = " "))
    #  df_tcl[1,5]<-paste0('mol material Opaque\n',
    #                      'mol addrep ',(i-1),'\n',
    #                      'mol modstyle 1 ',(i-1),' Licorice\n',
    #                      'mol modcolor 1 ',(i-1),' ColorID 16 \n')
    #if (nrow(df_interaction)>0){
    #  df_tcl[1,4]<-paste0('mol selection resname ',paste0(unique(df_interaction$resid.y),collapse = " "))
    #  df_tcl[1,5]<-paste0('mol modmaterial 1 ',(i-1),' Opaque\n',
    #                      'mol addrep ',(i-1),'\n',
    #                      'mol modstyle 1 ',(i-1),' Licorice \n',
    #                      'mol modcolor 1 ',(i-1),' Name\n')
    #  df_tcl[1,6]<-paste0('mol selection (resid ',paste0(unique(df_interaction$resno.x),collapse = " "),")")
    #  df_tcl[1,7]<-paste0(' mol modmaterial 2 ',(i-1),' Opaque\n',
    #                      'mol addrep ',(i-1),'\n',
    #                      'mol modstyle 2 ',(i-1),' Licorice 0.1 12 12 \n',
    #                      'mol modcolor 2 ',(i-1),' Type')
      
    #  for (p in 1:nrow(df_interaction)) {
    #    df_tcl[(p+1),1]<-paste0('set atomID1 [[atomselect ',(i-1),' "(resid ',df_interaction$resno.x[p],
    #                            ' and name ',df_interaction$elety.x[p],
    #                            " and resname ",df_interaction$resid.x[p],')"] list]')
    #    df_tcl[(p+1),2]<-paste0('set atomID2 [[atomselect ',(i-1),
    #                            ' "(x > ',df_interaction$x.y[p]-0.5,' and x < ',df_interaction$x.y[p]+0.5,
    #                            ' and y > ',df_interaction$y.y[p]-0.5,' and y < ',df_interaction$y.y[p]+0.5,
    #                            ' and z > ',df_interaction$z.y[p]-0.5,' and z < ',df_interaction$z.y[p]+0.5,')"] list]')
    #    
    #    df_tcl[(p+1),3]<-paste0('label add Bonds ',(i-1),'/$atomID1 ',(i-1),'/$atomID2')
    #  }
    #}
    df_tcl[is.na(df_tcl)]<-""
    write.csv(df_tcl,paste0("make_picture_tcl_center/",df_start$name[i],"_",df_start$group_number[i],".tcl"),row.names = F)
  }#else{print(i)}
#} 

#v_structure<-list.files("make_picture_tcl_center/")
df_tcl<-read.csv(paste0("make_picture_tcl_center/",df_start$name[1],"_",df_start$group_number[1],".tcl"),stringsAsFactors = F)
i<-2
if(nrow(df_start)>1){
for (i in 2:nrow(df_start)) {
  df_tcl_add<-read.csv(paste0("make_picture_tcl_center/",df_start$name[i],"_",df_start$group_number[i],".tcl"),stringsAsFactors = F)
  df_tcl<-rbind(df_tcl,df_tcl_add)
}
}
write.table(df_tcl,paste0("make_picture_tcl_center.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
