#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(stringr)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
#library(magick)
#install.packages("ggraph")
library(ggraph)

library(igraph)
#library(network)
setwd(part_start)
df_start<-read.csv("results/df_second_part.csv",stringsAsFactors = F)
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
df_color<-read.csv("start/vmd_coloring.csv",stringsAsFactors = F)
df_seq<-data.frame(matrix(nrow = max(c(df_start$first_part_finish,df_start$second_part_finish)),ncol = 2))
colnames(df_seq)<-c("amino","TMD")
df_seq$amino<-c(1:nrow(df_seq))
df_topology_add<-df_topology[df_topology$type%in%c("CTM","ETM"),]

df_topology<-df_topology[df_topology$type%in%c("Transmembrane"),]
df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))


df_topology<-rbind(df_topology,df_topology_add)
for (i in 1:nrow(df_topology)) {
    df_seq$TMD[df_seq$amino>=df_topology$seq_beg[i]&df_seq$amino<=df_topology$seq_end[i]]<-df_topology$type[i]
}
df_seq<-df_seq%>%filter(!is.na(TMD))
part<-paste0(part_start,"results/second_part/")
setwd(part)
i<-1
if (!dir.exists(paste0("make_picture_tcl_center"))){dir.create(paste0("make_picture_tcl_center"))}
i<-1
#make_picture_tcl_center
for (i in 1:nrow(df_start)) {
    
    pdb_name<-paste0(part,"structure/",df_start$name[i],"/first_part_",df_start$first_part_group_number[i],
                     "_second_part_",df_start$second_part_group_number[i],".pdb")
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 2+nrow(df_color)))
    df_tcl[1,1]<-paste0(#'cd ', part_name,"complex_structure_center/\n\n",
        'mol new {',pdb_name,'} type {pdb}\n')
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
    #    a<-c(df_start$first_part_start[i] :  df_start$first_part_finish[i])
    #    a<-a[a%in%df_seq$amino]
    #    b<-c(df_start$second_part_start[i] :  df_start$second_part_finish[i])
    #    b<-b[b%in%df_seq$amino]
    df_tcl[1,3]<-paste('mol modselect 0 ',i-1,' all \n',
                       'mol modmaterial 0 ',(i-1),' Transparent\n',
                       'mol modstyle 0 ' ,i-1, ' NewCartoon\n',
                       'mol modcolor 0 ' ,i-1, ' Type \n\n')#,
    toster<-str_detect(pattern=df_color$domain,df_start$name[i])
    df_colored<-df_color[toster,]
    for (p in 1:nrow(df_colored)) {
        
        a<-c(df_colored$start[p] :  df_colored$finish[p])
        a<-a[a%in%df_seq$amino]
        df_tcl[1,p+3]<-paste('mol selection  protein and',
                             ' resid ',paste0(a,collapse = " "),
                             ' \n',
                             'mol addrep ',(i-1),'\n',
                             'mol modmaterial ',p, ' ',(i-1),' Opaque\n',
                             'mol modstyle ',p, ' ',i-1, ' NewCartoon\n',#,
                             'mol modcolor ',p, ' ',i-1, ' ColorID ',df_colored$colour[p],' \n')
        df_tcl[2,p+3]<-paste('mol selection  protein and',
                             ' resid ',paste0(df_colored$selected_amino[p],collapse = " "),
                             'and type CA \n',
                             'mol addrep ',(i-1),'\n',
                             'mol modmaterial ',p+nrow(df_colored), ' ',(i-1),' Opaque\n',
                             'mol modstyle ',p+nrow(df_colored), ' ',i-1, ' Surf\n',#,
                             'mol modcolor ',p+nrow(df_colored), ' ',i-1, ' ColorID ',df_colored$colour[p],' \n')
        
    }
    write.csv(df_tcl,paste0("make_picture_tcl_center/",df_start$name[i],"_first_part_",df_start$first_part_group_number[i],
                            "_second_part_",df_start$second_part_group_number[i],".tcl"),row.names = F)
}

df_tcl<-read.csv(paste0("make_picture_tcl_center/",df_start$name[1],"_first_part_",df_start$first_part_group_number[1],
                        "_second_part_",df_start$second_part_group_number[1],".tcl"),stringsAsFactors = F)
i<-2
if(nrow(df_start)>1){
    for (i in 2:nrow(df_start)) {
        df_tcl_add<-read.csv(paste0("make_picture_tcl_center/",df_start$name[i],"_first_part_",df_start$first_part_group_number[i],
                                    "_second_part_",df_start$second_part_group_number[i],".tcl"),stringsAsFactors = F)
        df_tcl<-rbind(df_tcl,df_tcl_add)
    }
}
write.table(df_tcl,paste0("make_picture_tcl_center.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
