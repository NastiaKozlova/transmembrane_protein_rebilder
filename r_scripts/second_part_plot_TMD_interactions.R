#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(magick)
library(igraph)
#library(network)
setwd(part_start)
df_start<-read.csv("results/df_second_part.csv",stringsAsFactors = F)
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
df_seq<-data.frame(matrix(nrow = max(c(df_start$second_part_finish,df_start$second_part_finish)),ncol = 2))
colnames(df_seq)<-c("amino","TMD")
df_seq$amino<-c(1:nrow(df_seq))
df_topology<-df_topology%>%filter(type=="Transmembrane")
df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))
for (i in 1:nrow(df_topology)) {
    df_seq$TMD[df_seq$amino>=df_topology$seq_beg[i]&df_seq$amino<=df_topology$seq_end[i]]<-df_topology$type[i]
}
df_seq<-df_seq%>%filter(!is.na(TMD))
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring2/"))){dir.create(paste0(part_start,"results/second_part/ring2/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring/"))){dir.create(paste0(part_start,"results/second_part/ring/"))}
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/TMD_interactions/"))){dir.create(paste0(part_start,"results/second_part/TMD_interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/TMD_interactions_plot/"))){dir.create(paste0(part_start,"results/second_part/TMD_interactions_plot/"))}

if (!dir.exists(paste0(part_start,"results/second_part/plot_merge/"))){dir.create(paste0(part_start,"results/second_part/plot_merge/"))}
if (!dir.exists(paste0(part_start,"results/second_part/plot_arrange/"))){dir.create(paste0(part_start,"results/second_part/plot_arrange/"))}

i<-1


part<-paste0(part_start,"results/second_part/")
setwd(part)
i<-1
df_start<-df_start%>%mutate(group_number=paste0("first_part_",first_part_group_number,"_second_part_",second_part_group_number))

for (i in 1:nrow(df_start)) {
    
    df_ring<-read.csv(paste0("TMD_interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),stringsAsFactors = F)
    df_ring<-df_ring%>%filter(TMD_1!=TMD_2)
    
    df_ring<-df_ring%>%select(TMD_1,TMD_2,bonds_quantity)
    
    
    df_pdb1<-df_ring%>%select(TMD_1)
    df_pdb2<-df_ring%>%select(TMD_2)
    colnames(df_pdb2)<-colnames(df_pdb1)
    df_pdb<-rbind(df_pdb1,df_pdb2)
    df_pdb<-unique(df_pdb)
    df_pdb<-df_pdb%>%mutate(part="first")
    
    df_pdb$part[as.numeric(df_pdb$TMD_1)>6]<-"second"
    df_pdb$part[as.numeric(df_pdb$TMD_1)==7]<-"TM7"
    df_pdb$part[is.na(as.numeric(df_pdb$TMD_1))]<-"not TM"
    df_pdb$part<-as.factor(df_pdb$part)
    
    g<-graph_from_data_frame(d = df_ring, vertices = df_pdb,
                             directed = FALSE)#%%>%
    g <- set_vertex_attr(g, 'color' , value = df_pdb$part)
    g <- set_edge_attr(g, 'bonds' , value = df_ring$bonds_quantity)
    png(paste0("TMD_interactions_plot/",df_start$name[i],"_",df_start$group_number[i],".png"), width = 1000, height = 1000)
    
    plot.igraph(g,
                edge.width = E(g)$bonds,
                vertex.label.color = "black",
                size=100,
                layout = layout_nicely(g),
                vertex.color = V(g)$color,
    )
    dev.off()  
}
