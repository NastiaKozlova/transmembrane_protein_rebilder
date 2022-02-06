#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)


library(GGally)
library(network)
setwd(part_start)
df_start<-read.csv("results/df_second_part.csv",stringsAsFactors = F)
df_start<-df_start%>%mutate(group_number=paste0("first_part_",first_part_group_number,"_second_part_",second_part_group_number))
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
df_seq<-data.frame(matrix(nrow = max(c(df_start$first_part_finish,df_start$second_part_finish,df_start$third_part_finish)),ncol = 2))
colnames(df_seq)<-c("amino","TMD")
df_seq$amino<-c(1:nrow(df_seq))
df_topology<-df_topology%>%filter(type=="Transmembrane")

df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))
i<-1
for (i in 1:nrow(df_topology)) {
  df_seq$TMD[df_seq$amino>=df_topology$seq_beg[i]&df_seq$amino<=df_topology$seq_end[i]]<-df_topology$type[i]
}
df_seq<-df_seq%>%filter(!is.na(TMD))
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring2/"))){dir.create(paste0(part_start,"results/second_part/ring2/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring/"))){dir.create(paste0(part_start,"results/second_part/ring/"))}
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/TMD_interactions/"))){dir.create(paste0(part_start,"results/second_part/TMD_interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/plot/"))){dir.create(paste0(part_start,"results/second_part/plot/"))}
if (!dir.exists(paste0(part_start,"results/second_part/plot/TMD_interactions"))){dir.create(paste0(part_start,"results/second_part/plot/TMD_interactions"))}
if (!dir.exists(paste0(part_start,"results/second_part/plot/interactions"))){dir.create(paste0(part_start,"results/second_part/plot/interactions"))}
i<-1


part<-paste0(part_start,"results/second_part/")
setwd(part)
i<-3
for (i in 1:nrow(df_start)) {
  
  df_ring<-read.csv(paste0("TMD_interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),stringsAsFactors = F)
  df_ring<-df_ring%>%mutate(TMD_1=paste0("TM",TMD_1))
  df_ring<-df_ring%>%mutate(TMD_2=paste0("TM",TMD_2))
  df_ring<-df_ring%>%filter(TMD_1!=TMD_2)
  df_pdb1<-df_ring%>%select(TMD_1, TMD_2, bond,           bonds_quantity)
  df_pdb2<-df_ring%>%select(TMD_2, TMD_1, bond,           bonds_quantity)
  colnames(df_pdb2)<-colnames(df_pdb1)
  df_pdb<-rbind(df_pdb1,df_pdb2)
  df_pdb<-df_pdb%>%select(TMD_1)
  df_pdb<-unique(df_pdb)
#  df_pdb<-df_pdb%>%mutate(type=paste0("TM",TMD_1))
  df_pdb<-unique(df_pdb)
  if(nrow(df_ring)>1){
    vlist<-list(edges =df_ring,vertices =df_pdb)
    
    rownames(vlist$vertices) <- vlist$vertices$TMD_1
    mm.net <- network(vlist$edges[, 1:4], directed = FALSE)
    mm.net %v% "type" <- as.character(
      vlist$vertices[ network.vertex.names(mm.net), "type"]
    )
    # type color palette
    mm.col <- c("EMD1"= "#ff69b4",
                "TM1"= "#4783B7",
                "EMD2" = "#ff69b4", 
                "loop1"= "#80B682",
                "EMD3" = "#ff69b4",
                "TM2" = "#4783B7",
                "EMD4"  = "#ff69b4",
                "TM3"  = "#4783B7",
                "EMD5" = "#ff69b4",
                "TM4" = "#4783B7",
                "EMD6" = "#ff69b4",
                "loop2" = "#80B682",
                "EMD7"  = "#ff69b4",
                "TM5" = "#4783B7",
                "EMD8"= "#ff69b4",
                "TM6" = "#4783B7",
                "EMD9" = "#ff69b4",
                "TM7"  = "#4783B7",
                "TM8"  = "#4783B7",
                "EMD10"= "#ff69b4",
                "EMD11"= "#ff69b4",
                "EMD12"= "#ff69b4")
    # create plot for ggnet2
    set.seed(10052006)
    p<-ggnet2(mm.net, #edge.color = mm.col[ mm.net %v% "type" ],
              label = T,
              #  node.label = "type",
              # labelon = TRUE, 
              #     label.color = mm.col[ mm.net %v% "type" ],
              edge.size = "bonds_quantity",
              size = 5, vjust = -0.6, mode = "kamadakawai", label.size = 5)
    #ggsave(p,filename=paste0("TMD_interactions.png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
    ggsave(p,filename = paste0("plot/TMD_interactions/",df_start$name[i],"_",df_start$group_number[i],".png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
  }
}
i<-1
for (i in 1:nrow(df_start)) {
  
  df_ring<-read.csv(paste0("interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),stringsAsFactors = F)
  df_ring<-df_ring%>%filter(abs(NodeId1-NodeId2)>5)
  df_ring<-df_ring%>%mutate(Node_1=paste(resid_1,NodeId1))
  df_ring<-df_ring%>%mutate(Node_2=paste(resid_2,NodeId2))
  df_ring<-df_ring%>%select(Node_1,Node_2)
  df_ring<-unique(df_ring)
  pdb<-read.pdb(paste0("structure/",df_start$name[i],"/",df_start$group_number[i],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%select(resid,resno)
  df_pdb<-unique(df_pdb)
  df_pdb<-df_pdb%>%mutate(type=paste(resid,resno))
  df_pdb<-df_pdb%>%select(type)
  if(nrow(df_ring)>1){
    vlist<-list(edges =df_ring,vertices =df_pdb)
    
    rownames(vlist$vertices) <- vlist$vertices$type
    mm.net <- network(vlist$edges[, 1:2], directed = FALSE)
    mm.net %v% "type" <- as.character(
      vlist$vertices[ network.vertex.names(mm.net), "type"]
    )
    # type color palette
    mm.col <- c("EMD1"= "#ff69b4",
                "TM1"= "#4783B7",
                "EMD2" = "#ff69b4", 
                "loop1"= "#80B682",
                "EMD3" = "#ff69b4",
                "TM2" = "#4783B7",
                "EMD4"  = "#ff69b4",
                "TM3"  = "#4783B7",
                "EMD5" = "#ff69b4",
                "TM4" = "#4783B7",
                "EMD6" = "#ff69b4",
                "loop2" = "#80B682",
                "EMD7"  = "#ff69b4",
                "TM5" = "#4783B7",
                "EMD8"= "#ff69b4",
                "TM6" = "#4783B7",
                "EMD9" = "#ff69b4",
                "TM7"  = "#4783B7",
                "TM8"  = "#4783B7",
                "EMD10"= "#ff69b4",
                "EMD11"= "#ff69b4",
                "EMD12"= "#ff69b4")
    # create plot for ggnet2
    set.seed(10052006)
    p<-ggnet2(mm.net, #edge.color = mm.col[ mm.net %v% "type" ],
              label = T,
              #  node.label = "type",
              # labelon = TRUE, 
              #     label.color = mm.col[ mm.net %v% "type" ],
              #edge.size = "bonds_quantity",
              size = 5, vjust = -0.6, mode = "kamadakawai", label.size = 5)
    #ggsave(p,filename=paste0("TMD_interactions.png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
    ggsave(p,filename = paste0("plot/interactions/",df_start$name[i],"_",df_start$group_number[i],".png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
  }
}
