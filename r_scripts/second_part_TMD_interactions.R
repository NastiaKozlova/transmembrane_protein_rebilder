#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)
setwd(part_start)

if (!dir.exists(paste0(part_start,"results/second_part/pictures/"))){dir.create(paste0(part_start,"results/second_part/pictures/"))}

df_start<-read.csv("results/df_second_part.csv",stringsAsFactors = F)
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
df_seq<-data.frame(matrix(nrow = max(c(df_start$first_part_finish,df_start$second_part_finish)),ncol = 2))
colnames(df_seq)<-c("amino","TMD")
df_seq$amino<-c(1:nrow(df_seq))
#df_color<-df_color%>%filter(domain=="ETM")
#df_topology_add<-df_topology%>%filter(seq_beg==df_color$start[1])

df_topology_add<-df_topology[df_topology$type%in%c("CTM","ETM"),]

df_topology<-df_topology[df_topology$type%in%c("Transmembrane"),]
df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))


df_topology<-rbind(df_topology,df_topology_add)
for (i in 1:nrow(df_topology)) {
  df_seq$TMD[df_seq$amino>=df_topology$seq_beg[i]&df_seq$amino<=df_topology$seq_end[i]]<-df_topology$type[i]
}
df_seq<-df_seq%>%filter(!is.na(TMD))
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring2/"))){dir.create(paste0(part_start,"results/second_part/ring2/"))}
if (!dir.exists(paste0(part_start,"results/second_part/ring/"))){dir.create(paste0(part_start,"results/second_part/ring/"))}
if (!dir.exists(paste0(part_start,"results/second_part/interactions/"))){dir.create(paste0(part_start,"results/second_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/second_part/TMD_interactions/"))){dir.create(paste0(part_start,"results/second_part/TMD_interactions/"))}
i<-1
df_start<-df_start%>%mutate(group_number=paste0("first_part_",first_part_group_number,"_second_part_",second_part_group_number))


df_start<-df_start%>%mutate(script=NA)
for (i in 1:nrow(df_start)) {
  df_start$script[i]<-paste0("cd ",part_start,"programs/dist/\n", part_start,"programs/dist/bin/Ring -i ",
                             part_start,"results/second_part/structure/",df_start$name[i],"/",df_start$group_number[i],".pdb > ",
                             part_start,"results/second_part/ring2/",df_start$name[i],"_",df_start$group_number[i],".txt\n",
                             "rm ",part_start,"results/second_part/structure/",df_start$name[i],"/",df_start$group_number[i],".pdb_fasta*\n",
                             "rm ",part_start,"results/second_part/structure/",df_start$name[i],"/",df_start$group_number[i],".pdb_modified\n")}
df_script<-df_start%>%select(script)
write.table(df_script,paste0(part_start,"results/second_part/ring2.txt"), row.names = F,quote=F,col.names = F)
system(paste0("chmod +x ",part_start,"results/second_part/ring2.txt"),ignore.stdout=T,wait = T)
system(paste0(part_start,"results/second_part/ring2.txt"),ignore.stdout=T,wait = T)

i<-1
for (i in 1:nrow(df_start)) {
  if (file.exists(paste0("results/second_part/ring2/",df_start$name[i],"_",df_start$group_number[i],".txt"))){
 #   if (!file.exists(paste0("results/second_part/ring/",df_start$name[i],"_",df_start$group_number[i],".txt"))){
      df_ring<-read_delim(paste0("results/second_part/ring2/",df_start$name[i],"_",df_start$group_number[i],".txt"), delim = "\t", skip = 11)
      df_ring<-df_ring[1:(which(df_ring$NodeId1%in%"NodeId")-1),]
      df_ring<-df_ring%>%select(NodeId1,Interaction,NodeId2,Distance,Angle,Energy)
      write.csv(df_ring,paste0("results/second_part/ring/",df_start$name[i],"_",df_start$group_number[i],".txt"),row.names = F)
 #   }
  }
}
i<-1
j<-1
for (i in 1:nrow(df_start)) {
  if (file.exists(paste0("results/second_part/ring/",df_start$name[i],"_",df_start$group_number[i],".txt"))){
    pdb<-read.pdb(paste0("results/second_part/structure/",df_start$name[i],"/",df_start$group_number[i],".pdb"))
    df_pdb<-pdb$atom
    df_pdb<-df_pdb%>%select(resid,resno)
    df_pdb<-unique(df_pdb)
    df_ring<-read.csv(paste0("results/second_part/ring/",df_start$name[i],"_",df_start$group_number[i],".txt"),stringsAsFactors = F) 
    for (j in 1:nrow(df_ring)) {
      df_ring$NodeId1[j]<-strsplit(df_ring$NodeId1[j],split = ":",fixed = T)[[1]][2]
      df_ring$NodeId2[j]<-strsplit(df_ring$NodeId2[j],split = ":",fixed = T)[[1]][2]
      df_ring$Interaction[j]<-strsplit(df_ring$Interaction[j],split = ":",fixed = T)[[1]][1]
    }  
    df_ring<-df_ring%>%select(NodeId1,NodeId2,Interaction)
    df_ring<-unique(df_ring)
    df_ring$NodeId1<-as.numeric(df_ring$NodeId1)
    df_ring$NodeId2<-as.numeric(df_ring$NodeId2)

    df_ring<-left_join(df_ring,df_pdb,by=c("NodeId1"="resno"))
    df_ring<-left_join(df_ring,df_pdb,by=c("NodeId2"="resno"))
    colnames(df_ring)<-c("NodeId1", "NodeId2", "Interaction", "resid_1",     "resid_2")


    write.csv(df_ring,paste0("results/second_part/interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),row.names = F)
    df_ring<-df_ring%>%select(NodeId1,NodeId2)
    df_ring<-unique(df_ring)
    df_ring<-left_join(df_ring,df_seq,by=c("NodeId1"="amino" ))
    df_ring<-left_join(df_ring,df_seq,by=c("NodeId2"="amino" ))
    colnames(df_ring)<-c("NodeId1","NodeId2", "TMD_1", "TMD_2")
    df_ring<-df_ring%>%filter(!is.na(TMD_1))
    df_ring<-df_ring%>%filter(!is.na(TMD_2))
    df_ring<-df_ring%>%mutate(bond=paste(TMD_1,TMD_2,sep = "-"))
    df_ring<-df_ring%>%group_by(bond)%>%mutate(bonds_quantity=n())
    df_ring<-ungroup(df_ring)
    df_ring<-df_ring%>%select(TMD_1,TMD_2,bond,bonds_quantity)
    df_ring<-unique(df_ring)
    write.csv(df_ring,paste0("results/second_part/TMD_interactions/",df_start$name[i],"_",df_start$group_number[i],".txt"),row.names = F)
  }
}
