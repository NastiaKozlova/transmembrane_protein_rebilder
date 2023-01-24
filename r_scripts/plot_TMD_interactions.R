#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

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
if (!dir.exists(paste0(part_start,"results/first_part/interactions/"))){dir.create(paste0(part_start,"results/first_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/first_part/ring2/"))){dir.create(paste0(part_start,"results/first_part/ring2/"))}
if (!dir.exists(paste0(part_start,"results/first_part/ring/"))){dir.create(paste0(part_start,"results/first_part/ring/"))}
if (!dir.exists(paste0(part_start,"results/first_part/interactions/"))){dir.create(paste0(part_start,"results/first_part/interactions/"))}
if (!dir.exists(paste0(part_start,"results/first_part/TMD_interactions/"))){dir.create(paste0(part_start,"results/first_part/TMD_interactions/"))}
if (!dir.exists(paste0(part_start,"results/first_part/TMD_interactions_plot/"))){dir.create(paste0(part_start,"results/first_part/TMD_interactions_plot/"))}

if (!dir.exists(paste0(part_start,"results/first_part/plot_merge/"))){dir.create(paste0(part_start,"results/first_part/plot_merge/"))}
i<-1


part<-paste0(part_start,"results/first_part/")
setwd(part)
i<-1
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
  
  df_pdb$part[df_pdb$TMD_1>6]<-"second"
  df_pdb$part<-as.factor(df_pdb$part)
  
  g<-graph_from_data_frame(d = df_ring, vertices = df_pdb,
        directed = FALSE)#%%>%
  g <- set_vertex_attr(g, 'color' , value = df_pdb$part)
  g <- set_edge_attr(g, 'bonds' , value = df_ring$bonds_quantity)
  png(paste0("TMD_interactions_plot/",df_start$name[i],"_",df_start$group_number[i],".png"), width = 500, height = 500) 
  plot(g,
       edge.width = E(g)$bonds,
       vertex.label.color = "black",
       layout = layout_nicely(g),
       vertex.color = V(g)$color,
       )
  dev.off()  
  values <- as.numeric(factor(V(g)$color))

  # Calculate the observed assortativity
  observed.assortativity <- assortativity(g, values)
  # Calculate the assortativity of the network randomizing the gender attribute 1000 times
  results<-100000
  v_assortativity <- vector('list', results)
#  v_assortativity.degree <- vector('list', 100000)
  for(j in 1:results){
    v_assortativity[[j]] <- assortativity(g, sample(values),directed = F)
#    v_assortativity.degree[[j]] <- assortativity.degree(g, directed = FALSE)
  }
  v_scale_x<-seq(from=-1,to=1,by=0.1)
  # Plot the distribution of assortativity values and add a red vertical line at the original observed value
  v_assortativity<-(unlist(v_assortativity))
#  v_assortativity.degree<-(unlist(v_assortativity.degree))
  df_results<-data.frame(matrix(nrow=results,ncol=2))
  colnames(df_results)<-c("assortativity","assortativity.degree")
  df_results$assortativity<-v_assortativity
#  df_results$assortativity.degree<-v_assortativity.degree
  p_assortativity<-ggplot(df_results)+
    geom_density(aes(x=assortativity))+#,after_stat(density)),binwidth=0.05)+
    geom_vline(xintercept=observed.assortativity, linetype="dashed", color = "red")+
    theme_bw()+
    scale_x_continuous(breaks=v_scale_x,labels=v_scale_x)

  TMD_interaction<-paste0(part,"TMD_interactions_plot/",df_start$name[i],"_",df_start$group_number[i],".png")
  TMD_pictures<-paste0(part,"pictures/",df_start$name[i],"_",df_start$group_number[i],".png")
  #file.exists(file_name)
  #logo_file <- system.file( file_name, package = "cowplot")
  
  p_TMD_interaction<-ggdraw() +  
    draw_image(TMD_interaction, scale = 1)
  p_TMD_pictures<-ggdraw() +  
    draw_image(TMD_pictures, scale = 1)
  
  p1<-plot_grid(p_TMD_interaction, p_TMD_pictures, labels = c('B', 'C'), label_size = 12,ncol=2)
  merge_plot<-plot_grid(p_assortativity, p1, labels = c('A', ''), label_size = 12,ncol=1)+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  ggsave(merge_plot,filename = paste0(part,"plot_merge/",df_start$name[i],"_",df_start$group_number[i],".png"), width = 20, height = 20, units = c("cm"), dpi = 200 )
}

