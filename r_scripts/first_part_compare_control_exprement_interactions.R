#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)

library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
df_control<-read.csv(paste0("results/control/interactions/start.txt"),stringsAsFactors = F)
df_control<-df_control%>%filter(abs(NodeId1-NodeId2)>5)

df_first_part<-read.csv("results/df_first_part.csv",stringsAsFactors = F)
i<-2
df_first_part<-df_first_part%>%mutate(file_name=paste0("first_part_",orientarion,"_",group_number))
df_first_part<-df_first_part%>%mutate(anti_control=NA)
df_first_part<-df_first_part%>%mutate(anti_experiment=NA)
df_first_part<-df_first_part%>%mutate(semi=NA)
i<-1
if(!dir.exists("results/first_part/plot/")){dir.create("results/first_part/plot/")}
if(!dir.exists("results/first_part/plot/analysis/")){dir.create("results/first_part/plot/analysis/")}

df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
#df_seq<-data.frame(matrix(nrow = max(c(df_start$first_part_finish,df_start$second_part_finish)),ncol = 2))
#colnames(df_seq)<-c("amino","TMD")
#df_seq$amino<-c(1:nrow(df_seq))
df_topology_add<-df_topology[df_topology$type%in%c("CTM","ETM"),]

df_topology<-df_topology[df_topology$type%in%c("Transmembrane"),]
df_topology<-df_topology%>%mutate(type=c(1:nrow(df_topology)))


df_topology<-rbind(df_topology,df_topology_add)
i<-3
for (i in 1:nrow(df_first_part)) {
  v_first_part<-c(#df_first_part$third_part_start[i]:df_first_part$third_part_finish[i],
                   df_first_part$first_part_start[i]:df_first_part$first_part_finish[i],
                   df_first_part$second_part_start[i]:df_first_part$second_part_finish[i] )
  
  df_first_part_conrol<-df_control[df_control$NodeId1%in%v_first_part,]
  df_first_part_conrol<-df_first_part_conrol[df_first_part_conrol$NodeId2%in%v_first_part,]

#    df_first_part_conrol<-df_control
  df_first_part_experiment<-read.csv(paste0("results/first_part/interactions/",df_first_part$name[i],"_",df_first_part$group_number[i],".txt"),stringsAsFactors = F)
  df_first_part_experiment<-df_first_part_experiment%>%filter(abs(NodeId1-NodeId2)>5)
  
  df_first_part_anti_control<-anti_join(df_first_part_experiment,df_first_part_conrol)
  df_first_part_anti_experiment<-anti_join(df_first_part_conrol,df_first_part_experiment)
  df_first_part_semi<-semi_join(df_first_part_conrol,df_first_part_experiment)
  
  df_first_part_semi<-df_first_part_semi%>%mutate(type="common")
  df_first_part_anti_control<-df_first_part_anti_control%>%mutate(type="only experiment")
  df_first_part_anti_experiment<-df_first_part_anti_experiment%>%mutate(type="only control")
  
  df_first_part_all<-rbind(df_first_part_semi,df_first_part_anti_control,df_first_part_anti_experiment)
  

  
  df_first_part$anti_control[i]<-nrow(df_first_part_anti_control)
  df_first_part$anti_experiment[i]<-nrow(df_first_part_anti_experiment)
  df_first_part$semi[i]<-nrow(df_first_part_semi) 
  df_first_part_all<-df_first_part_all
  df_first_part_all<-df_first_part_all%>%mutate(Domain_1=NA)
  df_first_part_all<-df_first_part_all%>%mutate(Domain_2=NA)
  df_first_part_all$Domain_1[df_first_part_all$NodeId1%in%c(df_first_part$first_part_start[i]:df_first_part$first_part_finish[i])]<-df_first_part$first_part_model[i]
  df_first_part_all$Domain_1[df_first_part_all$NodeId1%in%c(df_first_part$second_part_start[i]:df_first_part$second_part_finish[i])]<-df_first_part$second_part_model[i]
  df_first_part_all$Domain_2[df_first_part_all$NodeId2%in%c(df_first_part$first_part_start[i]:df_first_part$first_part_finish[i])]<-df_first_part$first_part_model[i]
  df_first_part_all$Domain_2[df_first_part_all$NodeId2%in%c(df_first_part$second_part_start[i]:df_first_part$second_part_finish[i])]<-df_first_part$second_part_model[i]
  df_first_part_all<-df_first_part_all%>%filter(Domain_1!=Domain_2)
  df_topology_graph<-df_topology%>%filter(seq_beg<=(max(v_first_part)))
  df_topology_graph<-df_topology_graph%>%filter(seq_end>=(min(v_first_part)))
  a<-round((max(v_first_part)-min(v_first_part))/10,-1)
  v_breaks<-seq(from=round(min(v_first_part)-1),to=round(max(v_first_part)-1),by=a)
  p<-ggplot()+
      labs(x="aminoacid number",y="aminoacid number",title =paste0(df_first_part$name[i],"_",df_first_part$file_name[i]) )+
      geom_rect(aes(xmin=seq_beg-0.5,xmax=seq_end+0.5,ymin=-Inf,ymax=Inf,alpha=0.5),data = df_topology_graph)+
      geom_rect(aes(ymin=seq_beg-0.5,ymax=seq_end+0.5,xmin=-Inf,xmax=Inf,alpha=0.5),data = df_topology_graph)+
      
      geom_text(aes(x=(seq_beg+seq_end)/2,y=min(v_first_part),label=type),data = df_topology_graph)+
      geom_text(aes(y=(seq_beg+seq_end)/2,x=min(v_first_part),label=type),data = df_topology_graph)+

    
#    geom_rect(aes(xmin=df_first_part$first_part_start[i],xmax=df_first_part$first_part_finish[i],ymin=-Inf,ymax=Inf,alpha=0.5))+
#    geom_rect(aes(xmin=df_first_part$second_part_start[i],xmax=df_first_part$second_part_finish[i],ymin=-Inf,ymax=Inf,alpha=0.5))+

    
#    geom_rect(aes(ymin=df_first_part$first_part_start[i],ymax=df_first_part$first_part_finish[i],xmin=-Inf,xmax=Inf,alpha=0.5))+
#    geom_rect(aes(ymin=df_first_part$second_part_start[i],ymax=df_first_part$second_part_finish[i],xmin=-Inf,xmax=Inf,alpha=0.5))+

      scale_x_continuous(breaks = v_breaks,labels = v_breaks,limits = c(min(v_first_part)-1,max(v_first_part)+1))+
      scale_y_continuous(breaks = v_breaks,labels = v_breaks,limits = c(min(v_first_part)-1,max(v_first_part)+1))+
    geom_point(aes(x=NodeId1,y=NodeId2,colour=type),data=df_first_part_all)+
    theme_bw()+theme(legend.position = "bottom")
  ggsave(p,filename = paste0("results/first_part/plot/analysis/",df_first_part$name[i],"_",df_first_part$file_name[i],".png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
}
