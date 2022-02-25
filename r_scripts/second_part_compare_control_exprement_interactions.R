#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)

library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
df_control<-read.csv(paste0("results/control/interactions/start.txt"),stringsAsFactors = F)
df_control<-df_control%>%filter(abs(NodeId1-NodeId2)>5)

df_second_part<-read.csv("results/df_second_part.csv",stringsAsFactors = F)
i<-2
df_second_part<-df_second_part%>%mutate(group_number=paste0("first_part_",first_part_group_number,
                                                            "_second_part_",second_part_group_number))
df_second_part<-df_second_part%>%mutate(file_name=paste0("first_part_",first_part_orientarion,"_",first_part_group_number,
                                                            "_second_part_",second_part_orientarion,"_",second_part_group_number))
df_second_part<-df_second_part%>%mutate(anti_control=NA)
df_second_part<-df_second_part%>%mutate(anti_experiment=NA)
df_second_part<-df_second_part%>%mutate(semi=NA)
i<-1
if(!dir.exists("results/second_part/plot/analysis/")){dir.create("results/second_part/plot/analysis/")}
for (i in 1:nrow(df_second_part)) {
  v_second_part<-c(df_second_part$third_part_start[i]:df_second_part$third_part_finish[i],
                 df_second_part$second_part_start[i]:df_second_part$second_part_finish[i],
                 df_second_part$first_part_start[i]:df_second_part$first_part_finish[i] )
  df_second_part_conrol<-df_control[df_control$NodeId1%in%v_second_part,]
  df_second_part_conrol<-df_second_part_conrol[df_second_part_conrol$NodeId2%in%v_second_part,]
  df_second_part_experiment<-read.csv(paste0("results/second_part/interactions/",df_second_part$name[i],"_",df_second_part$group_number[i],".txt"),stringsAsFactors = F)
  df_second_part_experiment<-df_second_part_experiment%>%filter(abs(NodeId1-NodeId2)>5)
  
  
  df_second_part_anti_control<-anti_join(df_second_part_experiment,df_second_part_conrol)
  df_second_part_anti_experiment<-anti_join(df_second_part_conrol,df_second_part_experiment)
  df_second_part_semi<-semi_join(df_second_part_conrol,df_second_part_experiment)
  
  df_second_part_semi<-df_second_part_semi%>%mutate(type="common")
  df_second_part_anti_control<-df_second_part_anti_control%>%mutate(type="only experiment")
  df_second_part_anti_experiment<-df_second_part_anti_experiment%>%mutate(type="only control")
  
  df_second_part_all<-rbind(df_second_part_semi,df_second_part_anti_control,df_second_part_anti_experiment)
  
  df_second_part$anti_control[i]<-nrow(df_second_part_anti_control)
  df_second_part$anti_experiment[i]<-nrow(df_second_part_anti_experiment)
  df_second_part$semi[i]<-nrow(df_second_part_semi)
  
  df_second_part_all<-df_second_part_all%>%mutate(Domain_1=NA)
  df_second_part_all<-df_second_part_all%>%mutate(Domain_2=NA)
  df_second_part_all$Domain_1[df_second_part_all$NodeId1%in%c(df_second_part$first_part_start[i]:df_second_part$first_part_finish[i])]<-df_second_part$first_part_model[i]
  df_second_part_all$Domain_1[df_second_part_all$NodeId1%in%c(df_second_part$second_part_start[i]:df_second_part$second_part_finish[i])]<-df_second_part$second_part_model[i]
  df_second_part_all$Domain_1[df_second_part_all$NodeId1%in%c(df_second_part$third_part_start[i]:df_second_part$third_part_finish[i])]<-df_second_part$third_part_finish[i]
  
  
  df_second_part_all$Domain_2[df_second_part_all$NodeId2%in%c(df_second_part$first_part_start[i]:df_second_part$first_part_finish[i])]<-df_second_part$first_part_model[i]
  df_second_part_all$Domain_2[df_second_part_all$NodeId2%in%c(df_second_part$second_part_start[i]:df_second_part$second_part_finish[i])]<-df_second_part$second_part_model[i]
  df_second_part_all$Domain_2[df_second_part_all$NodeId2%in%c(df_second_part$third_part_start[i]:df_second_part$third_part_finish[i])]<-df_second_part$third_part_finish[i]
  
  
  df_second_part_all<-df_second_part_all%>%filter(Domain_1!=Domain_2)
  
  p<-ggplot()+
    
    geom_rect(aes(xmin=df_second_part$first_part_start[i],xmax=df_second_part$first_part_finish[i],ymin=-Inf,ymax=Inf,alpha=0.5))+
    geom_rect(aes(xmin=df_second_part$second_part_start[i],xmax=df_second_part$second_part_finish[i],ymin=-Inf,ymax=Inf,alpha=0.5))+
    geom_rect(aes(xmin=df_second_part$third_part_start[i],xmax=df_second_part$third_part_finish[i],ymin=-Inf,ymax=Inf,alpha=0.5))+
    
    geom_rect(aes(ymin=df_second_part$first_part_start[i],ymax=df_second_part$first_part_finish[i],xmin=-Inf,xmax=Inf,alpha=0.5))+
    geom_rect(aes(ymin=df_second_part$second_part_start[i],ymax=df_second_part$second_part_finish[i],xmin=-Inf,xmax=Inf,alpha=0.5))+
    geom_rect(aes(ymin=df_second_part$third_part_start[i],ymax=df_second_part$third_part_finish[i],xmin=-Inf,xmax=Inf,alpha=0.5))+
    geom_point(aes(x=NodeId1,y=NodeId2,colour=type),data=df_second_part_all)+
    theme_bw()+theme(legend.position = "bottom")
  ggsave(p,filename = paste0("results/second_part/plot/analysis/",df_second_part$name[i],"_",df_second_part$file_name[i],".png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
}
