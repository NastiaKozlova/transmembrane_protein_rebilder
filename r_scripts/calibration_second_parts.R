#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#path_to_PatchDock<-paste0(part_start,"programs/PatchDock/")
#pach_dock_repeats<-1
library(dplyr)
library(bio3d)
library(readr)
library(ggplot2)

v_RMSD<-10

part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
df_start_ad<-read.csv(paste0(part_start,"results/df_first_part.csv"),stringsAsFactors = F)
df_start<-full_join(df_start_ad,df_start,by = c("name", "first_part_model",
                                                "first_part_start", "first_part_finish",
                                                "second_part_model", "second_part_start",
                                                "second_part_finish"))
setwd(part)
i<-1
df_start<-df_start%>%filter(!is.na(third_part_finish))
df_start<-df_start%>%filter(third_part_start>0)
df_start<-df_start%>%filter(third_part_finish>0)
for(w in 1:nrow(df_start)){
  if(!file.exists(paste0(df_start$name[w],"/fin.csv"))){
    df_start$name[w]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(name))
w<-1
if(nrow(df_start)>0){
  df_RMSD<-read.csv(paste0(part,df_start$name[w],"/add_domain/",df_start$group_number[w],"/df_RMSD_all.csv"),stringsAsFactors = F)
  if(nrow(df_start)>1){
    for (w in 2:nrow(df_start)){
      df_RMSD_add<-read.csv(paste0(part,df_start$name[w],"/add_domain/",df_start$group_number[w],"/df_RMSD_all.csv"),stringsAsFactors = F)
      df_RMSD<-rbind(df_RMSD,df_RMSD_add)
    }
  }
}
df_RMSD<-df_RMSD%>%select(RMSD)
df_RMSD<-df_RMSD%>%group_by(RMSD)%>%mutate(count=n())
df_RMSD<-unique(df_RMSD)
p<-ggplot(data=df_RMSD, aes(x=RMSD,y=count))+
  labs(x="RMSD, A")+
  geom_line()+
  scale_x_continuous(breaks = seq(from=0,to=10,by=0.5),labels =  seq(from=0,to=10,by=0.5))+
  theme_bw()
ggsave(p,filename = paste0(part_start,"results/calibrtion_first_parts.png"), width = 20, height = 15, units = c("cm"), dpi = 300 ) 