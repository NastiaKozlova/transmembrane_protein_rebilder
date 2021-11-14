#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
setwd(part_start)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
#w<-1
#protein<-1
setwd(part)
for(w in 1:nrow(df_start)){
  if(!file.exists(paste0(df_start$name[w],"/fin.csv"))){
    df_start$name[w]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(name))
df_start_all<-read.csv(paste0(df_start$name[1],"/fin.csv"),stringsAsFactors = F)
v_start<-length(list.files(paste0(df_start$name[1],"/structure")))
df_start_all<-df_start_all%>%mutate(name=df_start$name[1])
df_start_all<-left_join(df_start_all,df_start)
df_start_all<-df_start_all%>%filter(group_models>5)
df_start_all<-df_start_all%>%filter(group_models>v_start/1000)
df_start_all<-df_start_all%>%filter(group_models>=quantile(df_start_all$group_models,0.975))

df_start_all<-df_start_all%>%mutate(frequence=group_models/sum(df_start_all$group_models)*100)
if(nrow(df_start)>1){
  for (i in 2:nrow(df_start)) {
    df_start_add<-read.csv(paste0(df_start$name[i],"/fin.csv"),stringsAsFactors = F)
    v_start<-length(list.files(paste0(df_start$name[i],"/structure")))
    df_start_add<-df_start_add%>%mutate(name=df_start$name[i])
    df_start_add<-left_join(df_start_add,df_start)
    df_start_add<-df_start_add%>%filter(group_models>5)
    df_start_add<-df_start_add%>%filter(group_models>v_start/1000)
    df_start_add<-df_start_add%>%filter(group_models>=quantile(df_start_add$group_models,0.975))
    df_start_add<-df_start_add%>%mutate(frequence=group_models/sum(df_start_add$group_models)*100)
    df_start_all<-rbind(df_start_all,df_start_add)
  }
}
df_groups<-df_start_all
#df_groups<-df_groups%>%mutate(angle=angle*90/pi*2)
df_groups<-df_groups%>%mutate(angle_mem=angle_mem*90/pi*2)
df_groups<-df_groups%>%mutate(orientarion="between")
df_groups$orientarion[abs(df_groups$angle)<60]<-"as WT"
df_groups$orientarion[abs(df_groups$angle)>120]<-"inverted"
df_groups$orientarion[df_groups$orientarion=="as WT"&df_groups$RMSD<5]<-"WT"
df_groups<-df_groups%>%select(name,group_number,group_models,align_models,name,      
                              RMSD,bond_energy,bond_energy_fs,angle,angle_mem,
                              orientarion,min_RMSD,max_RMSD,frequence,
                              first_part_model,first_part_start,first_part_finish,
                              second_part_model,second_part_start, 
                              second_part_finish,group_number,bond_energy,bond_energy_fs)

df_start_all<-df_groups
df_start_all<-df_start_all%>%mutate(RMSD=round(RMSD,digits = 2))
df_start_all<-df_start_all%>%mutate(frequence=round(frequence,digits = 1))
df_start_all<-df_start_all%>%mutate(angle=round(angle,digits = 0))
df_start_all<-df_start_all%>%mutate(persent_align=round(align_models/group_models*100,digits = 1))
df_start_all<-df_start_all%>%mutate(orientarion="between")
df_start_all$orientarion[abs(df_start_all$angle)<60]<-"as WT"
df_start_all$orientarion[abs(df_start_all$angle)>120]<-"inverted"
df_start_all$orientarion[df_start_all$orientarion=="as WT"&df_start_all$RMSD<10]<-"WT"
#df_start_all<-df_start_all%>%select(name,group_number,group_models,align_models,name,      
#                                    RMSD,bond_energy,bond_energy_fs,angle,angle_mem,
#                                    orientarion,min_RMSD,max_RMSD)
df_start_all<-df_start_all%>%select(name,orientarion,RMSD,frequence, RMSD,persent_align,group_models, angle,
                                    first_part_model,first_part_start,first_part_finish,second_part_model,second_part_start, 
                                    second_part_finish,group_number,bond_energy,bond_energy_fs)
if(!dir.exists(paste0(part_start,"results/"))){dir.create(paste0(part_start,"results/"))}
if(!dir.exists(paste0(part_start,"results/first_part/"))){dir.create(paste0(part_start,"results/first_part/"))}
if(!dir.exists(paste0(part_start,"results/first_part/structure"))){dir.create(paste0(part_start,"results/first_part/structure"))}
i<-1
write.csv(df_start_all,paste0(part_start,"results/df_first_part.csv"),row.names = F)
for (i in 1:nrow(df_start_all)) {
  if(!dir.exists(paste0(part_start,"results/first_part/structure/",df_start_all$name[i]))){dir.create(paste0(part_start,"results/first_part/structure/",df_start_all$name[i]))}
  if(file.exists(paste0(part,df_start_all$name[i],"/fin_str/",df_start_all$group_number[i],".pdb"))){
    pdb<-read.pdb(paste0(part,df_start_all$name[i],"/fin_str/",df_start_all$group_number[i],".pdb"))
    write.pdb(pdb,paste0(part_start,"results/first_part/structure/",df_start_all$name[i],"/",df_start_all$name[i],"_",df_start_all$group_number[i],".pdb"))
  }else(print(i))
}
