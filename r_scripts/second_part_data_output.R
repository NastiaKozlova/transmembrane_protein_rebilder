#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
library(ggplot2)

setwd(part_start)

i<-2
sort_structures<-function(df_start,i){
  v_start<-length(list.files(paste0(df_start$name[i],"/add_domain/",df_start$first_part_group_number[i],"/structure/")))
  df_start_all<-read.csv(paste0(df_start$name[i],"/add_domain/",df_start$first_part_group_number[i],"/fin.csv"),stringsAsFactors = F)
  df_start_all<-df_start_all%>%mutate(angle=angle*90/pi*2)
  df_start_all<-df_start_all%>%mutate(angle_mem=angle_mem*90/pi*2)
  
  df_start_all<-df_start_all%>%mutate(persent_align=round(align_models/group_models*100,digits = 1))
  df_start_all<-df_start_all%>%mutate(orientarion="between")
  
  df_start_all$orientarion[abs(df_start_all$angle)<60]<-"as WT"
  df_start_all$orientarion[abs(df_start_all$angle)>120]<-"inverted"
  df_start_all$orientarion[df_start_all$orientarion=="as WT"&df_start_all$RMSD<5]<-"WT"
  df_start_all$orientarion[abs(df_start_all$angle_mem)<30]<-"between"
  
  test<-df_start$third_part_finish[i]-df_start$third_part_start[i]
  if(test>50){
    df_start_all<-df_start_all%>%filter(orientarion!="between")
  }
  df_start_all<-df_start_all%>%mutate(name=df_start$name[i])
  df_start_all<-df_start_all%>%select(name,orientarion,RMSD,persent_align,group_models,angle, group_number,bond_energy_fs,align_models)
  colnames(df_start_all)<-c("name","second_part_orientarion","second_part_RMSD","second_part_persent_align", "second_part_group_models", 
                            "second_part_angle", "second_part_group_number",  "second_part_bond_energy_fs","second_part_align_models")
  df_start_all<-df_start_all%>%mutate(first_part_group_number=df_start$first_part_group_number[i])
  df_start_all<-df_start_all%>%filter(second_part_group_models>5)
  df_start_all<-df_start_all%>%filter(second_part_group_models>v_start/1000)
  df_start_all<-df_start_all%>%filter(second_part_group_models>=quantile(df_start_all$second_part_group_models,0.95))

  if(length(df_start_all$second_part_orientarion[df_start_all$second_part_orientarion%in%"WT"])>0){
    v_energy_test<-max(df_start_all$second_part_bond_energy_fs[df_start_all$second_part_orientarion%in%"WT"])
    df_start_all<-df_start_all%>%filter(second_part_bond_energy_fs<=v_energy_test)
  }else{
    df_start_all<-df_start_all%>%filter(second_part_bond_energy_fs<=quantile(df_start_all$second_part_bond_energy_fs,probs = 0.25))
  }
  
  df_start_all<-df_start_all%>%mutate(second_part_frequence=second_part_group_models/sum(df_start_all$second_part_group_models)*100) 
  df_start_all<-df_start_all%>%mutate(second_part_frequence=round(second_part_frequence,digits = 1))
  df_start_all<-df_start_all%>%mutate(second_part_RMSD=round(second_part_RMSD,digits = 2))
  df_start_all<-df_start_all%>%mutate(second_part_angle=round(second_part_angle,digits = 0))
  df_start_all<-df_start_all%>%mutate(second_part_persent_align=round(second_part_align_models/second_part_group_models*100,digits = 1))
  df_start_all$second_part_align_models<-NULL
  df_start_all$second_part_group_models<-NULL
  return(df_start_all)
}
part<-paste0(part_start,"structure_prediction/")

start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
df_start_ad<-read.csv(paste0(part_start,"results/df_first_part.csv"),stringsAsFactors = F)
df_start_ad$group_models<-NULL
df_start<-full_join(df_start_ad,df_start,by = c("name", "first_part_model",
                                                "first_part_start", "first_part_finish",
                                                "second_part_model", "second_part_start",
                                                "second_part_finish"))
colnames(df_start)<-c( "name",
                       "first_part_orientarion", "first_part_RMSD" , "first_part_frequence", "first_part_persent_align" ,
                       "first_part_angle", "first_part_model",  "first_part_start",  "first_part_finish", 
                       "second_part_model", "second_part_start", "second_part_finish","first_part_group_number",      
                       "first_part_bond_energy_fs", "third_part_model",  "third_part_start", "third_part_finish" )

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

#check total quantity of patchdock compexes

df_start_all<-sort_structures(df_start,1)
i<-2
if(nrow(df_start)>1){
  for (i in 2:nrow(df_start)) {
    df_start_add<-sort_structures(df_start,i)
    
    df_start_all<-rbind(df_start_all,df_start_add)
  }
}



df_start_all<-left_join(df_start,df_start_all,by = c("name", "first_part_group_number"))

colnames(df_start_all)
df_start_all<-df_start_all%>%select(name,first_part_orientarion,second_part_orientarion,

                                      first_part_RMSD,second_part_RMSD,
                                      first_part_frequence,first_part_persent_align,first_part_angle,first_part_bond_energy_fs,
                                      second_part_frequence,second_part_persent_align,second_part_angle,second_part_bond_energy_fs,
                                      first_part_group_number,second_part_group_number,
                                      first_part_model,first_part_start,first_part_finish,
                                      second_part_model,second_part_start,second_part_finish,
                                      third_part_model,third_part_start,third_part_finish)

df_start_all<-df_start_all%>%select(name,
                                     first_part_orientarion,first_part_RMSD,first_part_persent_align, first_part_frequence,
                                     first_part_group_number, first_part_bond_energy_fs,first_part_angle,
                                     
                                     first_part_model, first_part_start, first_part_finish,
                                     
                                     second_part_orientarion,second_part_RMSD,second_part_persent_align, second_part_frequence,
                                     second_part_group_number, second_part_bond_energy_fs,second_part_angle,
                                     
                                     second_part_model, second_part_start, second_part_finish,
                                     third_part_model,third_part_start,third_part_finish)

if(dir.exists(paste0(part_start,"results/second_part/"))) {system(command = paste0("rm -r ",part_start,"results/second_part/"),ignore.stdout=T,wait = T)}
if(!dir.exists(paste0(part_start,"results/"))){dir.create(paste0(part_start,"results/"))}
if(!dir.exists(paste0(part_start,"results/second_part/"))){dir.create(paste0(part_start,"results/second_part/"))}
if(!dir.exists(paste0(part_start,"results/second_part/structure"))){dir.create(paste0(part_start,"results/second_part/structure"))}




i<-2
write.csv(df_start_all,paste0(part_start,"results/df_second_part.csv"),row.names = F)
for (i in 1:nrow(df_start_all)) {
  if(!dir.exists(paste0(part_start,"results/second_part/structure/",df_start_all$name[i]))){dir.create(paste0(part_start,"results/second_part/structure/",df_start_all$name[i]))}
  if(file.exists(paste0(part,df_start_all$name[i],"/add_domain/",df_start_all$first_part_group_number[i],
                        "/fin_str/",df_start_all$second_part_group_number[i],".pdb"))){
    pdb<-read.pdb(paste0(part,df_start_all$name[i],"/add_domain/",df_start_all$first_part_group_number[i],
                         "/fin_str/",df_start_all$second_part_group_number[i],".pdb"))
    pdb$atom$chain<-"A"
    write.pdb(pdb,paste0(part_start,"results/second_part/structure/",df_start_all$name[i],"/first_part_",df_start_all$first_part_group_number[i],
                                                                            "_second_part_",df_start_all$second_part_group_number[i],".pdb"))
  }else(print(i))
}
df_start_all<-df_start_all%>%mutate(first_second=paste0(first_part_model,"-",second_part_model))
df_start_all<-df_start_all%>%mutate(orientation=paste0(first_part_orientarion,"-",second_part_orientarion))
df_start_all<-df_start_all%>%mutate(first_part_group_number=as.character(first_part_group_number))
df_start_all<-df_start_all%>%mutate(second_part_group_number=as.character(second_part_group_number))
p<-ggplot(data=df_start_all)+
  #  geom_text(aes(x=first_part_group_number,y=second_part_group_number,label=orientation))+
  geom_text(aes(x=first_part_orientarion,y=second_part_orientarion,label=second_part_frequence))+
  facet_grid(first_second~third_part_model,labeller = labeller(.rows = label_both, .cols = label_both))+theme_bw()
ggsave(p,filename = paste0(part_start,"results/all_parts_rebilder.png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
