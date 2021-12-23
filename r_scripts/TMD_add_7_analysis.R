#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#path_to_PatchDock<-paste0(part_start,"programs/PatchDock/")
#pach_dock_repeats<-1
library(dplyr)
library(bio3d)
library(readr)
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
for(w in 1:nrow(df_start)){
  part_fin<-paste0(part,df_start$name[w],"/add_domain/",df_start$group_number[w],"/")
  if(!dir.exists(part_fin)){dir.create(part_fin)}
  setwd(part_fin)
  if (!dir.exists("RMSD/")){dir.create("RMSD/")}
  df_RMSD<-read.csv("df_interactions_fin.csv",stringsAsFactors = F)
  df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
  df_RMSD<-df_RMSD%>%select(number.x,number.y,RMSD)
  v_model<-unique(df_RMSD$number.x)
  b<-c()
  for (j in 1:length(v_model)) {
    if (!file.exists(paste0("RMSD/",v_model[j],".csv")) ){
      b<-c(b,v_model[j])
    }
  }
  v_model<-b
  v_model<-sort(v_model)
  df_RMSD_TEMP<-df_RMSD
  df_RMSD<-df_RMSD[df_RMSD$number.x%in%v_model,]
  if(length(v_model)>0){
    for (j in 1:length(v_model)) {
      df_RMSD_TEMP<-  df_RMSD%>%filter(number.x==v_model[j])
      pdb_1<-read.pdb(paste0("structure/out.txt.",v_model[j],".pdb"))
      pdb_1.int<-atom.select(pdb_1, "calpha")
      for (i in 1:nrow(df_RMSD_TEMP)) {
        pdb_2<-read.pdb(paste0("structure/out.txt.",df_RMSD_TEMP$number.y[i],".pdb"))
        pdb_2.int<-atom.select(pdb_2, "calpha")
        df_RMSD_TEMP$RMSD[i]<-rmsd(a=pdb_1,b=pdb_2,a.inds = pdb_1.int,b.inds = pdb_2.int,fit=T)
  
      }
      df_RMSD_TEMP<-df_RMSD_TEMP%>%filter(RMSD<10)
      write.csv(df_RMSD_TEMP,paste0("RMSD/",v_model[j],".csv"),row.names = F)
      df_RMSD<-  df_RMSD%>%filter(number.x!=v_model[j])
    }
  }
}
