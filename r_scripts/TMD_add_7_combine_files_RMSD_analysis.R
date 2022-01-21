#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#path_to_PatchDock<-paste0(part_start,"programs/PatchDock/")
#pach_dock_repeats<-1

v_RMSD<-10

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
for (w in 1:nrow(df_start)) {
  part_fin<-paste0(part,df_start$name[w],"/add_domain/",df_start$group_number[w],"/")
  if(!dir.exists(part_fin)){dir.create(part_fin)}
  setwd(part_fin)
  if (!dir.exists("RMSD/")){dir.create("RMSD/")}
  v_RMSD_files<-list.files(paste0("RMSD/"))
  if(length(v_RMSD_files)>0){
    df_RMSD<-read.csv(paste0("RMSD/",v_RMSD_files[1]),stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%filter(RMSD<v_RMSD)
    for (i in 2:length(v_RMSD_files)) {
      df_RMSD_add<-read.csv(paste0("RMSD/",v_RMSD_files[i]),stringsAsFactors = F)
      df_RMSD_add<-df_RMSD_add%>%filter(RMSD<v_RMSD)
      df_RMSD<-rbind(df_RMSD,df_RMSD_add)
    }
    write.csv(df_RMSD,"df_RMSD_all.csv",row.names = F)
    rm(df_RMSD)
  }
}
