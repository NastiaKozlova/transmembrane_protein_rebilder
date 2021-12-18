#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
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
  if(dir.exists(part_fin)){
    
    setwd(part_fin)
    if (!dir.exists("RMSD_test/")){dir.create("RMSD_test/")}
    df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%mutate(number=NA)
    for (i in 1:nrow(df_RMSD)) {
      df_RMSD$number[i]<-strsplit(df_RMSD$models[i],split = ".",fixed = T)[[1]][3]
    }
    df_RMSD<-df_RMSD%>%select(models,number,RMSD)
    df_RMSD<-df_RMSD%>%mutate(number=as.numeric(number))
    df_RMSD_fin<-df_RMSD%>%mutate(rigth_models=0)
    df_RMSD<-df_RMSD%>%mutate(c="C")
    for (i in 1:nrow(df_RMSD)) {
      df_RMSD_add<-df_RMSD%>%filter(number==df_RMSD$number[i])
      df_RMSD_test<-full_join(df_RMSD_add,df_RMSD,by="c")
      df_RMSD_test<-df_RMSD_test%>%filter(abs(RMSD.x-RMSD.y)<5)
      df_RMSD_fin$rigth_models[i]<-nrow(df_RMSD_test)
      df_RMSD_test<-df_RMSD_test%>%select(number.x,number.y)
      df_RMSD_test<-df_RMSD_test%>%filter(number.x!=number.y)
      write.csv(df_RMSD_test,paste0("RMSD_test/",df_RMSD$number[i],".csv"),row.names = F)
    }
    write.csv(df_RMSD_fin,"df_RMSD_fin.csv",row.names = F)
  }
}
