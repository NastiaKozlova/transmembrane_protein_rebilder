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
  setwd(part_fin)
  if(!dir.exists(part_fin)){dir.create(part_fin)}
  if(dir.exists(part_fin)){
    
    setwd(part_fin)
    if (!dir.exists("compare_interaction/")){dir.create("compare_interaction/")}
    #     if (!file.exists("df_interactions_all.csv") | !file.exists("df_RSMD_interactions_all.csv")){
    df_RMSD<-read.csv("df_RMSD_fin.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%mutate(number_interactions=NA)
    df_interactions<-read.csv(paste0("interactions/",df_RMSD$models[1]),stringsAsFactors = F)
    df_interactions<-df_interactions%>%mutate(models=df_RMSD$models[1])
    df_RMSD$number_interactions[1]<-nrow(df_interactions)
    for (i in 2:nrow(df_RMSD)) {
      if(file.exists(paste0("interactions/",df_RMSD$models[i]))){
        df_interactions_add<-read.csv(paste0("interactions/",df_RMSD$models[i]),stringsAsFactors = F)
        df_interactions_add<-df_interactions_add%>%mutate(models=df_RMSD$models[i])
        df_RMSD$number_interactions[i]<-nrow(df_interactions_add)
        df_interactions<-rbind(df_interactions,df_interactions_add)
      }
    }
    df_interactions<-left_join(df_interactions,df_RMSD,by="models")
    write.csv(df_interactions,"df_interactions_all.csv",row.names=F)
    write.csv(df_RMSD,"df_RSMD_interactions_all.csv",row.names=F)
    #     }
    df_RMSD<-read.csv("df_RSMD_interactions_all.csv",stringsAsFactors = F)
    df_interactions<-read.csv("df_interactions_all.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%filter(!is.na(number_interactions))
    print(Sys.time())
    i<-1
    for (i in 1:nrow(df_RMSD)) {
      if(!file.exists(paste0("compare_interaction/",df_RMSD$number[i],".csv"))){
        if(file.exists(paste0("RMSD_test/",df_RMSD$number[i],".csv"))){
          df_RMSD_test<-read.csv(paste0("RMSD_test/",df_RMSD$number[i],".csv"),stringsAsFactors = F)
          df_interactions_a<-df_interactions[df_interactions$number%in%df_RMSD_test$number.x,]
          df_interactions_b<-df_interactions[df_interactions$number%in%df_RMSD_test$number.y,]
          df_interactions_full<-left_join(df_interactions_a,df_interactions_b,by="full_amino_name")
          df_interactions_full<-df_interactions_full%>%filter(!is.na(number.x))
          df_interactions_full<-df_interactions_full%>%filter(!is.na(number.y))
          df_interactions_full<-df_interactions_full%>%group_by(number.y)%>%mutate(number_semi=n())
          df_interactions_full<-ungroup(df_interactions_full)
          df_interactions_full<-df_interactions_full%>%select(number.x, number.y, number_semi,number_interactions.x,number_interactions.y)
          df_interactions_full<-unique(df_interactions_full)
          df_interactions_full<-df_interactions_full%>%mutate(persent.x=number_semi/number_interactions.x*100)
          df_interactions_full<-df_interactions_full%>%mutate(persent.y=number_semi/number_interactions.y*100)
          df_interactions_full<-df_interactions_full%>%filter(persent.x>25)
          df_interactions_full<-df_interactions_full%>%filter(persent.y>25)
          df_interactions_full<-df_interactions_full%>%select(number.x, number.y, number_semi,persent.x,persent.y)
          df_interactions_full<-unique(df_interactions_full)
          write.csv(df_interactions_full,paste0("compare_interaction/",df_RMSD$models[i]),row.names = F)
        }
      }
    }
  }
}
