part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
w<-1
protein<-1
setwd(part)
i<-1
name<-1
for(name in 1:nrow(df_start)){
  df_start_add<-read.csv(paste0(df_start$name[name],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){
    df_start_add<-df_start_add%>%mutate(script=NA)
    part_fin<-paste0(part,df_start$name[name],"/add_domain/",df_start_add$group_number[protein],"/")
    setwd(part_fin)
    if (!dir.exists("compare_interaction/")){dir.create("compare_interaction/")}
    if (!file.exists("df_interactions_all.csv")){
      if (!file.exists("df_RSMD_interactions_all.csv")){
        df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
        df_RMSD<-df_RMSD%>%mutate(number_interactions=NA)
        df_interactions<-read.csv(paste0("interactions/",df_RMSD$models[1]),stringsAsFactors = F)
        df_interactions<-df_interactions%>%mutate(model=df_RMSD$models[1])
        df_RMSD$number_interactions[1]<-nrow(df_interactions)
        for (i in 2:nrow(df_RMSD)) {
          if(file.exists(paste0("interactions/",df_RMSD$models[i]))){
            df_interactions_add<-read.csv(paste0("interactions/",df_RMSD$models[i]),stringsAsFactors = F)
            df_interactions_add<-df_interactions_add%>%mutate(model=df_RMSD$models[i])
            df_RMSD$number_interactions[i]<-nrow(df_interactions_add)
            df_interactions<-rbind(df_interactions,df_interactions_add)
          }
        }
        write.csv(df_interactions,"df_interactions_all.csv",row.names=F)
        write.csv(df_RMSD,"df_RSMD_interactions_all.csv",row.names=F)
      }
    }
    if (!dir.exists("compare_interaction/")){dir.create("compare_interaction/")}
    df_RMSD<-read.csv("df_RSMD_interactions_all.csv",stringsAsFactors = F)
    df_interactions<-read.csv("df_interactions_all.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%filter(!is.na(number_interactions))
    df_RMSD<-df_RMSD%>%select("models")
    df_RMSD<-df_RMSD%>%mutate(merge=NA)
    print(Sys.time())
    for (i in 1:nrow(df_RMSD)) {
      if(!file.exists(paste0("compare_interaction/",df_RMSD$models[i]))){
        df_interactions_a<-df_interactions%>%filter(model==df_RMSD$models[i])
        df_interactions_b<-df_interactions%>%filter(model!=df_RMSD$models[i])
        df_interactions_full<-full_join(df_interactions_a,df_interactions_b,by="full_amino_name")
        df_interactions_full<-df_interactions_full%>%filter(!is.na(model.x))
        df_interactions_full<-df_interactions_full%>%filter(!is.na(model.y))
        df_interactions_full<-df_interactions_full%>%group_by(model.y)%>%mutate(number_semi=n())
        df_interactions_full<-df_interactions_full%>%select(model.x, model.y, number_semi)
        df_interactions_full<-unique(df_interactions_full)
        df_interactions_full<-df_interactions_full%>%mutate(persent=number_semi/nrow(df_interactions_a)*100)
        df_interactions_full<-df_interactions_full%>%filter(persent>10)
        df_interactions_1<-df_interactions_full
        df_interactions_full<-full_join(df_interactions_b,df_interactions_a,by="full_amino_name")
        df_interactions_full<-df_interactions_full%>%filter(!is.na(model.x))
        df_interactions_full<-df_interactions_full%>%filter(!is.na(model.y))
        df_interactions_full<-df_interactions_full%>%group_by(model.x)%>%mutate(number_semi=n())
        df_interactions_full<-df_interactions_full%>%select(model.x, model.y, number_semi)
        df_interactions_full<-unique(df_interactions_full)
        df_interactions_full<-df_interactions_full%>%mutate(persent=number_semi/nrow(df_interactions_a)*100)
        df_interactions_full<-df_interactions_full%>%filter(persent>10)
        df_interactions_2<-df_interactions_full
        df_interactions_write<-rbind(df_interactions_1,df_interactions_2)
        write.csv(df_interactions_write,paste0("compare_interaction/",df_RMSD$models[i]),row.names = F)
      }
    }
  }
}
print(Sys.time())
df_interactions<-read.csv("df_interactions_all.csv",stringsAsFactors = F)
df_interactions<-full_join(df_interactions,df_RMSD,by=c("model"="models"))
print(Sys.time())
i<-1
for (i in 1:nrow(df_RMSD_fin)) {
  if(!file.exists(paste0("compare_interaction_test/",df_RMSD_fin$number[i],".csv"))){
    df_RMSD<-read.csv(paste0("RMSD_test/",df_RMSD_fin$number[i],".csv"))
    df_interactions_a<-df_interactions[df_interactions$number%in%df_RMSD$number.x,]
    df_interactions_b<-df_interactions[df_interactions$number%in%df_RMSD$number.y,]
    df_interactions_full<-full_join(df_interactions_a,df_interactions_b,by="full_amino_name")
    df_interactions_full<-df_interactions_full%>%filter(!is.na(number.x))
    df_interactions_full<-df_interactions_full%>%filter(!is.na(number.y))
    df_interactions_full<-df_interactions_full%>%group_by(number.y)%>%mutate(number_semi=n())
    df_interactions_full<-df_interactions_full%>%select(number.x, number.y, number_semi)
    df_interactions_full<-unique(df_interactions_full)
    df_interactions_full<-df_interactions_full%>%group_by(number.y)%>%mutate(num_interactions.y=n())
    df_interactions_full<-df_interactions_full%>%(num_interactions.x=nrow(df_interactions_a))
    df_interactions_full<-df_interactions_full%>%mutate(persent_first=number_semi/num_interactions.x*100)
    df_interactions_full<-df_interactions_full%>%mutate(persent_second=number_semi/num_interactions.y*100)
    df_interactions_full<-df_interactions_full%>%filter(persent_first>10)
    df_interactions_full<-df_interactions_full%>%filter(persent_second>10)
    write.csv(df_interactions_full,paste0("compare_interaction_test/",df_RMSD_fin$number[i],".csv"),row.names = F)
    
  }
}
}
