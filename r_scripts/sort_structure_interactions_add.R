part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
setwd(part)
for(w in 1:nrow(df_start)){
  if(!file.exists(paste0(df_start$name[w],"/fin.csv"))){
    df_start$name[w]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(third_part_model))
reformat_data<-function(df_interactions,df_RMSD){
  df_interactions<-left_join(df_interactions,df_RMSD,by=c("number.x"="models"))
  df_interactions<-left_join(df_interactions,df_RMSD,by=c("model.y"="models"))
  df_interactions<-df_interactions%>%mutate(number_min=number.y)
  df_interactions<-df_interactions%>%mutate(number_max=number.x)
  df_interactions$number_min[df_interactions$number.x<df_interactions$number.y]<-df_interactions$number.x[df_interactions$number.x<df_interactions$number.y]
  df_interactions$number_max[df_interactions$number.x<df_interactions$number.y]<-df_interactions$number.y[df_interactions$number.x<df_interactions$number.y]
  df_interactions<-df_interactions%>%filter(number_min<number_max)
  df_interactions<-df_interactions%>%mutate(compation=paste0(number_min,"_",number_max))
  df_interactions<-df_interactions%>%select(number.x,number.y,compation,persent)
  df_interactions<-df_interactions%>%group_by(compation)%>%mutate(test=n())
  df_interactions<-df_interactions%>%filter(test==2)
  df_interactions<-df_interactions%>%filter(number.x==df_interactions$number.x[1])
  return(df_interactions)  
}  
w<-1
protein<-1
for(w in 1:nrow(df_start)){
  df_start_add<-read.csv(paste0(part,df_start$name[w],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){
    df_start_add<-df_start_add%>%mutate(script=NA)
    part_fin<-paste0(part,df_start$name[w],"/add_domain/",df_start_add$group_number[protein],"/")
    if(dir.exists(part_fin)){
      
      setwd(part_fin)
      if (!dir.exists("compare_TEMP/")){dir.create("compare_TEMP/")}
      df_RMSD<-read.csv("df_RMSD_fin.csv",stringsAsFactors = F)
      df_RMSD<-df_RMSD%>%mutate(number_paired_models=NA)
      v_sort<-list.files("compare_interaction/")
      df_RMSD<-df_RMSD[df_RMSD$models%in%v_sort,]
      for (j in 1:nrow(df_RMSD)) {
#        if(file.exists(paste0("compare_interaction/",df_RMSD$models[j]))){
          df_interactions<-read.csv(paste0("compare_interaction/",df_RMSD$models[j]),stringsAsFactors = F)
          df_interactions<-df_interactions%>%filter(persent.x>50)
          df_interactions<-df_interactions%>%filter(persent.y>50)
          df_RMSD$number_paired_models[j]<-nrow(df_interactions)
          df_interactions<-df_interactions%>%filter(number.x<number.y) 
          write.csv(df_interactions,paste0("compare_TEMP/",df_RMSD$number[j],".csv"),row.names = F)
 #       }
      }
      print(Sys.time())
      write.csv(df_RMSD,paste0("df_RMSD_compare_interactions.csv"),row.names = F)
      
      if (!dir.exists("RMSD_TEMP/")){dir.create("RMSD_TEMP/")}
      i<-1
      df_RMSD<-read.csv(paste0("df_RMSD_compare_interactions.csv"),stringsAsFactors =  F)
      df_RMSD<-df_RMSD%>%mutate(group=round(number,digits=(-3)))
      v_group<-unique(df_RMSD$group)
      group_name<-v_group[1]
      for (group_name in v_group) {
        df_RMSD_TEMP<-df_RMSD%>%filter(group==group_name) 
        df_interactions<-read.csv(paste0("compare_TEMP/",df_RMSD_TEMP$number[1],".csv"),stringsAsFactors = F)
        for (i in 2:nrow(df_RMSD_TEMP)) {
          df_interactions_add<-read.csv(paste0("compare_TEMP/",df_RMSD_TEMP$number[i],".csv"),stringsAsFactors = F)
          df_interactions<-rbind(df_interactions,df_interactions_add)
        }
        write.csv(df_interactions,paste0("RMSD_TEMP/",group_name,".csv"),row.names=F)
      }
      v_name<-list.files("RMSD_TEMP/")
      df_interactions<-read.csv(paste0("RMSD_TEMP/",v_name[1]),stringsAsFactors = F)
      for (i in 2:length(v_name)) {
        df_interactions_add<-read.csv(paste0("RMSD_TEMP/",v_name[i]),stringsAsFactors = F)
        df_interactions<-rbind(df_interactions,df_interactions_add)
      }
      rm(df_interactions_add)
      write.csv(df_interactions,"df_interactions_fin.csv",row.names=F)
      system(command = paste0("rm -r  ",part_fin,"/RMSD_TEMP/"),ignore.stdout=T,wait = T)
      system(command = paste0("rm -r  ",part_fin,"/compare_TEMP/"),ignore.stdout=T,wait = T)
    }
  }
}
