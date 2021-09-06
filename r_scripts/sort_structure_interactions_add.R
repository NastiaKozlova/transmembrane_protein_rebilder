part_start = commandArgs(trailingOnly=TRUE)
Sys.time()
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
reformat_data<-function(df_interactions,df_RMSD){
  df_interactions<-left_join(df_interactions,df_RMSD,by=c("model.x"="models"))
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
for(name in 1:nrow(df_start)){
  df_start_add<-read.csv(paste0(df_start$name[name],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){
    df_start_add<-df_start_add%>%mutate(script=NA)
    part_fin<-paste0(part,df_start$name[name],"/add_domain/",df_start_add$group_number[protein],"/")
    setwd(part_fin)
    if (!dir.exists("compare_TEMP/")){dir.create("compare_TEMP/")}
    df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%mutate(number=NA)
    for (i in 1:nrow(df_RMSD)) {
      df_RMSD$number[i]<-strsplit(df_RMSD$models[i],split = ".",fixed = T)[[1]][3]
    }
    df_RMSD<-df_RMSD%>%select(models,number)
    df_RMSD<-df_RMSD%>%mutate(number=as.numeric(number))
    df_RMSD<-df_RMSD%>%mutate(number_models=NA)
    print(Sys.time())
    for (j in 1:nrow(df_RMSD)) {
      df_interactions<-read.csv(paste0("compare_interaction/",df_RMSD$models[j]),stringsAsFactors = F)
      df_interactions<-df_interactions%>%filter(persent>50)
      df_interactions<-reformat_data(df_interactions,df_RMSD)
      df_RMSD$number_models[j]<-nrow(df_interactions)
      write.csv(df_interactions,paste0("compare_TEMP/",df_RMSD$number[j],".csv"),row.names = F)
    }
    write.csv(df_interactions,"df_interactions_all.csv",row.names=F)
    df_interactions<-read.csv("df_interactions_all.csv",stringsAsFactors = F)
    write.csv(df_RMSD,"df_RSMD_interactions_all.csv",row.names=F)
    
  }
  if (!dir.exists("RMSD_TEMP/")){dir.create("RMSD_TEMP/")}
    i<-1
    df_RMSD<-read.csv(paste0("df_RMSD_compare_interactions.csv"),stringsAsFactors =  F)
    df_RMSD<-df_RMSD%>%filter(number_models>quantile(df_RMSD$number_models,probs = 0.25))
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
  }
  for (name in 1:nrow(df_start)) {
    part<-paste0(parta,df_start$name[name])
    setwd(part)
    if (!dir.exists("RMSD_TEMP/")){dir.create("RMSD_TEMP/")}
    group_name<-v_group[1]
    v_name<-list.files("RMSD_TEMP/")
    df_interactions<-read.csv(paste0("RMSD_TEMP/",v_name[1]),stringsAsFactors = F)
    for (i in 2:length(v_name)) {
      df_interactions_add<-read.csv(paste0("RMSD_TEMP/",v_name[i]),stringsAsFactors = F)
      df_interactions<-rbind(df_interactions,df_interactions_add)
    }
    rm(df_interactions_add)
    df_interactions<-read.csv(paste0("df_interactions_fin.csv"),stringsAsFactors = F)
    df_interactions<-df_interactions%>%mutate(number_min=number.y)
    df_interactions<-df_interactions%>%mutate(number_max=number.x)
    df_interactions$number_min[df_interactions$number.x<df_interactions$number.y]<-df_interactions$number.x[df_interactions$number.x<df_interactions$number.y]
    df_interactions$number_max[df_interactions$number.x<df_interactions$number.y]<-df_interactions$number.y[df_interactions$number.x<df_interactions$number.y]
    df_interactions<-df_interactions%>%select(number_min,number_max)
    df_interactions<-df_interactions%>%filter(number_min<number_max)
    df_interactions<-df_interactions%>%mutate(compation=paste0(number_min,"_",number_max))
    df_interactions<-df_interactions%>%group_by(compation)%>%mutate(test=n())
    df_interactions<-unique(df_interactions)
    
    write.csv(df_interactions,"df_interactions_fin.csv",row.names=F)
  }
  for (name in 1:nrow(df_start)) {
    part<-paste0(parta,df_start$name[name])
    system(command = paste0("rm -r  ",part,"/RMSD_TEMP/"),ignore.stdout=T,wait = T)
    system(command = paste0("rm -r  ",part,"/compare_TEMP/"),ignore.stdout=T,wait = T)
  }
}
