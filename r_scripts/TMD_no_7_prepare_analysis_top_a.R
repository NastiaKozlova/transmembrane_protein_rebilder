part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name])
  setwd(part)
  models<-list.files("structure")
  df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
  df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
  #  df_RMSD<-df_RMSD%>%mutate(number=NA)
  for (i in 1:nrow(df_RMSD)) {
    df_RMSD$models[i]<-strsplit(df_RMSD$models[i],split = ".",fixed = T)[[1]][3]
  }
  df_RMSD<-df_RMSD%>%mutate(models=as.numeric(models))
  df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
  df_RMSD_all<-df_RMSD_all%>%filter(models.x<models.y)
  if(file.exists("df_RMSD_all.csv")){
    df_RMSD_all<-read.csv("df_RMSD_all.csv",stringsAsFactors=T)
  }else{
    part<-paste0(parta,df_start$name[name])
    setwd(part)
    models<-list.files("structure")
    df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
    #  df_RMSD<-df_RMSD%>%mutate(number=NA)
    for (i in 1:nrow(df_RMSD)) {
      df_RMSD$models[i]<-strsplit(df_RMSD$models[i],split = ".",fixed = T)[[1]][3]
    }
    df_RMSD<-df_RMSD%>%mutate(models=as.numeric(models))
    df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
    df_RMSD_all<-df_RMSD_all%>%filter(models.x<models.y)
    write.csv(df_RMSD_all,"df_RMSD_all.csv",row.names = F)}
}