print(Sys.time())
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1

for (name in 1:nrow(df_start)) { 
  part<-paste0(parta,df_start$name[name],"/")
  setwd(part) 
  
  df_RMSD<-read.csv(paste0("df_RMSD_compare_interactions.csv"),stringsAsFactors = F)
  df_RMSD<-df_RMSD%>%mutate(number=NA)
  for (i in 1:nrow(df_RMSD)) {
    df_RMSD$number[i]<-strsplit(df_RMSD$models[i],split = ".",fixed = T)[[1]][3]
  }
  df_RMSD<-df_RMSD%>%select(models,number)
  df_RMSD<-df_RMSD%>%mutate(number=as.numeric(number))
  i<-1
  df_interactions<-read.csv(paste0("compare_TEMP/",df_RMSD$number[i],".csv"),stringsAsFactors = F)
  for (i in 2:nrow(df_RMSD)) {
    df_interactions_add<-read.csv(paste0("compare_TEMP/",df_RMSD$number[i],".csv"),stringsAsFactors = F)
    df_interactions<-rbind(df_interactions,df_interactions_add)
  }
  wirte.csv(df_interactions,"df_interactions_fin.csv",row.names=F)
}
print(Sys.time())


