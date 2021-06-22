part_start = commandArgs(trailingOnly=TRUE)
Sys.time()
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
i<-1
for (name in 1:nrow(df_start)) {
  print(Sys.time())
  part<-paste0(parta,df_start$name[name],"/")
  setwd(part)  
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
