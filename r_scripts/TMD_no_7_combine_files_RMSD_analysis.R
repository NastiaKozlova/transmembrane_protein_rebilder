part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
v_RMSD<-10
name<-1
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name])
  setwd(part)
  if (!dir.exists("RMSD/")){dir.create("RMSD/")}
  print(paste0("start ",Sys.time()))
  v_RMSD_files<-list.files(paste0("RMSD/"))
  if(length(v_RMSD_files)>0){
    df_RMSD<-read.csv(paste0("RMSD/",v_RMSD_files[1]),stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%filter(RMSD<v_RMSD)
    for (i in 2:length(v_RMSD_files)) {
        df_RMSD_add<-read.csv(paste0("RMSD/",v_RMSD_files[i]),stringsAsFactors = F)
        df_RMSD_add<-df_RMSD_add%>%filter(RMSD<v_RMSD)
        df_RMSD<-rbind(df_RMSD,df_RMSD_add)
    
    }
    print(paste0("finish ",Sys.time()))
    write.csv(df_RMSD,"df_RMSD_all.csv",row.names = F)
  }
}
