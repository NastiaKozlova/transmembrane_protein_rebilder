part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
v_RMSD<-5
name<-1
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name])
  setwd(part)
  if (!dir.exists("RMSD/")){dir.create("RMSD/")}
  df_RMSD<-read.csv("df_interactions_fin.csv",stringsAsFactors = F)
  df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
  df_RMSD<-df_RMSD%>%select(number.x,number.y,RMSD)
  v_model<-unique(df_RMSD$number.x)
  b<-c()
  for (j in 1:length(v_model)) {
    if (!file.exists(paste0("RMSD/",v_model[j],".csv")) ){
      b<-c(b,v_model[j])
    }
  }
  v_model<-b
  v_model<-sort(v_model)
  df_RMSD_TEMP<-df_RMSD
  df_RMSD<-df_RMSD[df_RMSD$number.x%in%v_model,]
  if(length(v_model)>0){
    for (j in 1:length(v_model)) {
      df_RMSD_TEMP<-  df_RMSD%>%filter(number.x==v_model[j])
      pdb_1<-read.pdb(paste0("structure/out.txt.",v_model[j],".pdb"))
      pdb_1.int<-atom.select(pdb_1, "calpha")
      for (i in 1:nrow(df_RMSD_TEMP)) {
        pdb_2<-read.pdb(paste0("structure/out.txt.",df_RMSD_TEMP$number.y[i],".pdb"))
        pdb_2.int<-atom.select(pdb_2, "calpha")
        df_RMSD_TEMP$RMSD[i]<-rmsd(a=pdb_1,b=pdb_2,a.inds = pdb_1.int,b.inds = pdb_2.int,fit=T)
  
      }
      df_RMSD_TEMP<-df_RMSD_TEMP%>%filter(RMSD<10)
      write.csv(df_RMSD_TEMP,paste0("RMSD/",v_model[j],".csv"),row.names = F)
      df_RMSD<-  df_RMSD%>%filter(number.x!=v_model[j])
    }
  }
}
