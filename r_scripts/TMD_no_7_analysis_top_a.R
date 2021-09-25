part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
q<-1
p<-1
v_RMSD<-5
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name])
  setwd(part)
  df_RMSD<-read.csv("df_RMSD_start.csv",stringsAsFactors = F)
  df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
  df_RMSD<-df_RMSD%>%select(models,RMSD)
  df_RMSD_start<-df_RMSD%>%filter(is.na(models))
  df_RMSD_start<-df_RMSD_start%>%filter(!is.na(models))
  df_RMSD_start<-full_join(df_RMSD_start,df_RMSD_start,by="RMSD")
  df_RMSD_test<-df_RMSD
  q<-1
#  v_reject_models<-c()
#  if(file.exists("df_RMSD_TEMP.csv")){
#    df_RMSD_start<-read.csv("df_RMSD_TEMP.csv",stringsAsFactors = F)
#    v_model<-unique(df_RMSD_start$models.y)
#    for (q in 1:length(v_model)) {
#      df_reject_models<-df_RMSD_start[df_RMSD_start$models.y==df_RMSD$models[q],] 
#      if((nrow(df_reject_models)+q-1)<nrow(df_RMSD)/100){
#        v_reject_models<-c(v_reject_models,df_RMSD$models[q])
#      }
#    }
#    df_RMSD_test<-df_RMSD_test[!df_RMSD_test$models%in%df_RMSD_start$models.y,]
#  }else{
#    (v_model<-c())
#  }
#    v_reject_models
#    if (nrow)
  if(file.exists("df_RMSD_TEMP.csv")){
    df_RMSD_start<-read.csv("df_RMSD_TEMP.csv",stringsAsFactors = F)
    v_model<-df_RMSD_test$models[df_RMSD_test$models%in%df_RMSD_start$models.y]
    df_RMSD_test<-df_RMSD_test[!df_RMSD_test$models%in%df_RMSD_start$models.y,]
  }else{
    (v_model<-c())
  }  
  q<-1
  for (q in 1:nrow(df_RMSD)) {
    df_RMSD_TEMP<-df_RMSD_test%>%filter(models==df_RMSD_test$models[q])
    df_RMSD_add<-full_join(df_RMSD,df_RMSD_TEMP,by="RMSD")
    v_model<-unique(c(v_model,df_RMSD_test$models[1:q]))
    df_RMSD_add<-df_RMSD_add%>%filter(models.x!=models.y)
    df_RMSD_add<-df_RMSD_add[!df_RMSD_add$models.x%in%v_model,]
#    df_RMSD_add<-df_RMSD_add[!df_RMSD_add$models.x%in%v_reject_models,]

    pdb_1<-read.pdb(paste0("structure/",df_RMSD_TEMP$models[1]))
    for (p in 1:nrow(df_RMSD_add)) {
      pdb_2<-read.pdb(paste0("structure/",df_RMSD_add$models.x[p]))
      df_RMSD_add$RMSD[p]<-rmsd(a=pdb_1,b=pdb_2,fit=T)
    }
    df_RMSD_add<-df_RMSD_add%>%filter(RMSD<v_RMSD) 
#    df_RMSD_add_test<-df_RMSD_add%>%filter(RMSD<v_RMSD) 
    
#    if(nrow (df_RMSD_add)>(nrow(df_RMSD))/1000){
#      v_reject_models<-c(v_reject_models,unique(df_RMSD_add$))
#    }

#    if(nrow (df_RMSD_add_test)<(nrow(df_RMSD_add))/1000){
#      df_RMSD_test<-df_RMSD_test%>%filter(models!=)
#    }
#    df_RMSD_add<-df_RMSD_add_test

    df_RMSD_start<-rbind(df_RMSD_start,df_RMSD_add)
    write.csv(df_RMSD_start,"df_RMSD_TEMP.csv",row.names = F)
  }
  df_RMSD_add<-df_RMSD_start%>%select(models.y,RMSD,models.x)
  colnames(df_RMSD_add)<-colnames(df_RMSD_start)
  df_RMSD_start<-rbind(df_RMSD_start,df_RMSD_add)
  write.csv(df_RMSD_start,"df_RMSD_all.csv",row.names = F)
}
