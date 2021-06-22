#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
v_RMSD<-5
persent_rigth<-3
for (name in 1:nrow(df_start)) {
  parta<-paste0(part_start,"structure_prediction/",df_start$name[name],"/")
  setwd(parta)
  df_energy<-read.csv("df_complex_energy.csv",stringsAsFactors = F)
  #read list of all models structures
  df_RMSD_control<-read.csv("df_RMSD.csv",stringsAsFactors =  F)
  df_RMSD_control<-df_RMSD_control%>%mutate(number=NA)
  for (i in 1:nrow(df_RMSD_control)) {
    df_RMSD_control$number[i]<-strsplit(df_RMSD_control$models[i],split = ".",fixed = T)[[1]][3]
  }
  df_RMSD_control<-df_RMSD_control%>%mutate(number=as.numeric(number))
  df_RMSD_all<-read.csv("df_RMSD_all.csv",stringsAsFactors = F)
  df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_RMSD)
  
  
  df_RMSD_control<-df_RMSD_control[df_RMSD_control$number %in% unique(c(df_RMSD_all$number_min,df_RMSD_all$number_max)),]
  df_add<-df_RMSD_all%>%select(number_max,number_min,RMSD)
  colnames(df_add)<-colnames(df_RMSD_all)
  
  df_RMSD<-rbind(df_RMSD_all,df_add)
  df_RMSD<-unique(df_RMSD)
  df_RMSD<-df_RMSD%>%group_by(number_min)%>%mutate(number=n())
  df_RMSD<-ungroup(df_RMSD)                                             
  df_RMSD_sep<-df_RMSD
  df_RMSD_sep<-df_RMSD_sep%>%filter(RMSD<v_RMSD)
  df_RMSD_sep<-df_RMSD_sep%>%arrange(desc(number))
  
  df_RMSD<-df_RMSD%>%select(number_min,number)
  df_RMSD<-unique(df_RMSD)
  df_RMSD_add<-df_RMSD%>%filter(number==max(df_RMSD$number))
  df_RMSD<-df_RMSD%>%filter(number>2)
  df_RMSD<-rbind(df_RMSD,df_RMSD_add)
  df_RMSD<-unique(df_RMSD)
  df_RMSD<-df_RMSD%>%arrange(desc(number))
  df_RMSD<-df_RMSD%>%filter(number>(nrow(df_RMSD_control)/100))
  if (!dir.exists("din")) {dir.create("din")}

  i<-1
  for (i in 1:nrow(df_RMSD)) {
    if (!is.na(df_RMSD$number_min[i])) {
      df_RMSD_sep_test<-df_RMSD_sep%>%filter(number_min==df_RMSD$number_min[i])
      df_RMSD_sep_test_add<-data.frame(matrix(ncol=ncol(df_RMSD_sep_test),nrow = 1))
      colnames(df_RMSD_sep_test_add)<-colnames(df_RMSD_sep_test)
      df_RMSD_sep_test_add$number_min<-df_RMSD$number_min[i]
      df_RMSD_sep_test_add$number_max<-df_RMSD$number_min[i]
      df_RMSD_sep_test_add$RMSD<-0
      df_RMSD_sep_test<-rbind(df_RMSD_sep_test,df_RMSD_sep_test_add)
      df_RMSD_sep_test<-left_join(df_RMSD_sep_test,df_RMSD_control,by=c("number_max"="number"))
      colnames(df_RMSD_sep_test)<-c("number_min", "number_max", "RMSD_models", "number",
                                    "models", "RMSD_control", "max_length", "min_length")
      if(nrow(df_RMSD_sep_test)>(nrow(df_RMSD_control)/100*persent_rigth)){
        print(paste0(df_RMSD$number_min[i]," ",nrow(df_RMSD_sep_test)," ",(nrow(df_RMSD_control)/100*persent_rigth)))
        write.csv(df_RMSD_sep_test,paste0("din/grop_",i,".csv"),row.names = F)  
      }
    }
    df_RMSD_sep$number_min[df_RMSD_sep$number_min%in%df_RMSD_sep_test$number_max]<-NA
    df_RMSD_sep$number_min[df_RMSD_sep$number_max%in%df_RMSD_sep_test$number_max]<-NA
    df_RMSD_sep<-df_RMSD_sep%>%filter(!is.na(number_min))
    df_RMSD$number_min[df_RMSD$number_min%in%df_RMSD_sep_test$number_max]<-NA
  }
  
  v_groups<-list.files(paste0(parta,"din"))
  df_groups<-data.frame(matrix(ncol=7,nrow = length(v_groups)))
  colnames(df_groups)<-c("group_number","group_name","group_models","align_models","min_RMSD","max_RMSD","best_model")
  if(nrow(df_groups)>0){
    df_groups$group_number<-c(1:nrow(df_groups))
    df_groups$group_name<-v_groups
    i<-1 
#    if (!dir.exists(paste0(parta,"str"))) {dir.create(paste0(parta,"str"))}
    if (!dir.exists(paste0(parta,"fin_str/"))) {dir.create(paste0(parta,"fin_str/"))}
    j<-1
    for (j in 1:nrow(df_groups)) {

      df_RMSD<-read.csv(paste0("din/",df_groups$group_name[j]),stringsAsFactors = F)
#      if (!dir.exists(paste0(parta,"str/",j))) {dir.create(paste0(parta,"str/",j))}
      df_groups$min_RMSD[j]<-min(df_RMSD$RMSD_control)
      df_groups$max_RMSD[j]<-max(df_RMSD$RMSD_control)
      df_groups$group_models[j]<-nrow(df_RMSD)
      df_RMSD<-df_RMSD%>%filter(number_min==number_max)
      
      for (i in 1:nrow(df_RMSD)) {
        pdb<-read.pdb(paste0("structure/",df_RMSD$models[i]))
        write.pdb(pdb,paste0("fin_str/",j,"_",df_RMSD$models[i]))
      }
      df_groups$best_model[j]<-df_RMSD$number_min[1]
      df_groups$group_number[j]<-j
      df_RMSD<-df_RMSD%>%filter(RMSD_control<v_RMSD)
      df_groups$align_models[j]<-nrow(df_RMSD)
    }
    df_groups<-df_groups%>%mutate(name=df_start$name[name])
    df_groups<-left_join(df_groups,df_RMSD_control,by=c("best_model"="number"))
    df_groups<-left_join(df_groups,df_energy)#,by=c("models", "max_length",   "min_length","RMSD" ))
    write.csv(df_groups,file = paste0("fin.csv"),row.names = F)    
  }
}
for (i in 1:nrow(df_start)) {
  if(!file.exists(paste0(part_start,"structure_prediction/",df_start$name[1],"/fin.csv"))){
    df_start$name[i]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(name))
if(nrow(df_start)>0 ){
  df_groups<-read.csv(paste0(part_start,"structure_prediction/",df_start$name[1],"/fin.csv"),stringsAsFactors = F)
  if(nrow(df_start)>1){
    for (name in 2:nrow(df_start)) {
      df_groups_add<-read.csv(paste0(part_start,"structure_prediction/",df_start$name[name],"/fin.csv"),stringsAsFactors = F)
      df_groups<-rbind(df_groups,df_groups_add)
    }
  }
}
write.csv(df_groups,file = paste0(part_start,"structure_prediction/fin.csv"),row.names = F) 
for (name in 1:nrow(df_start)) {
  system(command = paste0("rm -r ",part_start,"structure_prediction/",df_start$name[name],"/compare_interaction"),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r ",part_start,"structure_prediction/",df_start$name[name],"/interactions"),ignore.stdout=T,wait = T)
}
