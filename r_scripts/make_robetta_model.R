part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
v_RMSD<-5
df_start<-read.csv("start/df_name.csv",stringsAsFactors = F)
if (!dir.exists("MD_compatirion/")){dir.create("MD_compatirion/")}
if (!dir.exists("MD_compatirion/stabilisation/")){dir.create("MD_compatirion/stabilisation/")}
if (!dir.exists("MD_compatirion/stabilisation/start_pdb")){dir.create("MD_compatirion/stabilisation/start_pdb")}
if (!dir.exists("MD_compatirion/stabilisation/RMSD")){dir.create("MD_compatirion/stabilisation/RMSD")}
for (i in 1:nrow(df_start)){
  if (!dir.exists(paste0("MD_compatirion/stabilisation/start_pdb/",df_start$structure[i]))){dir.create(paste0("MD_compatirion/stabilisation/start_pdb/",df_start$structure[i]))}
  pdbsplit(pdb.files=paste0("start/robetta/robetta_models_",df_start$prediction_name[i],".pdb"),
           path=paste0("MD_compatirion/stabilisation/start_pdb/",df_start$structure[i]),multi=T)
}
v_groups<-list.files(paste0("MD_compatirion/stabilisation/start_pdb/"))
setwd(paste0("MD_compatirion/stabilisation/"))
q<-1
for(q in 1:length(v_groups)){
  v_chains<-list.files(paste0("start_pdb/",v_groups[q]))
  df_RMSD<-data.frame(matrix(ncol=2,nrow = length(v_chains)))
  colnames(df_RMSD)<-c('model',"RMSD")
  df_RMSD$model<-v_chains
  df_RMSD_control<-df_RMSD
  df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
  df_RMSD_all<-df_RMSD_all%>%filter(model.x!=model.y)
  for (j in 1:nrow(df_RMSD_all)) {
    pdb_1<-read.pdb(paste0("start_pdb/",v_groups[q],"/",df_RMSD_all$model.x[j]))
    pdb_2<-read.pdb(paste0("start_pdb/",v_groups[q],"/",df_RMSD_all$model.y[j]))
    df_RMSD_all$RMSD[j]<-rmsd(pdb_1,pdb_2,fit = T)
  }
  write.csv(df_RMSD_all,paste0("RMSD/",v_groups[q],".csv"),row.names = F)
  
  df_RMSD_all<-read.csv(paste0("RMSD/",v_groups[q],".csv"),stringsAsFactors = F)
  df_RMSD<-df_RMSD_all%>%filter(RMSD<v_RMSD)
  df_RMSD<-df_RMSD%>%group_by(model.x)%>%mutate(number=n())
  df_RMSD<-ungroup(df_RMSD)                                             
  df_RMSD<-df_RMSD%>%select(model.x,number)
  df_RMSD<-unique(df_RMSD)
  df_RMSD_add<-df_RMSD%>%filter(number==max(df_RMSD$number))
  df_RMSD<-df_RMSD%>%filter(number>2)
  df_RMSD<-rbind(df_RMSD,df_RMSD_add)
  df_RMSD<-unique(df_RMSD)
  df_RMSD<-df_RMSD%>%arrange(desc(number))
  df_RMSD_sep<-df_RMSD_all%>%filter(RMSD<v_RMSD)
  #df_RMSD_sep<-df_RMSD_sep%>%arrange(desc(number))
  if (!dir.exists("din")) {dir.create("din")}
  if (!dir.exists(paste0("din/",v_groups[q]))) {dir.create(paste0("din/",v_groups[q]))}
  #  df_RMSD_control<-read.csv("df_RMSD.csv",stringsAsFactors =  F)
  i<-1
  for (i in 1:nrow(df_RMSD)) {
    if (!is.na(df_RMSD$model.x[i])) {
      df_RMSD_sep_test<-df_RMSD_sep%>%filter(model.x==df_RMSD$model.x[i])
      
      if (nrow(df_RMSD_sep_test)>(nrow(df_RMSD_control)/10)) {
        #        print (paste(df_RMSD$model.x[i],i,nrow(df_RMSD_sep_test)))
        
        df_RMSD_sep_test_add<-data.frame(matrix(ncol=ncol(df_RMSD_sep_test),nrow = 1))
        colnames(df_RMSD_sep_test_add)<-colnames(df_RMSD_sep_test)
        df_RMSD_sep_test_add$model.x<-df_RMSD$model.x[i]
        df_RMSD_sep_test_add$model.y<-df_RMSD$model.x[i]
        df_RMSD_sep_test_add$RMSD<-0
        df_RMSD_sep_test<-rbind(df_RMSD_sep_test,df_RMSD_sep_test_add)
        df_RMSD_sep_test<-left_join(df_RMSD_sep_test,df_RMSD_control,by=c("model.y"="model"))
        colnames(df_RMSD_sep_test)<-c("model.x", "model.y", "RMSD_models", "RMSD_control"  )
        write.csv(df_RMSD_sep_test,paste0("din/",v_groups[q],"/grop_",i,".csv"),row.names = F)    
      }
    }
    
    df_RMSD_sep$model.x[df_RMSD_sep$model.x%in%df_RMSD_sep_test$model.y]<-NA
    df_RMSD_sep$model.x[df_RMSD_sep$model.y%in%df_RMSD_sep_test$model.y]<-NA
    df_RMSD_sep<-df_RMSD_sep%>%filter(!is.na(model.x))
    df_RMSD$model.x[df_RMSD$model.x%in%df_RMSD_sep_test$model.y]<-NA
  }
  if (!dir.exists(paste0("str"))) {dir.create(paste0("str"))}
  v_group<-list.files(paste0("din/",v_groups[q]))
  for (j in 1:length(v_group)) {
    df_RMSD<-read.csv(paste0("din/",v_groups[q],"/grop_",j,".csv"),stringsAsFactors = F)
    if (!dir.exists(paste0("str/",v_groups[q]))) {dir.create(paste0("str/",v_groups[q]))}
    #    if (!dir.exists(paste0("str/",v_groups[i],"/",j))) {dir.create(paste0("str/",v_groups[i],"/",j))}
    structure_name<-unique(df_RMSD$model.x)
    
    for (p in 1:length(structure_name)) {
      pdb<-read.pdb(paste0("start_pdb/",v_groups[q],"/",structure_name[p]))
      write.pdb(pdb,paste0("str/",v_groups[q],"/",j,".pdb"))
    }
    
  }
}
