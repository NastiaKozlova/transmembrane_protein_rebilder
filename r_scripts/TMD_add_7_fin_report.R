part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
w<-1
protein<-1
setwd(part)
v_RMSD<-5
for(w in 1:nrow(df_start)){  
  df_start_add<-read.csv(paste0(df_start$name[w],"/fin.csv"),stringsAsFactors = F)
  part_fin<-paste0(part,df_start$name[w],"/add_domain/")
  setwd(part_fin)

  for(protein in 1:nrow(df_start_add)){
    
    #    pdb_1<-read.pdb(paste0(df_start$name[w],"/fin_str/",protein,"__",df_start_add$best_model[protein]))
    df_start_add<-df_start_add%>%mutate(script=NA)
    
    df_RMSD_all<-read.csv(paste0(df_start_add$group_number[protein],"/df_RMSD_all.csv"),stringsAsFactors = F)
    df_RMSD_control<-df_RMSD_all%>%select(models.x)
    df_RMSD_control<-unique(df_RMSD_control)
    df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_RMSD)
    df_RMSD<-df_RMSD_all%>%filter(RMSD<v_RMSD)
    df_RMSD<-df_RMSD%>%group_by(models.x)%>%mutate(number=n())
    df_RMSD<-ungroup(df_RMSD)                                             
    df_RMSD<-df_RMSD%>%select(models.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD_add<-df_RMSD%>%filter(number==max(df_RMSD$number))
    df_RMSD<-df_RMSD%>%filter(number>2)
    df_RMSD<-rbind(df_RMSD,df_RMSD_add)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    df_RMSD_sep<-df_RMSD_all%>%filter(RMSD<v_RMSD)
    #df_RMSD_sep<-df_RMSD_sep%>%arrange(desc(number))
    if (!dir.exists(paste0(df_start_add$group_number[protein],"/din"))) {dir.create(paste0(df_start_add$group_number[protein],"/din"))}
    #    df_RMSD_control<-read.csv("df_RMSD.csv",stringsAsFactors =  F)
    i<-1
    for (i in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$models.x[i])) {
        df_RMSD_sep_test<-df_RMSD_sep%>%filter(models.x==df_RMSD$models.x[i])
        if (nrow(df_RMSD_sep_test)>10) {
          df_RMSD_sep_test_add<-data.frame(matrix(ncol=ncol(df_RMSD_sep_test),nrow = 1))
          colnames(df_RMSD_sep_test_add)<-colnames(df_RMSD_sep_test)
          df_RMSD_sep_test_add$models.x<-df_RMSD$models.x[i]
          df_RMSD_sep_test_add$models.y<-df_RMSD$models.x[i]
          df_RMSD_sep_test_add$RMSD<-0
          df_RMSD_sep_test<-rbind(df_RMSD_sep_test,df_RMSD_sep_test_add)
          write.csv(df_RMSD_sep_test,paste0(df_start_add$group_number[protein],"/din/grop_",i,".csv"),row.names = F)    
        }
      }
      
      df_RMSD_sep$models.x[df_RMSD_sep$models.x%in%df_RMSD_sep_test$models.y]<-NA
      df_RMSD_sep$models.x[df_RMSD_sep$models.y%in%df_RMSD_sep_test$models.y]<-NA
      df_RMSD_sep<-df_RMSD_sep%>%filter(!is.na(models.x))
      df_RMSD$models.x[df_RMSD$models.x%in%df_RMSD_sep_test$models.y]<-NA
    }
    
    v_groups<-list.files(paste0(df_start_add$group_number[protein],"/din/"))
    df_groups<-data.frame(matrix(ncol=7,nrow = length(v_groups)))
    colnames(df_groups)<-c("group_number","group_name","group_models","align_models","min_RMSD","max_RMSD","best_model")
    df_groups$group_number<-c(1:nrow(df_groups))
    df_groups$group_name<-v_groups
    i<-1 
    if (!dir.exists(paste0(df_start_add$group_number[protein],"/str"))) {dir.create(paste0(df_start_add$group_number[protein],"/str"))}
    if (!dir.exists(paste0(df_start_add$group_number[protein],"/fin_str"))) {dir.create(paste0(df_start_add$group_number[protein],"/fin_str"))}
    j<-1
    for (j in 1:nrow(df_groups)) {
      df_RMSD<-read.csv(paste0(df_start_add$group_number[protein],"/din/",df_groups$group_name[j]),stringsAsFactors = F)
      if (!dir.exists(paste0(df_start_add$group_number[protein],"/str/",j))) {dir.create(paste0(df_start_add$group_number[protein],"/str/",j))}
      for (i in 1:nrow(df_RMSD)) {
        pdb<-read.pdb(paste0(df_start_add$group_number[protein],"/structure/",df_RMSD$models.y[i]))
        write.pdb(pdb,paste0(df_start_add$group_number[protein],"/str/",j,"/",df_RMSD$models.y[i]))
      }
      df_RMSD_best<-df_RMSD%>%filter(RMSD==0)
      for (i in 1:nrow(df_RMSD_best)) {
        pdb<-read.pdb(paste0(df_start_add$group_number[protein],df_RMSD_best$chain.x[i],"/structure/",df_RMSD_best$models.x[i]))
        write.pdb(pdb,paste0(df_start_add$group_number[protein],"/fin_str/",j,"_",df_RMSD_best$models.x[i]))
      }
      df_groups$best_model[j]<-df_RMSD$models.x[1]
      df_groups$group_number[j]<-j
      df_groups$min_RMSD[j]<-NA
      df_groups$max_RMSD[j]<-NA
      df_groups$align_models[j]<-NA
    }
    #  df_groups<-df_groups
    df_groups<-df_groups%>%mutate(name=df_start$name[w])
    write.csv(df_groups,file = paste0(df_start_add$group_number[protein],"/fin.csv"),row.names = F)       
  }
  df_groups<-read.csv(paste0(df_start_add$group_number[1],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){
    df_groups_add<-read.csv(paste0(df_start_add$group_number[protein],"/fin.csv"),stringsAsFactors = F)
    df_groups<-rbind(df_groups,df_groups_add)
  }
  write.csv(df_groups,file = paste0("fin.csv"),row.names = F) 
}
