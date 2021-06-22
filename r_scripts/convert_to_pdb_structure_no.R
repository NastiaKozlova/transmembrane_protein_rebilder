part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
i<-1
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name])
  setwd(part)
  if(file.exists("patchdock/params.txt")){file.remove("patchdock/params.txt")}
  models<-list.files("patchdock")
  df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
  colnames(df_RMSD)<-c("models","RMSD")
  df_RMSD$models<-models
  #  df_RMSD<-df_RMSD%>%mutate(chain=model_name[name])
  if (!dir.exists("structure/")){dir.create("structure/")}
  if (!dir.exists("interactions/")){dir.create("interactions/")}
  for (i in 1:nrow(df_RMSD)) {
    pdb<-read.pdb(paste0("patchdock/",df_RMSD$models[i]))
    df_pdb<-pdb$atom
    df_pdb<-df_pdb%>%filter(elety=="CA")
    
    df_pdb_test<-df_pdb%>%select(resno,elety, x, y, z)
    
    v_pdb<-(df_pdb$resno)
    df_pdb<-data.frame(matrix(ncol = 2,nrow = length(v_pdb)))
    colnames(df_pdb)<-c("start","finish")
    df_pdb$start<-v_pdb
    df_pdb$finish<-v_pdb
    for (j in 2:nrow(df_pdb)) {
      if ((df_pdb$finish[j-1]+1)==(df_pdb$finish[j])){
        df_pdb$start[j]<-df_pdb$start[j-1]
        df_pdb$start[j-1]<-NA
      }
    }
    df_pdb<-df_pdb%>%filter(!is.na(start))
    df_pdb<-df_pdb%>%arrange(start)
    df_pdb_TEMP<-df_pdb
    for (j in 1:(nrow(df_pdb))) {
      df_pdb_TEMP$finish[j]<-df_pdb$start[j]
    }
    for (j in 2:(nrow(df_pdb))) {
      df_pdb_TEMP$start[j]<-df_pdb$finish[j-1]
    }
    df_pdb_TEMP<-df_pdb_TEMP%>%arrange(start)
    df_pdbt<-df_pdb_TEMP%>%filter(finish>1)
    max_nrow<-nrow(df_pdbt)
    df_pdbt<-left_join(df_pdbt,df_pdb_test,by=c("start"="resno"))
    df_pdbt<-left_join(df_pdbt,df_pdb_test,by=c("finish"="resno"))
    df_pdbt<-df_pdbt%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
    df_pdbt<-df_pdbt%>%mutate(max_length=(finish-start-1)*3)
    df_pdbt<-df_pdbt%>%filter(length>2)
    df_pdbt<-df_pdbt%>%filter(length<max_length)

    df_RMSD$max_length[i]<-max(df_pdbt$length)
    df_RMSD$min_length[i]<-min(df_pdbt$length)
    if(max_nrow==nrow(df_pdbt)){
      pdb.int<-atom.select(pdb,resno = df_pdb$start[1]:df_pdb$finish[1])
      pdb_1<-trim.pdb(pdb,pdb.int)
      pdb.int<-atom.select(pdb,resno = df_pdb$start[2]:df_pdb$finish[2])
      pdb_2<-trim.pdb(pdb,pdb.int)
      pdb<-cat.pdb(pdb_1,pdb_2,renumber = F,rechain = F)
      write.pdb(pdb,paste0("structure/",df_RMSD$models[i]))
      bs1<-binding.site(a=pdb_1,b=pdb_2)
      bs2<-binding.site(a=pdb_2,b=pdb_1)
      bs<-unique(c(bs1$resnames,bs2$resnames))
      df_interactions<-data.frame(matrix(ncol=1,nrow = length(bs)))
      colnames(df_interactions)<-"full_amino_name"
      df_interactions$full_amino_name<-bs
      write.csv(df_interactions,paste0("interactions/",df_RMSD$models[i]),row.names = F)
    }
  }
  
  df_RMSD<-df_RMSD%>%filter(min_length>2)
  write.csv(df_RMSD,"df_RMSD_start.csv",row.names = F)
  df_RMSD<-read.csv("df_RMSD_start.csv",stringsAsFactors = F)
  pdb_com<-read.pdb("control.pdb")
  models<-list.files("structure")
  for (i in 1:nrow(df_RMSD)) {
    if (file.exists(paste0("structure/",df_RMSD$models[i]))){
      

    pdb<-read.pdb(paste0("structure/",df_RMSD$models[i])) 
    df_pdb<-pdb$atom
    v_pdb<-unique(df_pdb$resno)
    pdb.int<-atom.select(pdb,resno = v_pdb,elety="CA")
    pdb<-trim.pdb(pdb,pdb.int)
    pdb_com.int<-atom.select(pdb_com,resno = v_pdb,elety="CA")
    pdb_sep<-trim.pdb(pdb_com,pdb_com.int)
    df_RMSD$RMSD[i]<-rmsd(pdb_sep,pdb,fit=T)
    }
  }
  df_RMSD<-df_RMSD%>%filter(!is.na(RMSD))
  write.csv(df_RMSD,"df_RMSD.csv",row.names = F)
}
