#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
df_topology<-read.csv("start/df_topology.csv",stringsAsFactors = F)
name<-1
v_RMSD<-5
persent_rigth<-3
for (name in 1:nrow(df_start)) {
  parta<-paste0(part_start,"structure_prediction/",df_start$name[name],"/")
  setwd(parta)
  if(file.exists("df_complex_energy.csv")){
    df_energy<-read.csv("df_complex_energy.csv",stringsAsFactors = F)
    #read list of all models structures
    df_RMSD_control<-read.csv("df_RMSD.csv",stringsAsFactors =  F)
    df_RMSD_control<-df_RMSD_control%>%mutate(number=NA)
    for (i in 1:nrow(df_RMSD_control)) {
      df_RMSD_control$number[i]<-strsplit(df_RMSD_control$models[i],split = ".",fixed = T)[[1]][3]
    }
    df_RMSD_control<-df_RMSD_control%>%mutate(number=as.numeric(number))
    if(file.exists("df_RMSD_all.csv")){
      df_RMSD_all<-read.csv("df_RMSD_all.csv",stringsAsFactors = F)
      df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_RMSD)
      
      
      df_RMSD_control<-df_RMSD_control[df_RMSD_control$number %in% unique(c(df_RMSD_all$number.x,df_RMSD_all$number.y)),]
      df_add<-df_RMSD_all%>%select(number.y,number.x,RMSD)
      colnames(df_add)<-colnames(df_RMSD_all)
      
      df_RMSD<-rbind(df_RMSD_all,df_add)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%group_by(number.x)%>%mutate(number=n())
      df_RMSD<-ungroup(df_RMSD)                                             
      df_RMSD_sep<-df_RMSD
      df_RMSD_sep<-df_RMSD_sep%>%filter(RMSD<v_RMSD)
      df_RMSD_sep<-df_RMSD_sep%>%arrange(desc(number))
      
      df_RMSD<-df_RMSD%>%select(number.x,number)
      df_RMSD<-unique(df_RMSD)
      
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      if (dir.exists("din")) {system(command = paste0("rm -r ",paste0(parta,"din/")),ignore.stdout=T,wait = T)}
      if (!dir.exists("din")) {dir.create("din")}
      
      i<-1
      for (i in 1:nrow(df_RMSD)) {
        if (!is.na(df_RMSD$number.x[i])) {
          df_RMSD_sep_test<-df_RMSD_sep%>%filter(number.x==df_RMSD$number.x[i])
          df_RMSD_sep_test_add<-data.frame(matrix(ncol=ncol(df_RMSD_sep_test),nrow = 1))
          colnames(df_RMSD_sep_test_add)<-colnames(df_RMSD_sep_test)
          df_RMSD_sep_test_add$number.x<-df_RMSD$number.x[i]
          df_RMSD_sep_test_add$number.y<-df_RMSD$number.x[i]
          df_RMSD_sep_test_add$RMSD<-0
          df_RMSD_sep_test<-rbind(df_RMSD_sep_test,df_RMSD_sep_test_add)
          df_RMSD_sep_test<-left_join(df_RMSD_sep_test,df_RMSD_control,by=c("number.y"="number"))
          colnames(df_RMSD_sep_test)<-c("number.x", "number.y", "RMSD_models", "number",
                                        "models", "RMSD_control", "max_length", "min_length")
          if(nrow(df_RMSD_sep_test)>2 ){
            write.csv(df_RMSD_sep_test,paste0("din/grop_",i,".csv"),row.names = F)  
          }
        }
        df_RMSD_sep$number.x[df_RMSD_sep$number.x%in%df_RMSD_sep_test$number.y]<-NA
        df_RMSD_sep$number.x[df_RMSD_sep$number.y%in%df_RMSD_sep_test$number.y]<-NA
        df_RMSD_sep<-df_RMSD_sep%>%filter(!is.na(number.x))
        df_RMSD$number.x[df_RMSD$number.x%in%df_RMSD_sep_test$number.y]<-NA
      }
      
      v_groups<-list.files(paste0(parta,"din"))
      df_groups<-data.frame(matrix(ncol=7,nrow = length(v_groups)))
      colnames(df_groups)<-c("group_number","group_name","group_models","align_models","min_RMSD","max_RMSD","best_model")
      if(nrow(df_groups)>0){
        df_groups$group_number<-c(1:nrow(df_groups))
        df_groups$group_name<-v_groups
        if (dir.exists(paste0(parta,"fin_str/"))) {system(command = paste0("rm -r ",paste0(parta,"fin_str/")),ignore.stdout=T,wait = T)}
        if (!dir.exists(paste0(parta,"fin_str/"))) {dir.create(paste0(parta,"fin_str/"))}
        j<-1
        for (j in 1:nrow(df_groups)) {
          
          df_RMSD<-read.csv(paste0("din/",df_groups$group_name[j]),stringsAsFactors = F)
          df_groups$min_RMSD[j]<-min(df_RMSD$RMSD_control)
          df_groups$max_RMSD[j]<-max(df_RMSD$RMSD_control)
          df_groups$group_models[j]<-nrow(df_RMSD)
          df_RMSD_align<-df_RMSD%>%filter(RMSD_control<v_RMSD)
          df_groups$align_models[j]<-nrow(df_RMSD_align)
          df_RMSD<-df_RMSD%>%filter(number.x==number.y)
          df_groups$best_model[j]<-df_RMSD$number.x[1]
          df_groups$group_number[j]<-j
          
        }
        df_groups<-df_groups%>%mutate(name=df_start$name[name])
        df_groups<-left_join(df_groups,df_RMSD_control,by=c("best_model"="number"))
        df_groups<-left_join(df_groups,df_energy)#,by=c("models", "max_length",   "min_length","RMSD" ))
        write.csv(df_groups,file = paste0("fin_TEMP.csv"),row.names = F)   
      }
      
      if(file.exists(paste0("fin_TEMP.csv"))){
        df_groups<-read.csv(file = paste0("fin_TEMP.csv"),stringsAsFactors =  F)  

        df_groups<-df_groups%>%filter(group_models>5)
        df_groups<-df_groups%>%filter(group_models>nrow(df_RMSD_control)/1000)
        df_groups<-df_groups%>%filter(group_models>=quantile(df_groups$group_models,0.975))
        df_groups<-df_groups%>%arrange(desc(group_models))
        
        v_length_first<-(df_start$first_part_finish[name]-df_start$first_part_start[name]+1)
        v_length_second<-(df_start$second_part_finish[name]-df_start$second_part_start[name]+1)
        if (v_length_first<v_length_second) {
          v_moving_start<-df_start$first_part_start[name]
          v_moving_finish<-df_start$first_part_finish[name]
        }else{
          v_moving_start<-df_start$second_part_start[name]
          v_moving_finish<-df_start$second_part_finish[name]
        }
        df_topology_TEST<-df_topology%>%filter(seq_beg>=v_moving_start)
        if (df_topology_TEST$seq_beg[1]>min(df_topology$seq_beg)) {
          df_topology_add<-df_topology%>%filter(seq_end==(min(df_topology_TEST$seq_beg)-1))
          df_topology_TEST<-rbind(df_topology_add,df_topology_TEST)

        }
        df_topology_TEST<-df_topology_TEST%>%filter(seq_end<=v_moving_finish)
        df_topology_TEST<-df_topology_TEST%>%mutate(symbol=(1))
        for (i in 2:nrow(df_topology_TEST)) {
          if(df_topology_TEST$type[(i-1)]!="Cytoplasmic"){df_topology_TEST$symbol[i]<-(-1)}
        }
        df_topology_TEST<-df_topology_TEST%>%filter(type=="Transmembrane")
        df_groups<-df_groups%>%mutate(angle=NA)
        pdb<-read.pdb(paste0("control.pdb"))
        df_pdb<-pdb$atom
        df_pdb<-df_pdb%>%filter(elety=="CA")
        df_topology_TEST<-df_topology_TEST%>%mutate(x=NA)
        df_topology_TEST<-df_topology_TEST%>%mutate(y=NA)
        df_topology_TEST<-df_topology_TEST%>%mutate(z=NA)
        for (p in 1:nrow(df_topology_TEST)) {
          df_A1<-df_pdb%>%filter(resno==df_topology_TEST$seq_beg[p])
          df_A2<-df_pdb%>%filter(resno==df_topology_TEST$seq_end[p])
          v_length<-sqrt((df_A1$x[1]-df_A2$x[1])^2+(df_A1$y[1]-df_A2$y[1])^2+(df_A1$z[1]-df_A2$z[1])^2)
          df_topology_TEST$x[p]<-(df_A1$x[1]-df_A2$x[1])/v_length*df_topology_TEST$symbol[p]
          df_topology_TEST$y[p]<-(df_A1$y[1]-df_A2$y[1])/v_length*df_topology_TEST$symbol[p]
          df_topology_TEST$z[p]<-(df_A1$z[1]-df_A2$z[1])/v_length*df_topology_TEST$symbol[p]
        }
        v_control_x<-mean(df_topology_TEST$x)
        v_control_y<-mean(df_topology_TEST$y)
        v_control_z<-mean(df_topology_TEST$z)
        i<-1
        for (i in 1:nrow(df_groups)) {
          pdb<-read.pdb(paste0("structure/",df_groups$models[i]))
          write.pdb(pdb,paste0("fin_str/",df_groups$group_number[i],".pdb"))
          df_pdb<-pdb$atom
          df_pdb<-df_pdb%>%filter(elety=="CA")
          p<-1
          df_topology_TEST<-df_topology_TEST%>%mutate(x=NA)
          df_topology_TEST<-df_topology_TEST%>%mutate(y=NA)
          df_topology_TEST<-df_topology_TEST%>%mutate(z=NA)
          for (p in 1:nrow(df_topology_TEST)) {
            df_A1<-df_pdb%>%filter(resno==df_topology_TEST$seq_beg[p])
            df_A2<-df_pdb%>%filter(resno==df_topology_TEST$seq_end[p])
            v_length<-sqrt((df_A1$x[1]-df_A2$x[1])^2+(df_A1$y[1]-df_A2$y[1])^2+(df_A1$z[1]-df_A2$z[1])^2)
            df_topology_TEST$x[p]<-(df_A1$x[1]-df_A2$x[1])/v_length*df_topology_TEST$symbol[p]
            df_topology_TEST$y[p]<-(df_A1$y[1]-df_A2$y[1])/v_length*df_topology_TEST$symbol[p]
            df_topology_TEST$z[p]<-(df_A1$z[1]-df_A2$z[1])/v_length*df_topology_TEST$symbol[p]
          }
          v_topology_x<-mean(df_topology_TEST$x)
          v_topology_y<-mean(df_topology_TEST$y)
          v_topology_z<-mean(df_topology_TEST$z)
          v_tost<-v_topology_x*v_control_x+v_topology_y*v_control_y+v_topology_z*v_control_z
          v_tost_mem<-v_topology_x*0+v_topology_y*0+v_topology_z*1
          df_groups$angle[i]<-acos(v_tost)
          df_groups$angle_mem[i]<-asin(v_tost_mem)
        }
        df_groups<-df_groups%>%mutate(angle=angle*90/pi*2)
        df_groups<-df_groups%>%mutate(orientarion="between")
        df_groups$orientarion[abs(df_groups$angle)<45]<-"as WT"
        df_groups$orientarion[abs(df_groups$angle)>135]<-"inverted"
        df_groups$orientarion[df_groups$orientarion=="as WT"&df_groups$RMSD<10]<-"WT"
        df_groups<-df_groups%>%select(name,group_number,group_models,align_models,name,      
                                      RMSD,bond_energy,bond_energy_fs,angle,angle_mem,
                                      orientarion,min_RMSD,max_RMSD)
        write.csv(df_groups,file = paste0("fin.csv"),row.names = F) 
      }
    }
  }
}