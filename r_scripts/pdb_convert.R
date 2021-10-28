#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
path_to_PatchDock<-paste0(part_start,"programs/PatchDock/")
pach_dock_repeats<-1
library(dplyr)
library(bio3d)
library(readr)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
setwd(part)
i<-1
df_start<-df_start%>%filter(!is.na(third_part_finish))
df_start<-df_start%>%filter(third_part_start>0)
df_start<-df_start%>%filter(third_part_finish>0)
for(w in 1:nrow(df_start)){
  if(!file.exists(paste0(df_start$name[w],"/fin.csv"))){
    df_start$name[w]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(name))

w<-1
for(w in 1:nrow(df_start)){
  if(file.exists(paste0(df_start$name[w],"/fin.csv"))){
    v_start<-length(list.files(paste0(df_start$name[1],"/structure")))
    df_start_add<-read.csv(paste0(df_start$name[w],"/fin.csv"),stringsAsFactors = F)
    df_start_add<-df_start_add%>%mutate(name=df_start$name[1])
    df_start_add<-left_join(df_start_add,df_start)
    df_start_add<-df_start_add%>%filter(group_models>5)
    df_start_add<-df_start_add%>%filter(group_models>v_start/1000)
    df_start_add<-df_start_add%>%filter(group_models>=quantile(df_start_add$group_models,0.975))
    print(paste0(w," ",df_start$name[w]," ",nrow(df_start_add)))
    if(nrow(df_start_add)>0){
      for(protein in 1:nrow(df_start_add)){
        
        pdb_1<-read.pdb(paste0(df_start$name[w],"/fin_str/",df_start_add$group_number[protein],".pdb"))
        df_start_add<-df_start_add%>%mutate(script=NA)
        part_fin<-paste0(part,df_start$name[w],"/add_domain/")
        if(!dir.exists(paste0(part_fin)))    { dir.create(paste0(part_fin))}
        if(!dir.exists(paste0(part_fin,df_start_add$group_number[protein])))    { dir.create(paste0(part_fin,df_start_add$group_number[protein]))}
        if(!dir.exists(paste0(part_fin,df_start_add$group_number[protein],"/patchdock/")))    { dir.create(paste0(part_fin,df_start_add$group_number[protein],"/patchdock"))}
        write.pdb(pdb_1,paste0(part_fin,df_start_add$group_number[protein],"/receptor.pdb"))
        
        min_resno<-df_start$third_part_start[w]:df_start$third_part_finish[w]
        pdb_f.int<-atom.select(start,resno = min_resno)
        pdb_f<-trim(start,pdb_f.int)
        pdb_f$atom$chain <- "Z"
        write.pdb(pdb_f,paste0(part_fin,df_start_add$group_number[protein],"/ligand.pdb"))
        df_start_add$script[protein]<-paste0("cd ",part_fin,df_start_add$group_number[protein],"/patchdock/\n",
                                             path_to_PatchDock,"buildParamsFine.pl ","../receptor.pdb ", "../ligand.pdb 2.0 EI\n",
                                             path_to_PatchDock,"patch_dock.Linux params.txt out.txt\n")
        system(command=df_start_add$script[protein],ignore.stdout=T,wait = T,intern=F,ignore.stderr = T)
        if(file.exists(paste0(part_fin,df_start_add$group_number[protein],"/patchdock/out.txt"))){
          df_out<-read_table(paste0(part_fin,df_start_add$group_number[protein],"/patchdock/out.txt"),skip = 26,skip_empty_rows = T,col_names = F)
          df_out<-df_out%>%select(X1, X3, X5, X7, X9, X11, X13, X15, X17, X19, X21, X23)
          colnames(df_out)<-c("number", "score", "pen.", "Area", "as1", "as2", "as12", "ACE", "hydroph", "Energy", "cluster", "dist.")
          n_structure<-nrow(df_out)
          a<-paste0("cd ",part_fin,df_start_add$group_number[protein],"/patchdock/\n",
                    path_to_PatchDock,"transOutput.pl out.txt 0 ",nrow(df_out),"\n",
                    "rm out.txt\n",
                    "rm params.txt \n rm patch_dock.log\n",collapse = "")
          system(command=a,ignore.stdout=T,wait = T,intern=F,ignore.stderr = T)
          write.csv(df_out,paste0(part_fin,df_start_add$group_number[protein],"/out.csv"),row.names = F)
        }
      }
    }
  }
}
