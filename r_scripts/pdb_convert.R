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
df_start_ad<-read.csv(paste0(part_start,"results/df_first_part.csv"),stringsAsFactors = F)
df_start<-full_join(df_start_ad,df_start,by = c("name", "first_part_model",
                                                "first_part_start", "first_part_finish",
                                                "second_part_model", "second_part_start",
                                                "second_part_finish"))
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
  part_fin<-paste0(part,df_start$name[w],"/add_domain/")
  setwd(part_fin)
  if(!dir.exists(part_fin)){dir.create(part_fin)}
  pdb_1<-read.pdb(paste0(part_start,"results/first_part/structure/",df_start$name[w],"/",df_start$name[w],"_",df_start$group_number[w],".pdb"))
  df_start<-df_start%>%mutate(script=NA)
  if(!dir.exists(paste0(part_fin)))    { dir.create(paste0(part_fin))}
  if(!dir.exists(paste0(part_fin,df_start$group_number[w])))    { dir.create(paste0(part_fin,df_start$group_number[w]))}
  if(!dir.exists(paste0(part_fin,df_start$group_number[w],"/patchdock/")))    { dir.create(paste0(part_fin,df_start$group_number[w],"/patchdock"))}
  write.pdb(pdb_1,paste0(part_fin,df_start$group_number[w],"/receptor.pdb"))
  
  min_resno<-df_start$third_part_start[w]:df_start$third_part_finish[w]
  pdb_f.int<-atom.select(start,resno = min_resno)
  pdb_f<-trim(start,pdb_f.int)
  pdb_f$atom$chain <- "Z"
  write.pdb(pdb_f,paste0(part_fin,df_start$group_number[w],"/ligand.pdb"))
  df_start$script[w]<-paste0("cd ",part_fin,df_start$group_number[w],"/patchdock/\n",
                             path_to_PatchDock,"buildParamsFine.pl ","../receptor.pdb ", "../ligand.pdb 2.0 EI\n",
                             path_to_PatchDock,"patch_dock.Linux params.txt out.txt\n")
  system(command=df_start$script[w],ignore.stdout=T,wait = T,intern=F,ignore.stderr = T)
  if(file.exists(paste0(part_fin,df_start$group_number[w],"/patchdock/out.txt"))){
    df_out<-read_table(paste0(part_fin,df_start$group_number[w],"/patchdock/out.txt"),skip = 26,skip_empty_rows = T,col_names = F)
    df_out<-df_out%>%select(X1, X3, X5, X7, X9, X11, X13, X15, X17, X19, X21, X23)
    colnames(df_out)<-c("number", "score", "pen.", "Area", "as1", "as2", "as12", "ACE", "hydroph", "Energy", "cluster", "dist.")
    n_structure<-nrow(df_out)
    a<-paste0("cd ",part_fin,df_start$group_number[w],"/patchdock/\n",
              path_to_PatchDock,"transOutput.pl out.txt 0 ",nrow(df_out),"\n",
              "rm out.txt\n",
              "rm params.txt \n rm patch_dock.log\n",collapse = "")
    system(command=a,ignore.stdout=T,wait = T,intern=F,ignore.stderr = T)
    write.csv(df_out,paste0(part_fin,df_start$group_number[w],"/out.csv"),row.names = F)
  }
}


