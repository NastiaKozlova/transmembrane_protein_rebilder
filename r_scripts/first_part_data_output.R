#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
setwd(part_start)
part<-paste0(part_start,"structure_prediction/")
start<-read.pdb(paste0(part_start,"start/start.pdb"))
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
#w<-1
#protein<-1
setwd(part)
for(w in 1:nrow(df_start)){
  if(!file.exists(paste0(df_start$name[w],"/fin.csv"))){
    df_start$name[w]<-NA
  }
}
df_start<-df_start%>%filter(!is.na(name))
df_start_all<-read.csv(paste0(df_start$name[1],"/fin.csv"),stringsAsFactors = F)
df_start_all<-df_start_all%>%mutate(name=df_start$name[1])
df_start_all<-left_join(df_start_all,df_start)
df_start_all<-df_start_all%>%mutate(frequence=group_models/sum(df_start_all$group_models)*100)
if(nrow(df_start)>1){
  for (i in 2:nrow(df_start)) {
    df_start_add<-read.csv(paste0(df_start$name[i],"/fin.csv"),stringsAsFactors = F)
    df_start_add<-df_start_add%>%mutate(name=df_start$name[i])
    df_start_add<-left_join(df_start_add,df_start)
    df_start_add<-df_start_add%>%mutate(frequence=group_models/sum(df_start_add$group_models)*100)
    df_start_all<-rbind(df_start_all,df_start_add)
  }
}
df_start_all<-df_start_all%>%mutate(RMSD=round(RMSD,digits = 2))
df_start_all<-df_start_all%>%mutate(frequence=round(frequence,digits = 1))
df_start_all<-df_start_all%>%mutate(angle=round(angle,digits = 0))
df_start_all<-df_start_all%>%mutate(persent_align=round(align_models/group_models*100,digits = 1))
df_start_all<-df_start_all%>%select(name,orientarion,RMSD,frequence, RMSD,persent_align,group_models, angle,
                                    first_part_model,first_part_start,first_part_finish,second_part_model,second_part_start, 
                                    second_part_finish,group_number,bond_energy,bond_energy_fs)


if(!dir.exists(paste0(part_start,"first_part/"))){dir.create(paste0(part_start,"first_part/"))}
if(!dir.exists(paste0(part_start,"results/first_part/"))){dir.create(paste0(part_start,"results/first_part/"))}
if(!dir.exists(paste0(part_start,"results/first_part/structure"))){dir.create(paste0(part_start,"results/first_part/structure"))}
i<-1
write.csv(df_start_all,paste0(part_start,"results/df_first_part.csv"),row.names = F)
for (i in 1:nrow(df_start_all)) {
  if(!dir.exists(paste0(part_start,"results/first_part/structure/",df_start_all$name[i]))){dir.create(paste0(part_start,"results/first_part/structure/",df_start_all$name[i]))}
  if(file.exists(paste0(part,df_start_all$name[i],"/fin_str/",df_start_all$group_number[i],".pdb"))){
    pdb<-read.pdb(paste0(part,df_start_all$name[i],"/fin_str/",df_start_all$group_number[i],".pdb"))
    write.pdb(pdb,paste0(part_start,"results/first_part/structure/",df_start_all$name[i],"/",df_start_all$name[i],"_",df_start_all$group_number[i],".pdb"))
  }else(print(i))
}

#i<-1
#for (i in 1:nrow(df_start_all)) {
#  name_dir<-paste0(part_start,"results/first_part,",df_start_all$name[i])
#  pdb<-read.pdb(paste0(part_start,"results/first_part/",df_start_all$name[i],"/",df_start_all$name[i],"_",df_start_all$group_number[i],".pdb"))
#  pdb.int<-atom.select(pdb,resno = df_start_all$first_part_start[i]:df_start_all$first_part_finish[i])
#  pdb_1<-trim.pdb(pdb,pdb.int)
#  pdb.int<-atom.select(pdb,resno = df_start_all$second_part_start[i]:df_start_all$second_part_finish[i])
#  pdb_2<-trim.pdb(pdb,pdb.int)
#  bs1<-binding.site(a=pdb_1,b=pdb_2)
#  bs2<-binding.site(a=pdb_2,b=pdb_1)
#  bs<-unique(c(bs1$resnames,bs2$resnames))
#  df_interactions<-data.frame(matrix(ncol=1,nrow = length(bs)))
#  colnames(df_interactions)<-"full_amino_name"
#  df_interactions$full_amino_name<-bs
  
#}
#/interactions/


#write.csv(df_interactions,paste0("interactions/",df_RMSD$models[i]),row.names = F)