#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
path_to_PatchDock<-paste0(part_start,"programs/PatchDock/")
#pach_dock_repeats<-1
library(dplyr)
library(bio3d)
library(readr)
start<-read.pdb(paste0(part_start,"start/start.pdb"))
part<-paste0(part_start,"structure_prediction/")
df_start<-read.csv(paste0(part_start,"start/df_parts.csv"),stringsAsFactors = F)
df_start<-df_start%>%mutate(script=NA)

if(!dir.exists(part))    { dir.create(part)}
setwd(part)
chains<-unique(start$atom$chain)
if (length(chains)>1){
  df_chains<-data.frame(matrix(ncol=2,nrow = length(chains)))
  colnames(df_chains)<-c("chain","RMSD")
  df_chains$chain<-chains
  df_chains_all<-full_join(df_chains,df_chains,by="RMSD")
  df_chains_all<-df_chains_all%>%filter(chain.x!=chain.y)
  i<-1
  for (i in 1:nrow(df_chains_all)) {
    pdb_1.int<-atom.select(start,chain = df_chains_all$chain.x[i])
    pdb_2.int<-atom.select(start,chain = df_chains_all$chain.x[i])
    df_chains_all$RMSD[i]<-rmsd(a = start,b = start,a.inds = pdb_1.int,b.inds = pdb_2.int,fit = T)
  }
  df_chains_all<-df_chains_all%>%filter(RMSD>1)
  if (nrow(df_chains_all)==0){
    chains<-chains[1]
  }else(print("check your pdb structure, there are multiple different chains"))
}
pdb_1.int<-atom.select(start,chain = chains)
pdb_1<-trim(start,pdb_1.int)
i<-1
for (i in 1:nrow(df_start)) {
  if(!dir.exists(paste0(df_start$name[i])))    { dir.create(paste0(df_start$name[i]))}
  if(!dir.exists(paste0(df_start$name[i],"/patchdock")))    { dir.create(paste0(df_start$name[i],"/patchdock"))}
  write.pdb(pdb_1,paste0(df_start$name[i],"/control.pdb"))
  resno_1<-df_start$first_part_start[i]:df_start$first_part_finish[i]
  resno_2<-df_start$second_part_start[i]:df_start$second_part_finish[i]
  if(length(resno_1)<length(resno_2)){
    min_resno<-resno_1
    max_resno<-resno_2
  }else{
    min_resno<-resno_2
    max_resno<-resno_1
  }
  pdb_s.int<-atom.select(pdb_1,resno = max_resno)
  pdb_f.int<-atom.select(pdb_1,resno = min_resno)
  pdb_s<-trim(pdb_1,pdb_s.int)
  pdb_f<-trim(pdb_1,pdb_f.int)
  #  TMD7<-trim(pdb_1,TMD7.int)
  pdb_f$atom$chain <- "A"
  pdb_s$atom$chain <- "A"
  write.pdb(pdb_f,paste0(df_start$name[i],"/ligand.pdb"))
  write.pdb(pdb_s,paste0(df_start$name[i],"/receptor.pdb"))
  
  df_start$script[i]<-paste0("cd ",part,df_start$name[i],"/patchdock/\n",
                             path_to_PatchDock,"buildParamsFine.pl ","../receptor.pdb ", "../ligand.pdb 2.0 EI\n",
                             path_to_PatchDock,"patch_dock.Linux params.txt out.txt\n")
  system(command=df_start$script[i],ignore.stdout=T,wait = T,intern=F,ignore.stderr = T,show.output.on.console = F,minimized = T)
  df_out<-read_table(paste0(part,df_start$name[i],"/patchdock/out.txt"),skip = 26,skip_empty_rows = T,col_names = F)
  df_out<-df_out%>%select(X1, X3, X5, X7, X9, X11, X13, X15, X17, X19, X21, X23)
  colnames(df_out)<-c("number", "score", "pen.", "Area", "as1", "as2", "as12", "ACE", "hydroph", "Energy", "cluster", "dist.")
  n_structure<-nrow(df_out)
  a<-paste0("cd ",part,df_start$name[i],"/patchdock/\n",
                             path_to_PatchDock,"transOutput.pl out.txt 0 ",n_structure,"\n",
                             "rm out.txt\n",
                             "rm params.txt \n rm patch_dock.log\nrm cdrs3\n",collapse = "")
  system(command=a,ignore.stdout=T,wait = T,intern=F,ignore.stderr = T)
  write.csv(df_out,paste0(part,df_start$name[i],"/out.csv"),row.names = F)
}
