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
w<-1
protein<-1
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

i<-1
for(w in 1:nrow(df_start)){
  df_start_add<-read.csv(paste0(df_start$name[w],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){

    pdb_1<-read.pdb(paste0(df_start$name[w],"/fin_str/",protein,"_",df_start_add$models[protein]))
    df_start_add<-df_start_add%>%mutate(script=NA)
    part_fin<-paste0(part,df_start$name[w],"/add_domain/")
    if(!dir.exists(paste0(part_fin)))    { dir.create(paste0(part_fin))}
    if(!dir.exists(paste0(part_fin,df_start_add$group_number[protein])))    { dir.create(paste0(part_fin,df_start_add$group_number[protein]))}
    if(!dir.exists(paste0(part_fin,df_start_add$group_number[protein],"/patchdock/")))    { dir.create(paste0(part_fin,df_start_add$group_number[protein],"/patchdock"))}
    write.pdb(pdb_1,paste0(part_fin,df_start_add$group_number[protein],"/receptor.pdb"))
    min_resno<-df_start$third_part_start[w]:df_start$third_part_finish[w]
    pdb_f.int<-atom.select(start,resno = min_resno,chain = chains)
    pdb_f<-trim(start,pdb_f.int)
    pdb_f$atom$chain <- "Z"
    write.pdb(pdb_f,paste0(part_fin,df_start_add$group_number[protein],"/ligand.pdb"))
    df_start_add$script[protein]<-paste0("cd ",part_fin,df_start_add$group_number[protein],"/patchdock/\n",
                                         path_to_PatchDock,"buildParamsFine.pl ","../receptor.pdb ", "../ligand.pdb 2.0 EI\n",
                                         path_to_PatchDock,"patch_dock.Linux params.txt out.txt\n")
    system(command=df_start_add$script[protein],ignore.stdout=T,wait = T,intern=F,ignore.stderr = T,show.output.on.console = F,minimized = T)
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
