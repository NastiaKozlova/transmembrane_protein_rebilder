part_start = commandArgs(trailingOnly=TRUE)
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
i<-1
for(w in 1:nrow(df_start)){
  df_start_add<-read.csv(paste0(df_start$name[w],"/fin.csv"),stringsAsFactors = F)
  for(protein in 1:nrow(df_start_add)){
    part_fin<-paste0(part,df_start$name[w],"/add_domain/",df_start_add$group_number[protein],"/")
    setwd(part_fin)
    models<-list.files("structure")
    df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
    df_RMSD_temp<-df_RMSD%>%mutate(RMSD=NA)
    df_RMSD_all<-full_join(df_RMSD_temp,df_RMSD_temp,by="RMSD")
    df_RMSD_all<-df_RMSD_all%>%filter(models.x!=models.y)
    df_RMSD_all<-df_RMSD_all%>%select(models.x,models.y,RMSD)
    if(file.exsit("df_RMSD_all.csv")){
      df_RMSD_all<-read.csv("df_RMSD_all.csv",stringsAsFactors=T)
    }
    v_test<-length(df_RMSD_all$RMSD[!is.na(df_RMSD_all$RMSD)])
    for (q in (v_test+1):nrow(df_RMSD_all)) {
      pdb_1<-read.pdb(paste0("structure/",df_RMSD_all$models.x[q]))
      pdb_2<-read.pdb(paste0("structure/",df_RMSD_all$models.y[q]))
      df_RMSD_all$RMSD[q]<-rmsd(a=pdb_1,b=pdb_2,fit=T)
      if(q%%10000==0){
          write.csv(df_RMSD_all,"df_RMSD_all.csv",row.names = F)
      }
    } 
    write.csv(df_RMSD_all,"df_RMSD_all.csv",row.names = F)
  }
}
