library(dplyr)
library(ggplot2)
library(bio3d)
part_start = commandArgs(trailingOnly=TRUE)
part<-paste0(part,"structure_prediction/")
setwd(part)
v_system<-c("pachdock")
df_system<-data.frame(matrix(ncol=2,nrow = length(v_system)))
colnames(df_system)<-c("system","total_models")
df_system<-df_system%>%mutate(system=v_system)
if (!dir.exists("din")) {dir.create("din")}
i<-1
j<-1

for (i in 1:nrow(df_system)) {
  df_system$total_models[i]<-length(list.files(paste0(v_system[i],"/topSolutions")))
  v_groups<-list.files(paste0(df_system$system[i],"/din"))
  df_groups<-data.frame(matrix(ncol=8,nrow = length(v_groups)))
  colnames(df_groups)<-c("system","group_number","group_name","group_models","align_models","min_RMSD","max_RMSD","best_model")
  df_groups$system<-df_system$system[i]
  df_groups$group_number<-c(1:nrow(df_groups))
  df_groups$group_name<-v_groups
  for (j in 1:nrow(df_groups)) {
    df_RMSD<-read.csv(paste0(df_system$system[i],"/din/",df_groups$group_name[j]),stringsAsFactors = F)
    df_groups$best_model[j]<-df_RMSD$models.x[1]
    df_groups$min_RMSD[j]<-min(df_RMSD$RMSD_control)
    df_groups$max_RMSD[j]<-max(df_RMSD$RMSD_control)
    df_groups$group_models[j]<-nrow(df_RMSD)
    df_RMSD<-df_RMSD%>%filter(RMSD_control<5)
    df_groups$align_models[j]<-nrow(df_RMSD)
  }
  df_groups<-left_join(df_groups,df_system,by="system")
  write.csv(df_groups,file = paste0("din/",df_system$system[i],".csv"),row.names = F)                  
}
i<-1
df_fin<-read.csv(file = paste0("din/",df_system$system[i],".csv"),stringsAsFactors =  F) 
if (nrow(df_system)>1) {
  for (i in 2:nrow(df_system)) {
    df_groups<-read.csv(file = paste0("din/",df_system$system[i],".csv"),stringsAsFactors =  F) 
    df_fin<-rbind(df_fin,df_groups)
    
  }
}
write.csv(df_fin,"fin_analysis.csv",row.names = F)
if (!dir.exists("fin_structure")){dir.create("fin_structure")}
i<-1
i<-2
for (i in 1:nrow(df_fin)) {
  pdb<-read.pdb(paste0(df_fin$system[i],"/str/",df_fin$group_number[i],"/",df_fin$best_model[i]))
  write.pdb(pdb,paste0("fin_structure/",df_fin$system[i],"_",df_fin$group_number[i],".pdb"))
}