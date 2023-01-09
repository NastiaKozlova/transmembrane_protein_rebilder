part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
df_start<-df_start%>%mutate(test=F)
v_RMSD<-10
name<-1
for (name in 1:nrow(df_start)) {
  part<-paste0(parta,df_start$name[name],"/")
  setwd(part)
  if (!dir.exists("RMSD/")){dir.create("RMSD/")}
  print(paste0("start ",Sys.time()))
  if (file.exists("df_RMSD_all.csv")){
    df_start$test[name]<-T
#    df_RMSD<-read.csv(paste0(part,"df_RMSD_all.csv"),stringsAsFactors =  F)

  }
}
df_start<-df_start%>%filter(test)
part<-paste0(parta,df_start$name[1],"/")
df_RMSD<-read.csv(paste0(part,"df_RMSD_all.csv"),stringsAsFactors =  F)
if (nrow(df_start)>1){
  for (name in 2:nrow(df_start)) {
    part<-paste0(parta,df_start$name[name],"/")
    df_RMSD_add<-read.csv(paste0(part,"df_RMSD_all.csv"),stringsAsFactors =  F)
    df_RMSD<-rbind(df_RMSD,df_RMSD_add)
  }
}
p<-ggplot(data=df_RMSD, aes(x=RMSD))+
  labs(x="RMSD, A")+
  geom_density()+
  scale_x_continuous(breaks = seq(from=0,to=10,by=0.5),labels =  seq(from=0,to=10,by=0.5))+
  theme_bw()
ggsave(p,filename = paste0(part_start,"results/calibrtion_first_parts.png"), width = 20, height = 15, units = c("cm"), dpi = 300 ) 