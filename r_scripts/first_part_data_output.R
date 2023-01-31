#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)
library(ggplot2)
setwd(part_start)
df_plot_name<-read.csv("start/domain_name.csv",stringsAsFactors =  F)
prot_name<-strsplit(part_start,split = "/")[[1]]
prot_name<-prot_name[length(prot_name)-1]
i<-2
i<-1
sort_structures<-function(df_start,i){
  v_start<-length(list.files(paste0(df_start$name[i],"/structure")))
  df_start_all<-read.csv(paste0(df_start$name[i],"/fin.csv"),stringsAsFactors = F)
  df_start_all<-df_start_all%>%mutate(angle=angle*90/pi*2)
  df_start_all<-df_start_all%>%mutate(angle_mem=angle_mem*90/pi*2)
  
  df_start_all<-df_start_all%>%mutate(persent_align=round(align_models/group_models*100,digits = 1))
  df_start_all<-df_start_all%>%mutate(orientarion="between")
  
  df_start_all$orientarion[abs(df_start_all$angle)<60]<-"as WT"
  df_start_all$orientarion[abs(df_start_all$angle)>120]<-"inverted"
  df_start_all$orientarion[df_start_all$orientarion=="as WT"&df_start_all$RMSD<5]<-"WT"
  df_start_all$orientarion[abs(df_start_all$angle_mem)<30]<-"between"
  
  test<-min(c(df_start$first_part_finish[i]-df_start$first_part_start[i],df_start$second_part_finish[i]-df_start$second_part_start[i]))
  if(test>50){
    df_start_all<-df_start_all%>%filter(orientarion!="between")
  }
  df_start_all<-df_start_all%>%mutate(name=df_start$name[i])
  df_start_all<-left_join(df_start_all,df_start)
  df_start_all<-df_start_all%>%filter(group_models>5)
  df_start_all<-df_start_all%>%filter(group_models>v_start/1000)
  df_start_all<-df_start_all%>%filter(group_models>=quantile(df_start_all$group_models,0.95))
  if(length(df_start_all$orientarion[df_start_all$orientarion%in%"WT"])>0){
    v_energy_test<-max(df_start_all$bond_energy_fs[df_start_all$orientarion%in%"WT"])
    df_start_all<-df_start_all%>%filter(bond_energy_fs<=v_energy_test)
  }else{
    df_start_all<-df_start_all%>%filter(bond_energy_fs<=quantile(df_start_all$bond_energy_fs,probs = 0.25))
  }
  
  df_start_all<-df_start_all%>%mutate(frequence=group_models/sum(df_start_all$group_models)*100)
  df_start_all<-df_start_all%>%mutate(RMSD=round(RMSD,digits = 2))
  df_start_all<-df_start_all%>%mutate(angle=round(angle,digits = 0))
  return(df_start_all)
}
part<-paste0(part_start,"structure_prediction/")
#start<-read.pdb(paste0(part_start,"start/start.pdb"))
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

#check total quantity of patchdock compexes

df_start_all<-sort_structures(df_start,1)

if(nrow(df_start)>1){
  for (i in 2:nrow(df_start)) {
    df_start_add<-sort_structures(df_start,i)
    
    df_start_all<-rbind(df_start_all,df_start_add)
  }
}


df_start_all<-df_start_all%>%mutate(RMSD=round(RMSD,digits = 2))
df_start_all<-df_start_all%>%mutate(frequence=round(frequence,digits = 1))
df_start_all<-df_start_all%>%mutate(angle=round(angle,digits = 0))
df_start_all<-df_start_all%>%mutate(persent_align=round(align_models/group_models*100,digits = 1))

df_start_all<-df_start_all%>%select(name,orientarion,RMSD,frequence, RMSD,persent_align,group_models, angle,
                                    first_part_model,first_part_start,first_part_finish,second_part_model,second_part_start, 
                                    second_part_finish,group_number,bond_energy_fs)
if(dir.exists(paste0(part_start,"results/"))) {system(command = paste0("rm -r ",part_start,"results/"),ignore.stdout=T,wait = T)}
if(!dir.exists(paste0(part_start,"results/"))){dir.create(paste0(part_start,"results/"))}
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
df_start_all<-df_start_all%>%mutate(first_second=paste0(first_part_model,"-",second_part_model))
df_start_all<-df_start_all%>%mutate(first_part_group_number=as.character(group_number))
df_start_all$angle<-NULL
df_angle<-data.frame(matrix(ncol=2,nrow=4))
colnames(df_angle)<-c("orientation","angle")
df_angle$orientation<- c("between","WT","as WT","inverted")
df_angle$angle<- c(90,0,0,180)
df_start_all<-left_join(df_start_all,df_angle,by=c("orientarion"="orientation"))

df_start_all<-df_start_all%>%mutate(plot_name=paste(group_number))
df_start_all<-df_start_all%>%mutate(first_part_center=(first_part_start+first_part_finish)/2)
df_start_all<-df_start_all%>%mutate(second_part_center=(second_part_start+second_part_finish)/2)
v_min<-min(df_start_all$first_part_start,df_start_all$second_part_start)
v_max<-max(df_start_all$first_part_finish,df_start_all$second_part_finish)
df_start_all<-left_join(x = df_start_all,y = df_plot_name,by=c("first_part_model"="part_name"))
df_start_all<-left_join(x = df_start_all,y = df_plot_name,by=c("second_part_model"="part_name"))

df_start_all<-left_join(x = df_start_all,y = df_plot_name,by=c("first_part_model"="part_name"))
df_start_all<-df_start_all%>%mutate(first_plot_name=fin_name)
df_start_all$fin_name<-NULL
df_start_all<-left_join(x = df_start_all,y = df_plot_name,by=c("second_part_model"="part_name"))
df_start_all<-df_start_all%>%mutate(second_plot_name=fin_name)
df_start_all$fin_name<-NULL
df_start_all<-df_start_all%>%mutate(name=paste0(first_plot_name,"-",second_plot_name))
df_start_all<-df_start_all%>%mutate(statistic=paste0("frequence ",frequence,"\nRMSD ",RMSD,"\n",persent_align))
p<-ggplot(data=df_start_all)+
  labs(x="number of aminoaids",y="")+
  geom_text(aes(x=v_min,y=plot_name,label=statistic,angle=0))+
  geom_text(aes(x=first_part_center,y=plot_name,label=first_plot_name,angle=0,color="1"))+
  geom_text(aes(x=second_part_center,y=plot_name,label=second_plot_name,angle=angle,color="2"))+
  scale_y_discrete(breaks = NULL,labels = NULL)+
  scale_x_continuous(breaks = NULL,labels = NULL,limits = c((v_min-20),v_max))+
  facet_grid(name~.,scales = "free", space = "free")+theme_bw()
ggsave(p,filename = paste0(part_start,"results/",prot_name,"_first_part_rebilder.png"), width = 20, height = 20, units = c("cm"), dpi = 200 ) 
