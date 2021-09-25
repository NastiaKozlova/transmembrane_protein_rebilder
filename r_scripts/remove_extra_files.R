#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
v_RMSD<-10
persent_rigth<-3
for (name in 1:nrow(df_start)) {
  parta<-paste0(part_start,"structure_prediction/",df_start$name[name],"/")
  system(command = paste0("rm -r  ",parta,"compare_interaction/ "),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r  ",parta,"compare_TEMP/ "),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r  ",parta,"interactions/ "),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r  ",parta,"RMSD/ "),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r  ",parta,"RMSD_TEMP/ "),ignore.stdout=T,wait = T)
  system(command = paste0("rm -r  ",parta,"RMSD_test/ "),ignore.stdout=T,wait = T)
}
