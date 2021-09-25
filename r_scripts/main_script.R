part_start<-"path to transmembrane_protein_rebilder"
#install.packages("dplyr")
#install.packages("bio3d")
#install.packages("readr")
#install.packages("ggplot2")
#install.packages("cowplot")

#install readr
#prepare pdb and run PatchDOCK first time
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/pdb_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_to_pdb_structure_no.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_energy.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/filter_pre_intercation_check.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_interactions.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/sort_structure_interactions.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_analysis.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_combine_files_RMSD_analysis.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_fin_report.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_data_output.R ",part_start),ignore.stdout=T,wait = T)
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/remove_extra_files.R ",part_start),ignore.stdout=T,wait = T)

#prepare pdb and run PatchDOCK second time
#change 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/pdb_convert.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_to_pdb_structure_add.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_energy_add.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_interactions_add.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/sort_structure_interactions.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_analysis_top_a.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_combine_files_RMSD_analysis.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_fin_report.R ",part_start),ignore.stdout=T,wait = T)