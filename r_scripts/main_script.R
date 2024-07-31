part_start<-"path to transmembrane_protein_rebilder"
part_start<-paste0(getwd(),"/")
setwd(part_start)
#install.packages("dplyr")
#install.packages("bio3d")
#install.packages("readr")
#install.packages("ggplot2")
#install.packages("cowplot")

#install readr
#prepare pdb and run PatchDOCK first time
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/pdb_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_to_pdb_structure_no.R ",part_start),ignore.stdout=T,wait = T)
#energy calculation between protein parts
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_energy.R ",part_start),ignore.stdout=T,wait = T)
#first filtration based on difference of tested between (starting structure and analyzing structure)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/filter_pre_intercation_check.R ",part_start),ignore.stdout=T,wait = T)
#amino acids searching which interact between two investigation parts of protein 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_interactions.R ",part_start),ignore.stdout=T,wait = T)
#second sort of structures based on list of amino acids 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/sort_structure_interactions.R ",part_start),ignore.stdout=T,wait = T)
#calculation RMSD between presortes structures
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_analysis.R ",part_start),ignore.stdout=T,wait = T)
#combine all RMSD files into 1 large file
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_combine_files_RMSD_analysis.R ",part_start),ignore.stdout=T,wait = T)
#caluculate the most appropriate RMSD cut out to separate structures
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/calibration_first_parts.R ",part_start),ignore.stdout=T,wait = T)
#count statistics, combine structures and write results
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_no_7_fin_report.R ",part_start),ignore.stdout=T,wait = T)


#prepare first docking data output
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_data_output.R ",part_start),ignore.stdout=T,wait = T)


#calculate interprotein interaction using Ring2 and interactions between domains based on Ring2 data  
#will give an error, make manually
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)


#make interactions plot 
#need pictures

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_plot_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)


#TMD orinetation picture
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_TMD_orientation_picture.R ",part_start),ignore.stdout=T,wait = T)

#VMD run and make pictures

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_plot_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)


#prepare pdb and run PatchDOCK second time
#change 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/pdb_convert.R ",part_start),ignore.stdout=T,wait = T)
#convert structures after PatchDOCK run
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_to_pdb_structure_add.R ",part_start),ignore.stdout=T,wait = T)
#energy calculation between parts of the proteins (1+2 against 3)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_energy_add.R ",part_start),ignore.stdout=T,wait = T)
#first sort based on difference of RMSD between predicted structures and starting structure 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/filter_pre_intercation_check_add.R ",part_start),ignore.stdout=T,wait = T)
#calculate list of interprotein interactions
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/check_interactions_add.R ",part_start),ignore.stdout=T,wait = T)
#second sort based on interprotein interactions
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/sort_structure_interactions_add.R ",part_start),ignore.stdout=T,wait = T)
#Calculate RMSD between all probably similar pairs of proteins
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_analysis.R ",part_start),ignore.stdout=T,wait = T)
#merge RMSD files into one
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_combine_files_RMSD_analysis.R ",part_start),ignore.stdout=T,wait = T)
#Calibrate threshold of RMSD between different protein conformations
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/calibration_second_parts.R ",part_start),ignore.stdout=T,wait = T)
#select most propable structures
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/TMD_add_7_fin_report.R ",part_start),ignore.stdout=T,wait = T)
#output of all parts dockng
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_data_output.R ",part_start),ignore.stdout=T,wait = T)

#calculate interprotein interaction using Ring2 and interactions between domains based on Ring2 data
#will give an error, make manually
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_plot_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_TMD_orientation_picture.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_plot_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)

#count control (starting structure parameters)
#will give an error, make manually
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/control_TMD_interactions.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_part_compare_control_exprement_interactions.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_part_compare_control_exprement_interactions.R ",part_start),ignore.stdout=T,wait = T)
