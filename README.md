# transmembrane_protein_rebilder

## Introduction

This project are focused on alternative folds searching if it is possible. It rebuild structure of TMD by spliting structure of protein monomer and with docking create alternative folds. Main feature of this project is groupping and sorting algoritm which can identefy complex structures with higest quality for shourtest amount of time. Scripts are capable of combaning up to 3 fragments of protein structure. In the future its possible to extend to combine more parts. 

## Algoritm

r_scripts/pdb_prepare.R ",part_start - prepare parts of protein strructures and run docking using PatchDock __Doesn't run on Debian 10__ Segmentation error, problem with PathDock, programm is docking shourter part to there longer one 
r_scripts/convert_to_pdb_structure_no.R - arrange protein parts that parts of proteins are in the same order as in the starting structure and counting aminoacids which intercation beetween the parts and sort imblsible structures based on the length of gup between parts of proteins.
r_scripts/check_energy.R - counting energy interactions between parts of structure
r_scripts/filter_pre_intercation_check.R - 
r_scripts/check_interactions.R - 
sort_structure_interactions.R - 
TMD_no_7_analysis.R - RMD ana
TMD_no_7_combine_files_RMSD_analysis.R
TMD_no_7_fin_report.R
first_part_data_output.R

## Dependenses

For this project you need several programs
Download PathDock from https://bioinfo3d.cs.tau.ac.il/PatchDock/ and open arcave in the directory programs/PatchDock.
Create directory start and put there put structure start/start.pdb and there separation protein start/df_parts.csv
Open script r_scripts/main_script.R, change path to your dicrectiry  and run scipt

Output for first part will the structure_prediction/'name of your parts'/fin_structure which are determine in yours 
Final will be in the directory start/fin_prediction
