# transmembrane_protein_rebilder

## Introduction

This project is focused on alternative folds searching if it is possible. It rebuilds the structure of TMD by splitting the structure of the protein monomer and with docking create alternative folds. The main feature of this project is the grouping and sorting algorithm which can identify complex structures with highest quality for shortest amount of time. Scripts are capable of combining up to 3 fragments of protein structure. In the future, it's possible to extend to combine more parts. 

## Algorithm

r_scripts/pdb_prepare.R ",part_start - prepare parts of protein strructures and run docking using PatchDock __Doesn't run on Debian 10__ Segmentation error are occure, problem with PathDock, programm is docking shourter part to there longer one 

r_scripts/convert_to_pdb_structure_no.R - arrange protein parts that parts of proteins are in the same order as in the starting structure and counting aminoacids which intercation beetween the parts and sort imblsible structures based on the length of gup between parts of proteins.

r_scripts/check_energy.R - counting energy interactions between parts of structure

r_scripts/filter_pre_intercation_check.R - 

r_scripts/check_interactions.R - 

sort_structure_interactions.R - 

TMD_no_7_analysis.R - counting RMSD between all selected as simulart parts of structues

TMD_no_7_combine_files_RMSD_analysis.R - combine all RMSD files in one

TMD_no_7_fin_report.R - 

first_part_data_output.R - create output file 

## Dependencies

For this project, you need several programs
Download PathDock from https://bioinfo3d.cs.tau.ac.il/PatchDock/ and open the archive in the directory programs/PatchDock.
Create directory start and put there put structure start/start.pdb and their separation protein start/df_parts.csv
Open script r_scripts/main_script.R, change the path to your directory  and run the script

Output for the first part will the structure_prediction/'name of your parts'/fin_structure which are determined in yours 
The Final will be in the directory start/fin_prediction
