# transmembrane_protein_rebilder

## Introduction

This project is focused on alternative folds searching if it is possible. It rebuilds the structure of TMD by splitting the structure of the protein monomer and with docking create alternative folds. The main feature of this project is the grouping and sorting algorithm which can identify complex structures with highest quality for shortest amount of time. Scripts are capable of combining up to 3 fragments of protein structure. In the future, it's possible to extend to combine more parts. 

## Algorithm

Algoritm contains following parts:
- Isolation of the two parts of protein structure with are will be combine
- Docking of the first two parts of the protein (a receptor is a structure with a longher sequence, a ligand with a shorter sequence) using PathDock
- Reassembly of the protein structure in accordance with the numbering of amino acids in the structure
- Sorting the resulting complexes based on RMSD with the original structure (structures are grouped according to the principle similar structures are equally similar to control)
- Compile lists of amino acids that interact between parts of a protein
- Calculation of interaction energy of protein parts
- Sorting possible similar pairs of structures based on lists of amino acids interacting between protein parts, where the percentage of interacting amino acids is more than 50%
- All protein pairs selected as "probably similar" were tested for similarity based on RMSD
- On the basis of RMSD, the structures were sorted by the convergence method, according to the principle that the more similar structures to the one under study, the more likely it is.
- Structures were considered likely if there were more than 5 similar structures and more than 1/1000 of the original structures obtained after docking and more than 95% of the quartile by the number of similar structures
- The orientation of the resulting structures was determined based on the angle between the membrane and the average vector of the terminal transmembrane domains of the structure under study and the average vector of the terminal transmembrane domains of the studied and initial structures
- The structure was considered the original (WT) if the orientation of the transmembrane domains was the same as in the original structure and the RMSD with the initial structure was less than 5\AA
- If there is an initial structure among the probable structures, then the probable structures are those whose interaction energy is less than or equal to the WT structure with the maximum energy
- If there is no WT structure among the probable structures, then structures with energies below the 25% quartile are considered probable, out of all probable structures

r_scripts/main_script.R - script which you should start after putting all files into they places.
r_scripts/pdb_prepare.R ",part_start - prepare parts of protein structures and run docking using PatchDock __Doesn't run on Debian 10__ You will see error "Segmentation error are appeared, problem with PathDock, program is docking shorter part to there longer one" 

## Dependencies

For this project, you need several programs
Download PathDock from https://bioinfo3d.cs.tau.ac.il/PatchDock/ and open the archive in the directory programs/PatchDock.

# Folder organisation to run file 

Open script r_scripts/main_script.R, change the path to your directory  and run the script
Create directory start and put there put structure start/start.pdb and their separation protein start/df_parts.csv
CHARMM36 force fields in to start/toppar
protein toplogy in start/df_topology.csv
if you neeed to rename domain names on the plots  put renamed plots in to start/domain_name.csv
File examples are present in start_example folder
/home/nastia/Dissertation/structure_rebilder/xylE/transmembrane_protein_rebilder/start/df_parts.csv


Output for the first part will the structure_prediction/'name of your parts'/fin_structure which are determined in yours 
The Final will be in the directory start/fin_prediction