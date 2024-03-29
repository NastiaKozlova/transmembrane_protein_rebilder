part_start = commandArgs(trailingOnly=TRUE)
Sys.time()
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_start)
parta<-paste0(part_start,"structure_prediction/")
df_start<-read.csv("start/df_parts.csv",stringsAsFactors = F)
name<-1
i<-1
for (name in 1:nrow(df_start)) {
  
  part<-paste0(parta,df_start$name[name],"/")
  setwd(part)  
  if (!dir.exists("compare_interaction/")){dir.create("compare_interaction/")}
  df_RMSD<-read.csv("df_RMSD.csv",stringsAsFactors = F)
  if(!dir.exists(paste0("Energy/"))){dir.create(paste0("Energy/"))}
  if(!dir.exists(paste0("Energy/structure_psf"))){dir.create(paste0("Energy/structure_psf"))}
  if(!dir.exists(paste0("Energy/tcl"))){dir.create(paste0("Energy/tcl"))}
  if(!dir.exists(paste0("Energy/namd_tcl"))){dir.create(paste0("Energy/namd_tcl"))}
  if(!dir.exists(paste0("Energy/full"))){dir.create(paste0("Energy/full"))}
  if(!dir.exists(paste0("Energy/first"))){dir.create(paste0("Energy/first"))}
  if(!dir.exists(paste0("Energy/second"))){dir.create(paste0("Energy/second"))}
  if(!dir.exists(paste0("Energy/bond_energy"))){dir.create(paste0("Energy/bond_energy"))}
  #  if(!dir.exists(paste0(name[part],"/MD/stabilisation"))){dir.create(paste0(name[part],"/MD/stabilisation"))}
  #  if(!dir.exists(paste0(name[part],"/MD/stabilisation/protein"))){dir.create(paste0(name[part],"/MD/stabilisation/protein"))}
  #  if(!dir.exists(paste0(name[part],"/MD/stabilisation/pdb"))){dir.create(paste0(name[part],"/MD/stabilisation/pdb"))}
  #  if(!dir.exists(paste0(name[part],"/MD/stabilisation/quench"))){dir.create(paste0(name[part],"/MD/stabilisation/quench"))}
  #  if(!dir.exists(paste0(name[part],"/MD/stabilisation/dcd"))){dir.create(paste0(name[part],"/MD/stabilisation/dcd"))}
  for (i in 1:nrow(df_RMSD)) {
    df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
    df_psfgen[1,1]<-paste0('cd ',part,'\n',
                           'mol delete all\n',  
                           'package require psfgen \n',  
                           'resetpsf\n',
                           'topology ',part_start,'start/toppar/top_all36_prot.rtf\n',
                           'topology ',part_start,'start/toppar/toppar_water_ions_namd.str\n',  
                           'pdbalias residue HIS HSE\n',
                           'pdbalias atom ILE CD1 CD\n',
                           'segment U { pdb  structure/',df_RMSD$models[i],'\n',
                           '}\n',
                           '\ncoordpdb structure/',df_RMSD$models[i],' U\n',
                           'regenerate angles dihedrals\n',  
                           'guesscoord\n',  
                           'writepdb Energy/structure_psf/',df_RMSD$models[i],'.pdb\n',  
                           'writepsf Energy/structure_psf/',df_RMSD$models[i],'.psf\n\n\n\nexit now')
    write.table(df_psfgen,paste0('Energy/tcl/psfgen_',df_RMSD$models[i],'.tcl'),col.names = F,row.names = F,quote = F)
    system(command = paste0("vmd -dispdev text -e ",part,'Energy/tcl/psfgen_',df_RMSD$models[i],'.tcl'),ignore.stdout=T,wait = T) 
  }
  df_RMSD<-df_RMSD%>%mutate(first_energy=NA)
  df_RMSD<-df_RMSD%>%mutate(second_energy=NA)
  df_RMSD<-df_RMSD%>%mutate(total_energy=NA)
  for (i in 1:nrow(df_RMSD)) {
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('cd ',part,'\n\npackage require namdenergy')
    df_tcl[1,2]<-paste0('mol new {Energy/structure_psf/',df_RMSD$models[i],'.psf} type {psf}')
    df_tcl[1,3]<-paste0('mol addfile {Energy/structure_psf/',df_RMSD$models[i],'.pdb} type {pdb}') 
    df_tcl[1,4]<-paste0('set sel2 [atomselect top "protein"]')
    df_tcl[1,5]<-paste0('set sel1 [atomselect top "protein and resid>=',df_start$first_part_start[name],' and resid<=',df_start$first_part_finish[name],'"]')
    df_tcl[1,6]<-paste0('set sel3 [atomselect top "protein and resid>=',df_start$second_part_start[name],' and resid<=',df_start$second_part_finish[name],'"]')
    df_tcl[1,7]<-paste0('namdenergy -sel $sel2  -bond -angl -dihe -impr -conf -vdw -elec -nonb -all -cutoff 12 -skip 0 -ofile Energy/full/',df_RMSD$models[i],
                        ' -switch 10 -exe ',part_start,'programs/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_start,'start/toppar/par_all36_carb.prm -par ',
                        part_start,'start/toppar/par_all36_cgenff.prm -par ',part_start,'start/toppar/par_all36_lipid.prm -par ',part_start,
                        'start/toppar/par_all36m_prot.prm -par ',part_start,'start/toppar/par_all36_na.prm -par ',
                        part_start,'start/toppar/par_all36_prot.prm -par ',part_start,'start/toppar/toppar_water_ions_namd.str')
    df_tcl[1,8]<-paste0('namdenergy -sel $sel1  -bond -angl -dihe -impr -conf -vdw -elec -nonb -all -cutoff 12 -skip 0 -ofile Energy/first/',df_RMSD$models[i],
                        ' -switch 10 -exe ',part_start,'programs/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_start,'start/toppar/par_all36_carb.prm -par ',
                        part_start,'start/toppar/par_all36_cgenff.prm -par ',part_start,'start/toppar/par_all36_lipid.prm -par ',part_start,
                        'start/toppar/par_all36m_prot.prm -par ',part_start,'start/toppar/par_all36_na.prm -par ',
                        part_start,'start/toppar/par_all36_prot.prm -par ',part_start,'start/toppar/toppar_water_ions_namd.str')
    df_tcl[1,9]<-paste0('namdenergy -sel $sel3  -bond -angl -dihe -impr -conf -vdw -elec -nonb -all -cutoff 12 -skip 0 -ofile Energy/second/',df_RMSD$models[i],
                        ' -switch 10 -exe ',part_start,'programs/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_start,'start/toppar/par_all36_carb.prm -par ',
                        part_start,'start/toppar/par_all36_cgenff.prm -par ',part_start,'start/toppar/par_all36_lipid.prm -par ',part_start,
                        'start/toppar/par_all36m_prot.prm -par ',part_start,'start/toppar/par_all36_na.prm -par ',
                        part_start,'start/toppar/par_all36_prot.prm -par ',part_start,'start/toppar/toppar_water_ions_namd.str')
    df_tcl[1,10]<-paste0('namdenergy -sel $sel3 $sel1  -vdw -elec -nonb -cutoff 12 -skip 0 -ofile Energy/bond_energy/',df_RMSD$models[i],
                        ' -switch 10 -exe ',part_start,'programs/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_start,'start/toppar/par_all36_carb.prm -par ',
                        part_start,'start/toppar/par_all36_cgenff.prm -par ',part_start,'start/toppar/par_all36_lipid.prm -par ',part_start,
                        'start/toppar/par_all36m_prot.prm -par ',part_start,'start/toppar/par_all36_na.prm -par ',
                        part_start,'start/toppar/par_all36_prot.prm -par ',part_start,'start/toppar/toppar_water_ions_namd.str')
    df_tcl[1,11]<-'mol delete all'
    df_tcl[1,12]<-'\n\nexit now'
    write.table(df_tcl,file =paste0('Energy/namd_tcl/',df_RMSD$models[i],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
    print(paste0('vmd -dispdev text -e ',part,'Energy/namd_tcl/',df_RMSD$models[i],'.tcl'))
    
    system(command = paste0('vmd -dispdev text -e ',part,'Energy/namd_tcl/',df_RMSD$models[i],'.tcl'),ignore.stdout=T,wait = T) 
    if(file.exists(paste0("Energy/first/",df_RMSD$models[i]))){
      df_first<-read.table(paste0("Energy/first/",df_RMSD$models[i]), sep="", header=T, na.strings ="", stringsAsFactors= F)
      df_second<-read.table(paste0("Energy/second/",df_RMSD$models[i]), sep="", header=T, na.strings ="", stringsAsFactors= F)
      df_full<-read.table(paste0("Energy/full/",df_RMSD$models[i]), sep="", header=T, na.strings ="", stringsAsFactors= F)
      df_bond_energy<-read.table(paste0("Energy/bond_energy/",df_RMSD$models[i]), sep="", header=T, na.strings ="", stringsAsFactors= F)
      df_RMSD$first_energy[i]<-df_first$Total[1]
      df_RMSD$second_energy[i]<-df_second$Total[1]
      df_RMSD$total_energy[i]<-df_full$Total[1]
      df_RMSD$bond_energy[i]<-df_bond_energy$Total[1]
    }
  }
  df_RMSD<-df_RMSD%>%mutate(bond_energy_fs=total_energy-first_energy-second_energy)
  write.csv(df_RMSD,"df_complex_energy.csv",row.names = F)
  system(command = paste0("rm -r ",part,"Energy"),ignore.stdout=T,wait = T)
}
