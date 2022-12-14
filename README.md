# FindContacts

<!---  --->
## Description
FindContacts is a Tcl script that adds functionalities to the [VMD software](https://www.ks.uiuc.edu/Research/vmd/). It especially allows to analyze and visualize interactions in a single molecule or in a complex. Several functionalities have been added, including the creation and visualization of interactions network, neighborings residues search, 2D-RMSD map, etc. But the script also integrates improvements to existing features in VMD.

<img src="network.png" width="30%" height="30%">


<!---  --->
## What can you do with FindContacts ?
* Make networks and statistic tables of _Hydrogen Bonds_ and _Salt Bridges_ and _Hydrophobic Interactions_ for protein and nucleic. But only HBonds interactions is available for other molecules.
* Compare interactions between multiple simulations.
* Draw 3D interaction network in VMD. Enhance interacting residues as 3D spheres in VMD.
* Calculate a contact map betwen tow selections (based on distance).
* Calculate the RMSD-2D map.
* Many other things :smile:
<!--- --->


<!---  --->
## Runing the script into VMD
1. Load your MD simulation.
2. Open VMD TkConsole and source the script. You can do this easily by writing _source_ then _script name_: **source findcontacts.tcl**

Please consider that FindContacts work only on Linux or MacOS, not Windows.
<!--- --->


<!---  --->
## Technical details
* The script is written in Tcl.
* Always use " " to write each argument.
* All output files are column based files in _.dat_ format. It correspond to text file. So they can be read, plot, etc ... very easily.
* The script always use the TOP molecule, else for the FC-rmsd-molecules command.
* Find-Contacts use equivalent selection variables as [GetContacts](https://getcontacts.github.io/).
* You can create scripts that use FindContacts commands. :heart: To do this, you only need to add the command line to source FindContacts at the beginning of your script. :thumbsup:
* You can use a graphing utility to display an heatmap. For exemple, in Gnuplot you can wrote: _plot 'file' using 1:2:3 with image_

### Dependency
**CPPTRAJ** in [AmberTools](https://ambermd.org/AmberTools.php). Optional: it is only necessary for only for clustering.

### Known bugs
* Sometimes the script causes the VMD window to be unresponsive. Just wait for the script to finish working.
* Sometimes an error message written in red appears in the TkConsole. It is only that VMD bugs due to too much memory usage. To solve this problem, you just have to restart VMD and start the analysis from where it was.
<!--- --->


<!--- --->
## Usage overview
|Command                     |Function|
|----------------------------|--------|
|FC-help                     |Display this help.|
|FC-receptor-ligand          |Interactions networks and tables of receptor/ligand, along the simulation.|
|FC-protein-changes          |Interactions networks and tables of protein changes, along the simulation.|
|FC-network-trace-bonds      |Draw the interaction network bonds.|
|FC-network-trace-resids     |Darw resid in interaction network as spheres.|
|FC-interactions-similitudes |Find the conserved interactions between multiples simulations, from network data.|
|FC-neighbors                |Give the frequence of residues in _selection 1_ within CutOff of _selection 2_.|
|FC-comcom                   |Measure the COM-COM distance of two selections, along the simulation.|
|FC-polrot                   |Gives the polar coordinates of _selection 2_ with respect to _selection 1_|
|FC-map-contacts             |Calculate the contact map between tow selections, along the simulation.|
|FC-sasa                     |Calculate SASA of a selection, along the simulation.|
|FC-sasa-rl                  |Calculate SASA and Interface for receptor/ligand, along the simulation.|
|FC-rmsf                     |Calculate the RMSF of a selection, along the simulation.|
|FC-rmsd                     |Calculate the RMSD of a selection, along the simulation.|
|FC-2drmsd                   |Calculate the 2D-RMSD of a selection.|
|FC-2drmsd-molecules         |Calculate the RMSD between tow molecules, along their simulations.|
|FC-rgyr                     |Calculate the radius of gyration of a selection, along the simulation.|
|FC-drgyr                    |Calculate the difference in radius of gyration of a selection.|
|FC-drgyr-molecules          |Calculate the DRGYR between tow molecules, along their simulations.|
|FC-clustering               |Make a clustering of the simulation. **Requires CPPTRAJ**|
|FC-renumber                 |Renumer resid of a selection from a given number.|
|FC-reconstruct              |Center the selection and remove atom coordinate jumps.|

### New selection variables added
|Variable|Equivalent VMD command|
|--------|----------------------|
|protein_SB_Pos |((resname HIS HSD HSE HSP HIE HIP HID and name ND1 NE2) or (resname LYS and name NZ) or (resname ARG and name NH1 NH2))|
|protein_SB_Neg |((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))|
|protein_SB_All |((resname HIS HSD HSE HSP HIE HIP HID and name ND1 NE2) or (resname LYS and name NZ) or (resname ARG and name NH1 NH2) or (resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))|
|protein_Hyd |(hydrophobic and not backbone and type C C1 C2 CA CB CC CE CI CK CQ CR CT CW)|
|nucleic_SB  |(name OP1 OP2)|

To know all the selection variables used natively by VMD, go to its [documentation](https://www.ks.uiuc.edu/Research/vmd/current/ug/node90.html).
<!--- --->


<!--- --->
## Using the commands
<!--- ............................................................. --->
### FC-receptor-ligand
This command compute interaction search and analysis of H-Bonds, Salt bridges and Hydrophobic interactions between Receptor and Ligand, along the simulation. At the end you will get files to draw the interaction-type network, and others corresponding to interaction-type statistics for the _min_frequence_. 
The _min_frequence_ parameter select residue pairs with a frequence >= _min_frequence_. The value must be between 0 and 1, and can take three decimals.
#### Usage
FC-receptor-ligand "OutputFileName" "receptor" "ligand" "min_frequence"
#### Outfiles
* The outfiles are:
  - OutputFileName-HBonds-Network-Frequences.dat
  - OutputFileName-HBonds-TableStats.dat
  - OutputFileName-Hydrophobic-Network-Frequences.dat
  - OutputFileName-Hydrophobic-TableStats.dat
  - OutputFileName-SaltBridges-Network-Frequences.dat
  - OutputFileName-SaltBridges-TableStats.dat
* "Network" files are used to draw network bonds and residues, but also to get interactions similitudes between multiple simulations. Type TableStats is a table with statistical analysis of interactions.
* In the TableStats files, interaction are selected by the frequence (>= minimum frequency) and by the mean of the atomic distance (<= 4 angstroms).
#### Warnings
* It may happen that the statistical table files are empty. This means that the frequency of interactions are all lower than the minimum frequency chosen, or the mean of distance is > 4 angstroms. In this case you should reduce the value of the minimum frequency.
* For protein-X interaction search, please always use "receptor" for your (protein/peptide/amino acids) selection.
* Interactions available:

|System          |Type of interactions available|
|----------------|------------------------------|
|protein-protein |H-Bonds, Salt bridges and Hydrophobic|
|protein-nucleic |H-Bonds and Salt bridges available|
|nucleic-nucleic |only H-Bonds|
|others          |only H-Bonds|


<!--- ............................................................. --->
### FC-protein-changes
This command compute interaction search and analysis of H-Bonds, Salt bridges and Hydrophobic interactions of protein changes along the simulation. At the end you will get files to draw the interaction-type network, and others corresponding to interaction-type statistics for the _min_frequence_. 
The _min_frequence_ parameter select residue pairs with a frequence >= _min_frequence_. The value must be between 0 and 1, and take three decimals.
#### Usage
FC-protein-changes "OutputFileName" "min_frequence"
#### Outfiles
* The outfiles are: 
  - OutputFileName-HBonds-Network-Frequences.dat
  - OutputFileName-HBonds-TableStats.dat
  - OutputFileName-Hydrophobic-Network-Frequences.dat
  - OutputFileName-Hydrophobic-TableStats.dat
  - OutputFileName-SaltBridges-Network-Frequences.dat
  - OutputFileName-SaltBridges-TableStats.dat
* "Network" files are used to draw network bonds and residues, but also to get interactions similitudes between multiple simulations.
* "TableStats" files are tables with statistical analysis of interactions. Interactions are selected by the frequence (>= minimum frequency) and by the mean of the atomic distance (<= 4 angstroms).
#### Warnings
* It may happen that the statistical table files are empty. This means that the frequency of interactions are all lower than the minimum frequency chosen, or the mean of distance is > 4 angstroms. In this case you should reduce the value of the minimum frequency.


<!--- ............................................................. --->
### FC-network-trace-bonds
This command draw the network from the _Input File_, at the curent frame.
It take 4 parameters :

| |Parameter|Meaning|
|-|---------|-------|
|1| InputFileName |Name of the file to be read. (Ex: xxx-HBonds-TableStat-Strong.dat)|
|2| expansion     |Used to grow, or reduce, radius of network cylinder.|
|3| MinRadius     |Select minimum radius you want to display. It take three decimals.|
|4| color         |Color of the networ cylinder. All VMD color name are are available. (red , green , blue ,...).|

#### Usage
FC-network-trace-bonds "InputFileName" "MinRadius" "expansion" "color"
#### Technical note
You can clear all the drawed object with the command : _draw delete all_


<!--- ............................................................. --->
### FC-network-trace-resids

This command draw residues involved in the network from the (Input File) as spheres, at the curent frame.
It take 5 parameters :

| |Parameter|Meaning|
|-|---------|-------|
|1|InputFileName |Name of the file to be read. (Ex: xxx-HBonds-TableStat-Strong.dat)|
|2|expansion     |Used to grow, or reduce, radius of network cylinder.|
|3|MinRadius     |Select minimum radius you want to display. Ittake three decimals.|
|4|color         |Color of the networ cylinder. All VMD color name are are available. (red , green , blue ,...).|
|5|1/2           |Write 1 or 2. In the InputFile residues are written in _resid in column 1_ and _resid in column 2_ format.|

#### Usage
FC-network-trace-resids "InputFileName" "MinRadius" "expansion" "color" "1/2"
#### Technical note
You can clear all the drawed object with the command : _draw delete all_


<!--- ............................................................. --->
### FC-comcom
This command measure the COM-COM distance evolution along the simulation.
#### Usage
FC-rmsd "OutFileName" "selection_1" "selection_2"
#### Outfiles
* The outfile will be _OutFileName-COMCOM.dat_.So if you write _MyFile_ as file name, the file will be _MyFile-COMCOM.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame and COM-COM-distance


<!--- ............................................................. --->
### FC-polrot
This command measure the polar coordinates of _selection 2_ with respect to _selection 1_. This the polar
#### Usage
FC-polrot "OutFileName" "selection_1" "selection_2"
#### Outfiles
* The outfile is _OutFileName-POLROT.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-POLROT.dat_
* The outfile is column based file. Columns are respectively: Frame, R, Theta angle (polar), Psi angle (polar), Theta angle (degree), Psi angle (degree).


<!--- ............................................................. --->
### FC-sasa
This command calculate the SASA of the selection. You can aslo specify a restriction parameter, but it's not necessary. It's similar to the sasa.tcl script available on the VMD website. << The restrict option can be used to prevent internal protein voids or pockets from affecting the surface area results. >>
#### Usage
FC-sasa "OutFileName" "selection" "restriction"
#### Outfiles
* The outfile is _OutFileName-SASA.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-SASA.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame and SASA


<!--- ............................................................. --->
### FC-sasa-rl
This command calculate the SASA for receptor and ligand interaction.
It compute _SASA Receptor_ , _SASA Ligand_ , _SASA Receptor/Ligand complex_ and _Interface of Receptor/Ligand_
#### Usage
FC-sasa-rl "OutFileName" "receptor" "ligand"
#### Outfiles
* The outfiles are:
  - OutFileName-Receptor-SASA.dat
  - OutFileName-Ligand-SASA.dat
  - OutFileName-Complex-SASA.dat
  - OutFileName-RLInterface.dat
* The outfile is column based file. Column 1 and 2 are respectively Frame and SASA


<!--- ............................................................. --->
### FC-rmsf
This command calculate the RMSF for each atom in the selection. The alignment is automaticaly done before RMSF calculations.
### Esemples of possible selection
* Protein selection : _protein and name CA_
* Nucleic selection : _nucleic and name N1_
####
FC-rmsf "OutFileName" "selection"
#### Outfiles
* The outfile will be _OutFileName-RMSF.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-RMSF.dat_
* The outfile is column based file. Column 1 and 2 are respectively AtomNumber and RMSF


<!--- ............................................................. --->
### FC-rmsd
This command calculate the RMSD of the selection. The alignment is automaticaly done during RMSD calculations.
#### Usage
FC-rmsd "OutFileName" "selection"
#### Outfiles
* The outfile will be _OutFileName-RMSD.dat_ .So if you write _MyFile_ as file name, the file will be _MyFile-RMSD.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame and RMSD


<!--- ............................................................. --->
### FC-2drmsd
This command calculate the 2D-RMSD of the selection. The alignment is automaticaly done during 2D-RMSD calculations.
You can use a graphing utility like Gnuplot to display an heatmap.
For exemple, in Gnuplot you can wrote:     plot 'file' using 1:2:3 with image
#### Usage
FC-2drmsd "OutFileName" "selection" "stride"
#### Outfiles
* The outfile will be _OutFileName-2dRMSD.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-2dRMSD.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame_A and Frame_B. Column 3 is the RMSD value.


<!--- ............................................................. --->
### FC-rmsd-molecules
This command calculate the RMSD between the selection into tow trajectories. It's an 2D-RMSD like process.
The alignment of the tow molecule is automaticaly done during the calculations
For exemple, you can plot an heatmap or make a histogram distribution.
#### Usage
FC-rmsd-molecules "OutFileName" "Molecule_1_ID" "Molecule_2_ID" "selection" "stride"
#### Outfiles
* The outfile will be _OutFileName-RMSDmolecules.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-RMSDmolecules.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame_A and Frame_B. Column 3 is the RMSD value.


<!--- ............................................................. --->
### FC-rgyr
This command calculate the radius of gyration of the selection.
#### Usage
FC-rgyr "OutFileName" "selection"
#### Outfiles
* The outfile will be _OutFileName-RMSD.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-RMSD.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame and RGYR.


<!--- ............................................................. --->
### FC-dgyr
This command calculate the difference in radius of gyration of the selection.
#### Usage
FC-drgyr "OutFileName" "selection" "stride"
#### Outfiles
* The outfile will be _OutFileName-2dRMSD.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-DRGYR.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame_A and Frame_B. Column 3 is the DRGYR value.


<!--- ............................................................. --->
### FC-drgyr-molecules
This command calculate the DRGYR between the selection into tow trajectories. It's an 2D-RGYR like process.
For exemple,y ou can plot an heatmap or make a histogram distribution.
#### Usage
FC-drgyr-molecules "OutFileName" "Molecule_1_ID" "Molecule_2_ID" "selection" "stride"
#### Outfiles
* The outfile will be _OutFileName-DRGYRmolecules.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-DRGYRmolecules.dat_
* The outfile is column based file. Column 1 and 2 are respectively Frame_A and Frame_B. Column 3 is the DRGYR value.


<!--- ............................................................. --->
### FC-map-contacts
This command measure the distance between atoms in _selection 1_ and _selection2_ for all frames in the MD.
Here, you can use a graphing utility to display an heatmap, or make an animated heatmap along the MD.
#### Usage
FC-map-contacts "OutFileName" "selection 1" "selection 2"
#### Esemple of possible selections

|Type of selection  |VMD command|
|-------------------|-----------|
|Protein selection  |protein and name CA|
|                   |resid 26 to 251 and name CA|
|Nucleic selection  |nucleic and name N1|
|                   |resid 1 to 24 and name N1|
|Ligand  selection  |resname LIG|
|                   |resid 124|

#### Outfiles
* The outfile will be _OutFileName-MapContacts.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-MapContacts.dat_
* The outfile is column based file. Column 1 and 2 are Resid (for protein or nucleic) or Atom Number (for ligand). Column >= 3 are distance value for one frame. Ex: column 3 is all distance for frame 0, column 4 for frame 1, etc.
#### Warnings
* It always use all defined atoms in the selection. So, for ligand molecule you will get all atoms. In that case, atoms are numbered according to their order in the selection, but not according to their index.
* For protein selection, always use : _name CA_
* For nucleic selection, always use : _name N1_ or _name C1'_


<!--- ............................................................. --->
### FC-clustering
This command use CPPTRAJ to calculate a cluster of structures in the simulation. So, it can take a long time.
The command automatically recognizes the parameter file and the trajectory files loaded in VMD. 
Ensure you have load Amber, or AmberTools, in the terminal were you run VMD. <br>
It take 6 drawing parameters :
| |Parameter|Meaning|
|-|---------|-------|
|1|FolderName|name of the CPTRAJ work dir. The final name will be FolderName-cluster|
|2|mask      |mask selection, in CPPTRAJ format. See this [page](https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/)|
|3|start     |frame to begin reading.|
|4|stop      |frame to stop reading at.|
|5|stride    |stride for reading in trajectory frames. Stride is applied to all trajectory files loaded (the first and the others).|
|6|first/all |"first" applies _start_value stop_value stride_value_ to the first trajectory file loaded, and (1 last stride_value) to the others. <br> "all" applies _start_value stop_value stride_value_ to the all trajectory files loaded.|

#### Usage
FC-clustering "FolderName" "mask" "start" "stop" "stride" "first/all"
#### Outfiles
* The outfiles are the file to run CPPTRAJ, and all output file of this software.
#### Warnings
* Please take into account that VMD has a command _measure cluster_, but which is less efficient than the analysis by CPPTRAJ.


<!--- ............................................................. --->
### FC-neighbors
This command give the frequence of residues in _selection_1_ within _CutOff_ of _selection_2_.
#### Usage
FC-neighbors "OutFileName" "selection_1" "selection_2" "CutOff"
#### Outfiles
* The outfile will be _OutFileName-NEIGHBORS.dat_. So if you write _MyFile_ as file name, the file will be _MyFile-NEIGHBORS.dat_
* The outfile is column based file. Column 1 and 2 are respectively Residue , Count in the MD, Frequence


<!--- ............................................................. --->
### FC-renumber
This command renumber the IDs of the residues of the _selection_, starting from the _startIDnumber_.
#### Usage
FC-renumber "selection" "startIDnumber"
#### Warnings
* It is possible that the Sequence Viewer in VMD does not show the correct numbers after using this command.


<!--- ............................................................. --->
### FC-reconstruct
This command uses a wrap/unwrap protocol to cleanly remove large jumps in atom coordinates of the selection (and the errors this causes in the bonds), and put the selection in the center of the cell.
You can chose the parallelepiped or orthorhombic option to << wrap the atoms into the unitcell parallelepiped or the corresponding orthorhombic box with the same volume and center as the _non-orthrhombic_ unitcell. >>
#### Usage
FC-reconstruct "selection" "parallelepiped | orthorhombic"
#### Warnings
* If you don't know which option to choose, you can use parallelepiped as defult.
* Using FC-reconstruct command may cause deformation of the water box.
