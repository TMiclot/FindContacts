# FindContacts

<!---  --->
## Description
FindContacts is a Tcl script that adds functionalities to the [VMD software](https://www.ks.uiuc.edu/Research/vmd/). It especially allows to analyze and visualize interactions in a single molecule or in a complex. Several functionalities have been added, including the creation and visualization of interactions network, neighborings residues search, 2D-RMSD map, etc. But the script also integrates improvements to existing features in VMD.

<!---  --->
## What can you do with FindContacts ?
* Make networks and statistic tables of _Hydrogen Bonds_ and _Salt Bridges_ and _Hydrophobic Interactions_ for protein and nucleic. But only HBonds interactions is available for other molecules.
* Compare interactions between multiple simulations.
* Draw 3D interaction network in VMD. Enhance interacting residues as 3D spheres in VMD.
* Calculate a contact map betwen tow selections (based on distance).
* Calculate the RMSD-2D map.
* Many other things :smile:

<!---  --->
## Runing the script into VMD
1. Load your MD simulation.
2. Open VMD TkConsole and source the script. You can do this easily by writing _source_ then _script name_: **source findcontacts.tcl**

Please consider that FindContacts work only on Linux or MacOS, not Windows.

<!---  --->
## Technical details
* The script is written in Tcl.
* Always use " " to write each argument.
* All output files are in .dat , it correspond to text file. So they can be read, plot, etc ... very easily.
* The script always use the TOP molecule, else for the FC-rmsd-molecules command.
* You can create scripts that use FindContacts commands. :heart: To do this, you only need to add the command line to source FindContacts at the beginning of your script. :thumbsup:
### Dependency
**cpptraj** in [AmberTools](https://ambermd.org/AmberTools.php). Optional: it is only necessary for only for clustering.

<!---  --->
## Known bugs
* Sometimes the script causes the VMD window to be unresponsive. Just wait for the script to finish working.
* Sometimes an error message written in red appears in the TkConsole. It is only that VMD bugs due to too much memory usage. To solve this problem, you just have to restart VMD and start the analysis from where it was.

<!--- --->
## Usage overview
|Command                   |Function|
FC-help                     |Display this help.
FC-receptor-ligand          |Interactions networks and tables of receptor/ligand, along the simulation.
FC-protein-changes          |Interactions networks and tables of protein changes, along the simulation.
FC-network-trace-bonds      |Draw the interaction network bonds.
FC-network-trace-resids     |Darw resid in interaction network as spheres.
FC-interactions-similitudes |Find the conserved interactions between multiples simulations, from network data.
FC-neighbors                |Give the frequence of residues in selection_1 within CutOff of selection_2.
FC-comcom                   |Measure the COM-COM distance of two selection, along the simulation.
FC-polrot                   |Gives the polar coordinates of selection 2 with respect to selection 1
FC-map-contacts             |Calculate the contact map between tow selections, along the simulation.
FC-sasa                     |Calculate SASA of a selection, along the simulation.
FC-sasa-rl                  |Calculate SASA and Interface for receptor/ligand, along the simulation.
FC-rmsf                     |Calculate the RMSF of a selection, along the simulation.
FC-rmsd                     |Calculate the RMSD of a selection, along the simulation.
FC-2drmsd                   |Calculate the 2D-RMSD of a selection.
FC-2drmsd-molecules         |Calculate the RMSD between tow molecules, along their simulations.
FC-rgyr                     |Calculate the radius of gyration of a selection, along the simulation.
FC-drgyr                    |Calculate the difference in radius of gyration of a selection.
FC-drgyr-molecules          |Calculate the DRGYR between tow molecules, along their simulations.
FC-clustering               |Make a clustering of the simulation. (CPPTRAJ)
FC-renumber                 |Renumer resid of a selection from a given number.
FC-reconstruct              |Center the selection and remove atom coordinate jumps.
