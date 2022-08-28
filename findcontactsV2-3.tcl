###########################  INFORMATIONS  ##########################
# NAME        : FindContacts
# AUTHOR      : Tom MICLOT
# DATE        : 2022 AUG 27
# VERSION     : 2.3
# STATUT      : Stable
# LANGUAGE    : Tcl
# LICENCE     : MIT License
# DESCRIPTION : FindContacts is a script that adds functionality to
#               the VMD software. It especially allows to analyze and
#               visualize interactions in a protein or in a
#               complex. Several functionalities have been added,
#               including the creation and visualization of interaction
#               network, neighboring residue search, 2D-RMSD map, etc.
#               But the script also integrates improvements to existing
#               features in VMD, for example the possibility to export
#               directly the RMSF in a file.
# RUNNING     : Open VMD TkConsole and write    source findcontacts.tcl
#               It work only on Linux or MacOS, not Windows.
# DEPENDENCY  : cpptraj (in AmberTools, only for clustering)
#####################################################################
#####################################################################













#=====================================================================
#=====================================================================
puts "
███████╗██╗███╗   ██╗██████╗        ██████╗ ██████╗ ███╗   ██╗████████╗ █████╗  ██████╗████████╗███████╗
██╔════╝██║████╗  ██║██╔══██╗      ██╔════╝██╔═══██╗████╗  ██║╚══██╔══╝██╔══██╗██╔════╝╚══██╔══╝██╔════╝
█████╗  ██║██╔██╗ ██║██║  ██║█████╗██║     ██║   ██║██╔██╗ ██║   ██║   ███████║██║        ██║   ███████╗
██╔══╝  ██║██║╚██╗██║██║  ██║╚════╝██║     ██║   ██║██║╚██╗██║   ██║   ██╔══██║██║        ██║   ╚════██║
██║     ██║██║ ╚████║██████╔╝      ╚██████╗╚██████╔╝██║ ╚████║   ██║   ██║  ██║╚██████╗   ██║   ███████║
╚═╝     ╚═╝╚═╝  ╚═══╝╚═════╝        ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═╝   ╚══════╝

v2.3 \t by Tom MICLOT                                                             
"
puts "\t Info) Find-Contacts is loading ..."
#=====================================================================
#=====================================================================















#=====================================================================
#=====================================================================
proc ::Set_Global_Variables {} {
    #
    puts "\t Info) Find-Contacts -- Initialization of variables : START."
    #
    global __TERMINAL__ __FRAME_INITIAL__ __FRAME_TOTAL__ __PROTEIN_SALT_BRIDGE_POSITIVES__ __PROTEIN_SALT_BRIDGE_NEGATIVES__ __PROTEIN_SALT_BRIDGE_ALL__ __PROTEIN_HYDROPHOBIC__ __NUCLEIC_SALT_BRIDGE__ protein_SB_Pos protein_SB_Neg protein_SB_All protein_Hyd nucleic_SB
    #
    set __TERMINAL__ [exec ls /usr/bin/ | grep terminal | head -1]
    #
    set __FRAME_INITIAL__ 0
    set __FRAME_TOTAL__ [molinfo top get numframes]
    #
    #
    #--------------------------------------------------------------------- default variables used in the script
    #
    set __PROTEIN_SALT_BRIDGE_POSITIVES__ "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) or ((resname LYS) and (name NZ)) or ((resname ARG) and (name NH1 NH2))"
    set __PROTEIN_SALT_BRIDGE_NEGATIVES__ "((resname ASP) and (name OD1 OD2)) or ((resname GLU) and (name OE1 OE2))"
    set __PROTEIN_SALT_BRIDGE_ALL__ "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) or ((resname LYS) and (name NZ)) or ((resname ARG) and (name NH1 NH2)) or ((resname ASP) and (name OD1 OD2)) or ((resname GLU) and (name OE1 OE2))"
    set __PROTEIN_HYDROPHOBIC__ "hydrophobic and not backbone and (type C C1 C2 CA CB CC CE CI CK CQ CR CT CW)"
    set __NUCLEIC_SALT_BRIDGE__ "(name OP1 OP2)" 
    #
    #--------------------------------------------------------------------- default variables used in the script -- converted for easy user usage
    #
    set protein_SB_Pos $__PROTEIN_SALT_BRIDGE_POSITIVES__
    set protein_SB_Neg $__PROTEIN_SALT_BRIDGE_NEGATIVES__
    set protein_SB_All $__PROTEIN_SALT_BRIDGE_ALL__
    set protein_Hyd $__PROTEIN_HYDROPHOBIC__
    set nucleic_SB $__NUCLEIC_SALT_BRIDGE__
    #
    #
    puts "\t Info) Find-Contacts -- Initialization of variables : END."
}
# end proc
#
# execute the procedure. This process is needed to get all global variable into other PROC
::Set_Global_Variables
#
#=====================================================================
#=====================================================================














#===================================================================== FONCTIONNEL
#===================================================================== 
proc ::sasa {__OUTPUT_FILE_NAME__ __MOLECULE__ __RESTRICTION__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-SASA.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame SASA"
    # update variables from global
    upvar 2 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 2 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __COUNTER__ 0
    #
    #-----------------------------------------------------------------
    #
    set __SEL_MOLECULE__ [atomselect top "$__MOLECULE__"]
    set __SEL_RESTRICTION__ [atomselect top "$__RESTRICTION__"]
    #
    #-----------------------------------------------------------------
    #
    puts "\t SASA -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # upadte atomselects to the the correct frame
        $__SEL_MOLECULE__ frame $__FRAME__
        $__SEL_RESTRICTION__ frame $__FRAME__
        # compute sasa
        set __SASA__ [measure sasa 1.4 $__SEL_MOLECULE__ -restrict $__SEL_RESTRICTION__]
        #
        puts $__OUTPUT_FILE__ "$__FRAME__ $__SASA__"
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\t SASA -- $__PERCENT__ %"
    }
    #
    #-----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "\t SASA -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc ::align {__SEL__} {
  upvar 2 __FRAME_TOTAL__ __FRAME_TOTAL__
  upvar 2 __FRAME_INITIAL__ __FRAME_INITIAL__
  #
  set __ALL_ATOMS__ [atomselect top all]
  set __REFERENCE__ [atomselect top $__SEL__ frame 0]
  set __SELECTION__ [atomselect top $__SEL__]
  #
  set __COUNTER__ 0
  #
  #------------------------------------------------------------------
  #
  puts "Align -- START."
  #
  for { set __FRAME__ $__FRAME_INITIAL__ } { $__FRAME__ < $__FRAME_TOTAL__ } { incr __FRAME__ } {
    $__SELECTION__ frame $__FRAME__
    $__ALL_ATOMS__ frame $__FRAME__
    $__ALL_ATOMS__ move [measure fit $__SELECTION__ $__REFERENCE__]
    # print the curent state
    set __COUNTER__ [expr $__COUNTER__ +1.0]
    set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
    puts "\tAlign -- $__PERCENT__ %"
  }
  # end for
  puts "Align -- END."
}
#end proc
#=====================================================================
#=====================================================================













#=====================================================================
#=====================================================================
proc ::get-list-atoms {__LIST_ATOMS_1__ __LIST_ATOMS_2__} {
    #----------------------------------------------------------------
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    #
    set __LIST_PAIR_ATOMS__ {}
    #
    set __ITEM__ {}
    set __COUNT__ {}
    set __COUNTERS__ {}
    #----------------------------------------------------------------
    if {[llength $__LIST_ATOMS_1__] == [llength $__LIST_ATOMS_2__]} {
        for {set j 0} {$j < [llength $__LIST_ATOMS_1__]} {incr j} {
        set __ATOM_1__ [lindex $__LIST_ATOMS_1__ $j]
        set __ATOM_2__ [lindex $__LIST_ATOMS_2__ $j]
        set __PAIR_ATOMS__ "$__ATOM_1__-$__ATOM_2__"
        set __LIST_PAIR_ATOMS__ [concat $__LIST_PAIR_ATOMS__ $__PAIR_ATOMS__]
        }
        # end for
    }
    # end if
    return $__LIST_PAIR_ATOMS__
}
# end proc
#=====================================================================
#=====================================================================














#=====================================================================
#=====================================================================
proc ::make-table {__OUTPUT_FILE_NAME__ __LIST_PAIR_ATOMS_ALLDYNAMIC__ __MIN_FREQ__} {
    #----------------------------------------------------------------
    # set outfiles
    set __OUTPUT_FILE_NAME_1__ "$__OUTPUT_FILE_NAME__-TableStat-All.dat"
    set __OUTPUT_FILE_1__ [open "$__OUTPUT_FILE_NAME_1__" w]
    puts $__OUTPUT_FILE_1__ "ResName1:ResID1:AtomName1:Atom1Idx ResName2:ResID2:AtomName2:Atom2Idx Mean StDev Min Max Count Frequence"
    #
    set __OUTPUT_FILE_NAME_2__ "$__OUTPUT_FILE_NAME__-TableStat-Strong.dat"
    set __OUTPUT_FILE_2__ [open "$__OUTPUT_FILE_NAME_2__" w]
    puts $__OUTPUT_FILE_2__ "ResName1:ResID1:AtomName1:Atom1Idx ResName2:ResID2:AtomName2:Atom2Idx Mean StDev Min Max Count Frequence"
    #----------------------------------------------------------------
    # 
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    # set variables
    set __ITEM__ {}
    set __COUNT__ {}
    set __COUNTERS__ {}
    set __ITEM_SUP_MIN_FREQ_LIST__ {}
    #
    #---------------------------------------------------------------- 
    #
    foreach __ITEM__ $__LIST_PAIR_ATOMS_ALLDYNAMIC__ {
         dict incr __COUNTERS__ $__ITEM__
    }
    # end foreach
    #
    dict for {__ITEM__ __COUNT__} $__COUNTERS__ {
        # calculate Frequence and round up to 4 decimal
        set __FREQUENCE__ [format "%.4f" [expr $__COUNT__ / ($__FRAME_TOTAL__ +0.0)]]
        #
        lassign [split $__ITEM__ - ] __ATOM_1__ __ATOM_2__
        #
        set __ATOM_1_SELECTION__ [atomselect top "index $__ATOM_1__"]
        set __ATOM_2_SELECTION__ [atomselect top "index $__ATOM_2__"]
        #
        # making resname
        set __ATOM_1_NAME__ [$__ATOM_1_SELECTION__ get name]
        set __ATOM_2_NAME__ [$__ATOM_2_SELECTION__ get name]
        #
        # making resid
        set __ATOM_1_RESID__ [$__ATOM_1_SELECTION__ get resid]
        set __ATOM_2_RESID__ [$__ATOM_2_SELECTION__ get resid]
        #
        # making resname
        set __ATOM_1_RESNAME__ [$__ATOM_1_SELECTION__ get resname]
        set __ATOM_2_RESNAME__ [$__ATOM_2_SELECTION__ get resname]
        #
        # making Stats af the distance between Atom1 and Atom2
        set __DISTANCES__ [measure bond "$__ATOM_1__  $__ATOM_2__" frame all]
        set __MEAN__ [format "%.4f" [vecmean $__DISTANCES__]]
        set __STDEV__ [format "%.4f" [vecstddev $__DISTANCES__]]
        set __MIN__ [format "%.4f" [lindex [lsort -real $__DISTANCES__] 0]] 
        set __MAX__ [format "%.4f" [lindex [lsort -real $__DISTANCES__] end]]
        #
        # write information in outputfile_1
        puts $__OUTPUT_FILE_1__ "$__ATOM_1_RESNAME__:$__ATOM_1_RESID__:$__ATOM_1_NAME__:$__ATOM_1__ $__ATOM_2_RESNAME__:$__ATOM_2_RESID__:$__ATOM_2_NAME__:$__ATOM_2__ $__MEAN__ $__STDEV__ $__MIN__ $__MAX__ $__COUNT__ $__FREQUENCE__"
        #
        # write information in outputfile_2 only if frequence >= min_freq && mean <= 4 angstrom
        if {$__FREQUENCE__ >= $__MIN_FREQ__ && $__MEAN__ <= 4.0} {
            puts $__OUTPUT_FILE_2__ "$__ATOM_1_RESNAME__:$__ATOM_1_RESID__:$__ATOM_1_NAME__:$__ATOM_1__ $__ATOM_2_RESNAME__:$__ATOM_2_RESID__:$__ATOM_2_NAME__:$__ATOM_2__ $__MEAN__ $__STDEV__ $__MIN__ $__MAX__ $__COUNT__ $__FREQUENCE__"
	}
	#end if
        #
        foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    }
    # end dict for
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE_1__
    close $__OUTPUT_FILE_2__
}
#end proc
#=====================================================================
#=====================================================================











#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-neighbors {__OUTPUT_FILE_NAME__ __SELECT_1__ __SELECT_2__ __CUTOFF__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-NEIGHBORS.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "resname.resid count freq"
    #
    set __SELECTION_1__ [atomselect top "$__SELECT_1__"]
    set __SELECTION_2__ [atomselect top "$__SELECT_2__"]
    #
    set __LIST_NAME_ID__ {}
    set __LIST_NAME_ID_ALLTRAJ__ {}
    #
    set __ITEM__ {}
    set __COUNT__ {}
    set __COUNTERS__ {}
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-neighbors -- START."
    #
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # upadte atomselects to the the correct frame
        foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
        set __SELECTION__ [atomselect top "same residue as $__SELECT_1__ within $__CUTOFF__ of $__SELECT_2__" frame $__FRAME__]
        # compute the fitting
        set __RES_NAME__ [$__SELECTION__ get resname]
        set __RES_ID__ [$__SELECTION__ get resid]
        #
        if {[llength $__RES_NAME__] == [llength $__RES_ID__]} {
            for {set j 0} {$j < [llength $__RES_NAME__]} {incr j} {
                set __NAME__ [lindex $__RES_NAME__ $j]
                set __ID__ [lindex $__RES_ID__ $j]
                set __NAME_ID__ "$__NAME__:$__ID__"
                set __LIST_NAME_ID__ [lsort -unique [concat $__LIST_NAME_ID__ $__NAME_ID__ ]]
            }
            # end for
        }
        #end if
        set __LIST_NAME_ID_ALLTRAJ__ [concat $__LIST_NAME_ID_ALLTRAJ__ $__LIST_NAME_ID__ ]
        set __LIST_NAME_ID__ {}
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\tFC-neighbors -- traj analysis $__PERCENT__ %"
    }
    #end for 
    #
    #----------------------------------------------------------------
    #
    puts "\tFC-neighbors -- residue list analysis start"
    #
    foreach __ITEM__ $__LIST_NAME_ID_ALLTRAJ__ {
         dict incr __COUNTERS__ $__ITEM__
    }
    # end foreach
    #
    dict for {__ITEM__ __COUNT__} $__COUNTERS__ {
        # calculate Frequence and round up to 4 decimal
        set __FREQUENCE__ [format "%.4f" [expr $__COUNT__ / ($__FRAME_TOTAL__ +0.0)]]
        # write information in outputfile
        puts $__OUTPUT_FILE__ "$__ITEM__ $__COUNT__ $__FREQUENCE__"
        #
	#end i
    }
    # end dict for
    puts "\tFC-neighbors -- residue list analysis end"
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
#    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-neighbors -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
# 
proc FC-receptor-ligand {__OUTPUT_FILE__ __RECEPTOR__ __LIGAND__ __MIN_FREQ__} {
    #
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    upvar 1 __PROTEIN_SALT_BRIDGE_POSITIVES__ __PROTEIN_SALT_BRIDGE_POSITIVES__
    upvar 1 __PROTEIN_SALT_BRIDGE_NEGATIVES__ __PROTEIN_SALT_BRIDGE_NEGATIVES__
    upvar 1 __PROTEIN_HYDROPHOBIC__ __PROTEIN_HYDROPHOBIC__
    upvar 1 __NUCLEIC_SALT_BRIDGE__ __NUCLEIC_SALT_BRIDGE__
    #
    set __COUNTER__ $__FRAME_INITIAL__
    #
    set __HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ {}
    set __SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ {}
    set __HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ {}
    #
    #-----------------------------------------------------------------
    #
    puts "FC-receptor-ligand -- START."
    #
    set __PROTEIN_NUCLEIC_CODE__ [list 1010 0101 1001 0110]
    set __RECEPTOR_LIGAND_TYPE_CODE__ [join [lindex [[atomselect top "$__RECEPTOR__"] get protein] 0][lindex [[atomselect top "$__RECEPTOR__"] get nucleic] 0][lindex [[atomselect top "$__LIGAND__"] get protein] 0][lindex [[atomselect top "$__LIGAND__"] get nucleic] 0]]
    #
    # check if receptor and ligand are protein or nucleic type
    #     if the code is 1010 receptor and ligand are protein type
    #     if the code is 0101 receptor and ligand are nucleic type
    #     if the code is 1001 receptor is protein type and ligand is nucleic type
    #     if the code is 0110 receptor is nucleic type and ligand is protein type
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        #
        # ........... display progress ...........
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\n\t Interaction search $__PERCENT__ %"
        # ........... HBONDS ...........
        #
        puts "\t\tH-Bonds search"
        #
        # HBonds search 1/2 : RECEPTOR is considered the donor and LIGAND is considered the acceptor.
        lassign [measure hbonds 3.5 70 [atomselect top "$__RECEPTOR__" frame $__FRAME__] [atomselect top "$__LIGAND__" frame $__FRAME__]] __RECEPTOR_ATOMS_COVAL_BOND_TO_H__ __RECEPTOR_ATOMS_H__ __RECEPTOR_ATOMS_H_ACCEPTORS__
        # HBonds search 1/2 : RECEPTOR is considered the donor and LIGAND is considered the acceptor.
        lassign [measure hbonds 3.5 70 [atomselect top "$__LIGAND__" frame $__FRAME__] [atomselect top "$__RECEPTOR__" frame $__FRAME__]] __LIGAND_ATOMS_COVAL_BOND_TO_H__ __LIGAND_ATOMS_H__ __LIGAND_ATOMS_H_ACCEPTORS__
        #
        # make one list with Atom_H and one list with Atom_H_Acceptors of receptor and ligands
        set __ATOMS_H__ [concat $__RECEPTOR_ATOMS_H__ $__LIGAND_ATOMS_H__]
        set __ATOMS_H_ACCEPTORS__ [concat $__RECEPTOR_ATOMS_H_ACCEPTORS__ $__LIGAND_ATOMS_H_ACCEPTORS__]
        #
        set __HBONDS_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__ATOMS_H__ $__ATOMS_H_ACCEPTORS__]
        set __HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ $__HBONDS_LIST_PAIR_ATOMS_FRAME__]
        #
        # ........... SALT BRIDGES ...........
        # Salt Bridges interaction search only if receptor and ligand ar not both nucleotide type
        if { [expr {"$__RECEPTOR_LIGAND_TYPE_CODE__" ne "0101"}] } {
            puts "\t\tSalt Bridges search"
            #
            lassign [measure contacts 4.0 [atomselect top "$__RECEPTOR__ and  $__PROTEIN_SALT_BRIDGE_POSITIVES__" frame $__FRAME__] [atomselect top "$__LIGAND__ and ( $__PROTEIN_SALT_BRIDGE_NEGATIVES__ or $__NUCLEIC_SALT_BRIDGE__ )" frame $__FRAME__]] __SALT_BRIDGE_LIST_ATOMS_1__ __SALT_BRIDGE_LIST_ATOMS_2__
            #
            set __SALTBRIDGES_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__SALT_BRIDGE_LIST_ATOMS_1__ $__SALT_BRIDGE_LIST_ATOMS_2__]
            set __SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ $__SALTBRIDGES_LIST_PAIR_ATOMS_FRAME__]
        }
        # end if
        #
        # ........... HYDROPHOBIC ...........
        # hydrophobic analysis only if receptor and ligand are both protein type
        if { [expr {"$__RECEPTOR_LIGAND_TYPE_CODE__" eq "1010"}] } {
            puts "\t\tHydrophobic search"
            lassign [measure contacts 3.9 [atomselect top "( $__RECEPTOR__ or $__LIGAND__ ) and $__PROTEIN_HYDROPHOBIC__" frame $__FRAME__]] __HYDROPHOBIC_LIST_ATOMS_1__ __HYDROPHOBIC_LIST_ATOMS_2__
            #
            set __HYDROPHOBICS_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__HYDROPHOBIC_LIST_ATOMS_1__ $__HYDROPHOBIC_LIST_ATOMS_2__]
            set __HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ $__HYDROPHOBICS_LIST_PAIR_ATOMS_FRAME__]
        }
        # end if
        #
        # ........... clear memory...........
        foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    }
    #end for

        #--------------------------------------------------------------------------
        puts "\n\t Interaction data processing."
        # Hbond processing
        #
        puts "\t\t HBonds analysis"
        ::make-table "$__OUTPUT_FILE__-HBonds" "$__HBONDS_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
        # Salt Bridge processing only if receptor and ligand are not both nucleic type
        #
        if { [expr {"$__RECEPTOR_LIGAND_TYPE_CODE__" ne "0101"}] } {
            puts "\t\t Salt Bridges analysis"
            #
            ::make-table "$__OUTPUT_FILE__-SaltBridges" "$__SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
        }
        #
        # hydrophobic processing only if receptor and ligand are both protein type
        if { [expr {"$__RECEPTOR_LIGAND_TYPE_CODE__" eq "1010"}] } {
            puts "\t\t Hydrophobic analysis"
            #
            ::make-table "$__OUTPUT_FILE__-Hydrophobics" "$__HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
        }
        #end if
    #
    #-----------------------------------------------------------------
    #
    puts "FC-receptor-ligand -- START."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#===================================================================== 
proc FC-protein-changes {__OUTPUT_FILE__ __MIN_FREQ__} {
    #
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    upvar 1 __PROTEIN_SALT_BRIDGE_POSITIVES__ __PROTEIN_SALT_BRIDGE_POSITIVES__
    upvar 1 __PROTEIN_SALT_BRIDGE_NEGATIVES__ __PROTEIN_SALT_BRIDGE_NEGATIVES__
    upvar 1 __PROTEIN_HYDROPHOBIC__ __PROTEIN_HYDROPHOBIC__
    #
    set __COUNTER__ $__FRAME_INITIAL__
    #
    set __HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ {}
    set __SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ {}
    set __HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ {}
    #
    #-----------------------------------------------------------------
    # 
    puts "FC-protein-changes -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        #
        # ........... display progress ...........
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\n\t Interaction search $__PERCENT__ %"
        # ........... HBONDS ...........
        #
        puts "\t\t1/3 H-Bonds search"
        #
        lassign [measure hbonds 3.5 70 [atomselect top "protein" frame $__FRAME__]] __ATOMS_COVAL_BOND_TO_H__ __ATOMS_H__ __ATOMS_H_ACCEPTORS__
        #
        set __HBONDS_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__ATOMS_H__ $__ATOMS_H_ACCEPTORS__]
        set __HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__HBONDS_LIST_PAIR_ATOMS_ALLFRAME__ $__HBONDS_LIST_PAIR_ATOMS_FRAME__]
        #
        # ........... SALT BRIDGES ...........
        # Salt Bridges interaction search only if receptor and ligand ar not both nucleotide type
        puts "\t\t2/3 Salt Bridges search"
        #
        lassign [measure contacts 4.0 [atomselect top "$__PROTEIN_SALT_BRIDGE_POSITIVES__" frame $__FRAME__] [atomselect top "$__PROTEIN_SALT_BRIDGE_NEGATIVES__" frame $__FRAME__]] __SALT_BRIDGE_LIST_ATOMS_1__ __SALT_BRIDGE_LIST_ATOMS_2__
        #
        set __SALTBRIDGES_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__SALT_BRIDGE_LIST_ATOMS_1__ $__SALT_BRIDGE_LIST_ATOMS_2__]
        set __SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__ $__SALTBRIDGES_LIST_PAIR_ATOMS_FRAME__]
        #
        # ........... HYDROPHOBIC ...........
        # hydrophobic analysis only if receptor and ligand are both protein type
        puts "\t\t3/3 Hydrophobic search"
        lassign lassign [measure contacts 3.9 [atomselect top "$__PROTEIN_HYDROPHOBIC__" frame $__FRAME__]] __HYDROPHOBIC_LIST_ATOMS_1__ __HYDROPHOBIC_LIST_ATOMS_2__
        #
        set __HYDROPHOBICS_LIST_PAIR_ATOMS_FRAME__ [::get-list-atoms $__HYDROPHOBIC_LIST_ATOMS_1__ $__HYDROPHOBIC_LIST_ATOMS_2__]
        set __HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ [concat $__HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__ $__HYDROPHOBICS_LIST_PAIR_ATOMS_FRAME__]
        #
        # ........... clear memory...........
        foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    }
    #end for

        #--------------------------------------------------------------------------
        puts "\n\t Interaction data processing."
        # Hbond processing
        puts "\t\t 1/3 HBonds analysis"
        ::make-table "$__OUTPUT_FILE__-HBonds" "$__HBONDS_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
        #
        # Salt Bridge processing only if receptor and ligand are not both nucleic type
        puts "\t\t 2/3 Salt Bridges analysis"
        ::make-table "$__OUTPUT_FILE__-SaltBridges" "$__SALTBRIDGES_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
        #
        # hydrophobic processing only if receptor and ligand are both protein type
        puts "\t\t 3/3 Hydrophobic analysis"
        #
        ::make-table "$__OUTPUT_FILE__-Hydrophobics" "$__HYDROPHOBICS_LIST_PAIR_ATOMS_ALLFRAME__" "$__MIN_FREQ__"
    #
    #
    #-----------------------------------------------------------------
    #
    puts "FC-protein-changes -- END."
}
#end proc
#=====================================================================
#=====================================================================










#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-renumber {__SELECT__ __MIN_NUMBER__} {
    #-----------------------------------------------------------------
    #
    puts "FC-renumber -- START"
    set __ID__ $__MIN_NUMBER__
    set __LIST_RESID_SELECT__ [lsort -unique -integer [[atomselect top "$__SELECT__"] get resid]]
    #
    foreach __ITEM__ $__LIST_RESID_SELECT__ {
        [atomselect top "resid $__ITEM__" frame all] set resid $__ID__
        set __ID__ [expr $__ID__ +1]
    }
    #-----------------------------------------------------------------
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-renumber -- END"
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-network-trace-bonds {__INPUT_FILE_NAME__ __FREQ_MIN__ __EXPANSION__ __COLOR__} {
    #-----------------------------------------------------------------
    # set color environment
    draw color $__COLOR__
    #
    set __MAX_LINE__ [exec wc -l "$__INPUT_FILE_NAME__" | cut -d " " -f1 ]
    #
    set p "p"
    #
    set __LIST_PAIR_RESIDS__ {}
    #
    #-----------------------------------------------------------------
    #
    puts "FC-network-trace-bonds -- START."
    puts "\t Search for residue pairs and their frequence in the file. Please wait ..."
    for {set __LINE_NUMBER__ 2} {$__LINE_NUMBER__ <= $__MAX_LINE__} {incr __LINE_NUMBER__} {
        # read file line by line
        set __LINE__ [exec sed -n "$__LINE_NUMBER__$p" $__INPUT_FILE_NAME__]
        #
        lassign [split $__LINE__ ] __RESNAME_RESID_1__ __RESNAME_RESID_2__
        lassign [split $__RESNAME_RESID_1__ :] __RESNAME_1__ __RESID_1__
        lassign [split $__RESNAME_RESID_2__ :] __RESNAME_2__ __RESID_2__
        # __FREQUENCE__ is the maximun frequence find for Resid1-Resid2 interaction
        set __FREQUENCE__ [exec cat $__INPUT_FILE_NAME__ | grep :$__RESID_1__: | grep :$__RESID_2__: | cut -d " " -f8 | sort | tail -1]
        #
        if {!([string first "$__RESID_1__-$__RESID_2__" $__LIST_PAIR_RESIDS__] != -1 || [string first "$__RESID_2__-$__RESID_1__" $__LIST_PAIR_RESIDS__] != -1)} {
            set __PAIR_RESID__ "$__RESID_1__-$__RESID_2__"
            set __LIST_PAIR_RESIDS__ [concat $__LIST_PAIR_RESIDS__ $__PAIR_RESID__]
            #
            if {$__FREQUENCE__ >= $__FREQ_MIN__} {
                # adding 0.00001 to ensure the radius is not 0
                set __RADIUS__ [expr ($__FREQUENCE__ * $__EXPANSION__) + 0.00001]
                # draw a cylinder between the tow geometric center of Resid_1 and Resid_2
                draw cylinder [measure center [atomselect top "resid $__RESID_1__"]] [measure center [atomselect top "resid $__RESID_2__"]] radius $__RADIUS__
                puts "\tDraw $__RESID_1__-$__RESID_2__ connection with radius = [format "%.4f" $__RADIUS__]"
            }
            #end if
        }
        #end if
    }
    #end for
    puts "FC-network-trace-bonds -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-network-trace-resid {__INPUT_FILE_NAME__ __FREQ_MIN__ __EXPANSION__ __COLOR__ __LOCATION__} {
    #-----------------------------------------------------------------
    # set color environment
    draw color $__COLOR__
    #
    set __MAX_LINE__ [exec wc -l "$__INPUT_FILE_NAME__" | cut -d " " -f1 ]
    #
    set p "p"
    #
    set __LIST_PAIR_RESIDS__ {}
    #
    #-----------------------------------------------------------------
    #
    puts "FC-network-trace-resid -- START."
    puts "\t Search for residue pairs and their weights in the file. Please wait ..."
    for {set __LINE_NUMBER__ 2} {$__LINE_NUMBER__ <= $__MAX_LINE__} {incr __LINE_NUMBER__} {
        # read file line by line
        set __LINE__ [exec sed -n "$__LINE_NUMBER__$p" $__INPUT_FILE_NAME__]
        #
        lassign [split $__LINE__ ] __RESNAME_RESID_1__ __RESNAME_RESID_2__
        lassign [split $__RESNAME_RESID_1__ :] __RESNAME_1__ __RESID_1
        lassign [split $__RESNAME_RESID_2__ :] __RESNAME_2__ __RESID_2
        # __FREQUENCE__ is the maximun frequence find for Resid1-Resid2 interaction
        set __FREQUENCE__ [exec cat $__INPUT_FILE_NAME__ | grep :$__RESID_1: | grep :$__RESID_2: | cut -d " " -f8 | sort | tail -1]
        #
        if {!([string first "$__RESID_1-$__RESID_2" $__LIST_PAIR_RESIDS__] != -1 || [string first "$__RESID_2-$__RESID_1" $__LIST_PAIR_RESIDS__] != -1)} {
            set __PAIR_RESID__ "$__RESID_1-$__RESID_2"
            set __LIST_PAIR_RESIDS__ [concat $__LIST_PAIR_RESIDS__ $__PAIR_RESID__]
            #
            if {$__FREQUENCE__ >= $__FREQ_MIN__} {
                # adding 0.00001 to ensure the radius is not 0
                set __RADIUS__ [expr ($__FREQUENCE__ * $__EXPANSION__) + 0.00001]
                # draw a cylinder between the tow geometric center of Resid_1 and Resid_2
                draw sphere [measure center [atomselect top "resid [set __RESID_$__LOCATION__]"]] radius $__RADIUS__
            puts "\tDraw [set __RESID_$__LOCATION__] sphere with radius = [format "%.4f" $__RADIUS__]"
            }
            #end if
        }
        #end if
        #
    }
    #end for
    #
    #-----------------------------------------------------------------
    #
    puts "FC-network-trace-resid -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-comcom {__OUTPUT_FILE_NAME__ __SELECT_1__ __SELECT_2__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-COMCOM.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame COM-COM-distance"
    #
    set __SELECTION_1__ [atomselect top "$__SELECT_1__"]
    set __SELECTION_2__ [atomselect top "$__SELECT_2__"]
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-comcom -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # upadte atomselects to the the correct frame
        $__SELECTION_1__ frame $__FRAME__
        $__SELECTION_2__ frame $__FRAME__
        # compute the fitting
        set __COM_SELECTION_1__ [measure center "$__SELECTION_1__" weight mass]
        set __COM_SELECTION_2__ [measure center "$__SELECTION_2__" weight mass]
        #
        set __COMCOM_DISTANCE__ [veclength [vecsub $__COM_SELECTION_1__ $__COM_SELECTION_2__]]
        # write data in the file
        puts $__OUTPUT_FILE__ "$__FRAME__ $__COMCOM_DISTANCE__"
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "FC-comcom -- $__PERCENT__ %"
    } 
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-comcom -- END."
}
#end proc
#=====================================================================
#=====================================================================












#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-polrot {__OUTPUT_FILE_NAME__ __SELECT_1__ __SELECT_2__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-POLROT.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame r theta psi theta(deg) psi(deg)"
    #
    set __SELECTION_1__ [atomselect top "$__SELECT_1__"]
    set __SELECTION_2__ [atomselect top "$__SELECT_2__"]
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set PI 3.1415926535897931
    #
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-polrot -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # update satomselect to te correct frame
        $__SELECTION_1__ frame $__FRAME__
        $__SELECTION_2__ frame $__FRAME__
        #
        # compute atoms translation to put $__SELECTION_1__ COM to the 0
        [atomselect top "all" frame $__FRAME__] moveby [vecscale -1.0 [measure center "$__SELECTION_1__" weight mass] ]
        #
        # get COM xyz coordinate
        set __COM_xyz__ [measure center "$__SELECTION_2__" weight mass]
        set __COM_x_coor__ [lrange $__COM_xyz__ 0 0]
        set __COM_y_coor__ [lrange $__COM_xyz__ 1 1]
        set __COM_z_coor__ [lrange $__COM_xyz__ 2 2]
        #
        # cartesian to sphere coordiantes
        set __R__ [expr sqrt( ($__COM_x_coor__**2) + ($__COM_y_coor__**2) + ($__COM_z_coor__**2) )]
        set __THETA__ [expr acos($__COM_z_coor__ / $__R__)]
        set __PSI__ [expr atan($__COM_x_coor__ / $__COM_y_coor__)]
        set __THETA_deg__ [expr ($__THETA__ / $PI) * 180]
        set __PSI_deg__ [expr ($__PSI__ / $PI) * 180]
        #
        # write data in the files
        puts $__OUTPUT_FILE__ "$__FRAME__ $__R__ $__THETA__ $__PSI__ $__THETA_deg__ $__PSI_deg__"
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "FC-polrot -- $__PERCENT__ %"
    }  
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-polrot -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
# Unfortunately the variable name  args  cannot be changed
proc FC-sasa {args} {
    #-----------------------------------------------------------------
    set __VAR_ARGUMENTS__ {}
    foreach __ARGUMENT__ $args { lappend __VAR_ARGUMENTS__ $__ARGUMENT__ }
    # set the output file
    set __OUTPUT_FILE_NAME__ [join [lrange $__VAR_ARGUMENTS__ 0 0 ]]
    # Check if Selection_1 is not empty and Selection_2 is empty --------------------- fonctionelle
    if {[expr {[lrange $__VAR_ARGUMENTS__ 1 1 ] ne ""}] && [expr {[lrange $__VAR_ARGUMENTS__ 2 2 ] eq ""}]} {
        set __SELECTION_1__ [join [lrange $__VAR_ARGUMENTS__ 1 1 ]]
        set __SELECTION_2__ [join [lrange $__VAR_ARGUMENTS__ 1 1 ]]
    # Check if Selection_1 and Selection_2 are not empty
    } elseif {[expr {[lrange $__VAR_ARGUMENTS__ 1 1 ] ne ""}] && [expr {[lrange $__VAR_ARGUMENTS__ 2 2 ] ne ""}]} {
        set __SELECTION_1__ [join [lrange $__VAR_ARGUMENTS__ 1 1 ]]
        set __SELECTION_2__ [join [lrange $__VAR_ARGUMENTS__ 2 2 ]]
    }
    #
    #-----------------------------------------------------------------
    #
    # execute the procedure ::sasa with the given arguments
    puts "FC-sasa -- START."
    ::sasa $__OUTPUT_FILE_NAME__ $__SELECTION_1__ $__SELECTION_2__
    puts "FC-sasa -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#===================================================================== 
proc FC-sasa-rl {__OUTPUT_FILE_NAME__ __RECEPTOR__ __LIGAND__} {
    #
    puts "FC-sasa-rl -- START."
    #
    set __COMPLEX__ "$__RECEPTOR__ or $__LIGAND__"
    #----------------------------------------------------------------
    # SASA receptor
    puts "\nFC-sasa-rl -- Step 1/4 : SASA Receptor"
    set __OUTPUT_FILE_1__ "$__OUTPUT_FILE_NAME__-Receptor"
    ::sasa $__OUTPUT_FILE_1__ $__RECEPTOR__ $__RECEPTOR__
    # SASA ligand
    puts "\nFC-sasa-rl -- Step 2/4 : SASA Ligand"
    set __OUTPUT_FILE_2__ "$__OUTPUT_FILE_NAME__-Ligand"
    ::sasa $__OUTPUT_FILE_2__ $__LIGAND__ $__LIGAND__
    # SASA Complex
    puts "\nFC-sasa-rl -- Step 3/4 : SASA Complex"
    set __OUTPUT_FILE_3__ "$__OUTPUT_FILE_NAME__-Complex"
    ::sasa $__OUTPUT_FILE_3__ $__COMPLEX__ $__COMPLEX__
    #
    #----------------------------------------------------------------
    # Temporary SASA file of the receptor with the "receptor/ligand interface hole"
    puts "\nFC-sasa-rl -- Step 4/4 : Interface Receptor/Ligand"
    set __OUTPUT_FILE_5__ [open "$__OUTPUT_FILE_NAME__-RInterface.dat" w]
    puts $__OUTPUT_FILE_5__ "frame SASA(Interface)"
    #
    set __FILE_TEMP_NAME__ "temporary"
    ::sasa $__FILE_TEMP_NAME__ $__COMPLEX__ $__RECEPTOR__
    #
    set __FILE_TEMP_MAX_LINE__ [exec wc -l "$__FILE_TEMP_NAME__-SASA.dat" | cut -d " " -f1 ]
    set __FILE_RECEPTOR_SASA_MAX_LINE__ [exec wc -l "$__OUTPUT_FILE_NAME__-Receptor-SASA.dat" | cut -d " " -f1 ]
    #
    set __COUNTER__ 1
    set p "p"
    #
    if {$__FILE_TEMP_MAX_LINE__ == $__FILE_RECEPTOR_SASA_MAX_LINE__} {
        #
        set __MAX_LINE__ $__FILE_RECEPTOR_SASA_MAX_LINE__
        #
        for {set __LINE_NUMBER__ 2} {$__LINE_NUMBER__ <= $__MAX_LINE__} {incr __LINE_NUMBER__} {
            #
            set __LINE_READ__ "$__LINE_NUMBER__$p"
            #
            set __SASA_FILE_TEMP__ [exec sed -n $__LINE_READ__ "$__FILE_TEMP_NAME__-SASA.dat" | cut -d " " -f2 ]
            set __SASA_RECEPTOR__ [exec sed -n $__LINE_READ__ "$__OUTPUT_FILE_NAME__-Receptor-SASA.dat" | cut -d " " -f2 ]
            set __FRAME__ [exec sed -n $__LINE_READ__ "$__OUTPUT_FILE_NAME__-Receptor-SASA.dat" | cut -d " " -f1 ]
            #
            set __INTERFACE__ [expr {$__SASA_RECEPTOR__ - $__SASA_FILE_TEMP__}]
            # write data in the file
            puts $__OUTPUT_FILE_5__ "$__FRAME__ $__INTERFACE__"
            # print the curent state
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__MAX_LINE__) *100]]
            puts "\t Interface -- $__PERCENT__ %"
        }
        #end for
    } else {
        puts "FC-sasa-rl -- ERROR in the number of lines in the files."
        close $__OUTPUT_FILE_5__
    }
    #end if
    #
    #----------------------------------------------------------------
    #
    exec rm -f temporary-SASA.dat
    close $__OUTPUT_FILE_5__
    #
    puts "FC-sasa-rl -- END."
}
#end proc
#=====================================================================
#=====================================================================
















#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-rmsf {__OUTPUT_FILE_NAME__ __SELECT__} {
    #-----------------------------------------------------------------
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-RMSF.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "AtomNumber RMSF"
    #
    #----------------------------------------------------------------
    #
    ::align "$__SELECT__"
    #
    #----------------------------------------------------------------
    #
    puts "\nFC-rmsf -- START."
    #
    set __ATOM__ 0
    foreach val [measure rmsf [atomselect top "$__SELECT__"]] {
        incr __ATOM__
        puts $__OUTPUT_FILE__ "$__ATOM__ $val"
    }
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-rmsf -- END."
}
#end proc
#=====================================================================
#=====================================================================














#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-rmsd {__OUTPUT_FILE_NAME__ __SELECTION__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-RMSD.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame RMSD"
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __REFERENCE__ [atomselect top "$__SELECTION__" frame $__FRAME_INITIAL__]
    set __COMPARE__ [atomselect top "$__SELECTION__"]
    set __ALL_ATOMS__ [atomselect top all]
    #
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-rmsd -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # upadte atomselects to the the correct frame
        $__COMPARE__ frame $__FRAME__
        $__ALL_ATOMS__ frame $__FRAME__
        # do the alignment
        $__ALL_ATOMS__ move [measure fit $__COMPARE__ $__REFERENCE__]
        # compute the RMSD
        set __RMSD__ [measure rmsd $__COMPARE__ $__REFERENCE__]
        # print the RMSD
        puts $__OUTPUT_FILE__ "$__FRAME__ $__RMSD__"
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\tFC-rmsd -- $__PERCENT__ %"
    } 
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-rmsd -- END."
}
# proc END
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-2drmsd {__OUTPUT_FILE_NAME__ __SELECTION__ __STRIDE__} {
    #-----------------------------------------------------------------
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-2dRMSD.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame_A frame_B RMSD"
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __REFERENCE__ [atomselect top "$__SELECTION__"]
    set __COMPARE__ [atomselect top "$__SELECTION__"]
    set __ALL_ATOMS__ [atomselect top all]
    #
    set __LOOPS_MAX__ [expr (($__FRAME_TOTAL__/ $__STRIDE__)+1)*(($__FRAME_TOTAL__ / $__STRIDE__)+1)]
        # the +1 is mandatory to round up, otherwise the calculations are wrong
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-2drmsd -- START."
    #
    for {set __FRAME_A__ $__FRAME_INITIAL__} {$__FRAME_A__ < $__FRAME_TOTAL__} {incr __FRAME_A__ $__STRIDE__} {
        $__REFERENCE__ frame $__FRAME_A__
        for {set __FRAME_B__ $__FRAME_INITIAL__} {$__FRAME_B__ < $__FRAME_TOTAL__} {incr __FRAME_B__ $__STRIDE__} {
            # upadte atomselects to the the correct frame
            $__COMPARE__ frame $__FRAME_B__
            $__ALL_ATOMS__ frame $__FRAME_B__
            # do the alignment
            $__ALL_ATOMS__ move [measure fit $__COMPARE__ $__REFERENCE__]
            # compute the RMSD
            set __RMSD__ [measure rmsd $__COMPARE__ $__REFERENCE__]
            # print the RMSD
            puts $__OUTPUT_FILE__ "$__FRAME_A__ $__FRAME_B__ $__RMSD__"
            # print the curent state
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__LOOPS_MAX__) *100]]
            puts "\tFC-2drmsd -- $__PERCENT__ %"
        } 
    }
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-2drmsd -- END."
}
# proc END
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-2drmsd-molecules {__OUTPUT_FILE_NAME__ __MOLECULE_1_ID__ __MOLECULE_2_ID__ __SELECTION__ __STRIDE__} {
    #-----------------------------------------------------------------
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-RMSDmolecules.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame_A frame_B RMSD"
    #
    set __REFERENCE__ [atomselect $__MOLECULE_1_ID__ "$__SELECTION__"]
    set __COMPARE__ [atomselect $__MOLECULE_2_ID__ "$__SELECTION__"]
    set __ALL_ATOMS_MOLECULE_2__ [atomselect $__MOLECULE_2_ID__ all]
    #
    # update variables from global
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __FRAME_TOTAL_MOLECULE_1__ [molinfo $__MOLECULE_1_ID__ get numframes]
    set __FRAME_TOTAL_MOLECULE_2__ [molinfo $__MOLECULE_2_ID__ get numframes]
    #
    set __LOOPS_MAX__ [expr (($__FRAME_TOTAL_MOLECULE_1__ / $__STRIDE__)+1) * (($__FRAME_TOTAL_MOLECULE_2__ / $__STRIDE__)+1)]
        # the +1 is mandatory to round up, otherwise the calculations are wrong
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    if {$__FRAME_TOTAL_MOLECULE_1__ == $__FRAME_TOTAL_MOLECULE_2__} {
    puts "FC-rmsd-molecules -- START."
    #
    set __FRAME_TOTAL__ $__FRAME_TOTAL_MOLECULE_1__
    for {set __FRAME_A__ $__FRAME_INITIAL__} {$__FRAME_A__ < $__FRAME_TOTAL__} {incr __FRAME_A__ $__STRIDE__} {
        $__REFERENCE__ frame $__FRAME_A__
        for {set __FRAME_B__ $__FRAME_INITIAL__} {$__FRAME_B__ < $__FRAME_TOTAL__} {incr __FRAME_B__ $__STRIDE__} {
            # upadte atomselects to the the correct frame
            $__COMPARE__ frame $__FRAME_B__
            $__ALL_ATOMS_MOLECULE_2__ frame $__FRAME_B__
            # do the alignment
            $__ALL_ATOMS_MOLECULE_2__ move [measure fit $__COMPARE__ $__REFERENCE__]
            # compute the RMSD
            set __RMSD__ [measure rmsd $__COMPARE__ $__REFERENCE__]
            # print the RMSD
            puts $__OUTPUT_FILE__ "$__FRAME_A__ $__FRAME_B__ $__RMSD__"
            # print the curent state
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__LOOPS_MAX__) *100]]
            puts "\tFC-rmsd-molecules -- $__PERCENT__ %"
        } 
        #end for
    }
    #end for
    } else {
        puts "FC-rmsd-molecules -- ERROR : trajectories for molecules $__MOLECULE_1_ID__ and $__MOLECULE_2_ID__ must have the same frame number !"
    }
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-rmsd-molecules -- END."
}
# proc END
#=====================================================================
#=====================================================================















#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-rgyr {__OUTPUT_FILE_NAME__ __SELECTION__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-RGYR.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame rgyr"
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-rgyr -- START."
    #
    for {set __FRAME__ $__FRAME_INITIAL__} {$__FRAME__ < $__FRAME_TOTAL__} {incr __FRAME__} {
        # upadte atomselects to the the correct frame and calculate rgyr 
        set __RGYR__ [measure rgyr [atomselect top "$__SELECTION__" frame $__FRAME__] weight mass]
        # print the RMSD
        puts $__OUTPUT_FILE__ "$__FRAME__ $__RGYR__"
        # print the curent state
        set __COUNTER__ [expr $__COUNTER__ +1.0]
        set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__FRAME_TOTAL__) *100]]
        puts "\tFC-rgyr -- $__PERCENT__ %"
    } 
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-rgyr -- END."
}
# proc END
#=====================================================================
#=====================================================================















#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-drgyr {__OUTPUT_FILE_NAME__ __SELECTION__ __STRIDE__} {
    #-----------------------------------------------------------------
    # set outfile
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-DRGYR.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame drgyr"
    #
    # update variables from global
    upvar 1 __FRAME_TOTAL__ __FRAME_TOTAL__
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __LOOPS_MAX__ [expr (($__FRAME_TOTAL__/ $__STRIDE__)+1)*(($__FRAME_TOTAL__ / $__STRIDE__)+1)]
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-drgyr -- START."
    #
    for {set __FRAME_A__ $__FRAME_INITIAL__} {$__FRAME_A__ < $__FRAME_TOTAL__} {incr __FRAME_A__ $__STRIDE__} {
        # upadte atomselects to the the correct frame and calculate rgyr
        set __RGYR_FRAME_A__ [measure rgyr [atomselect top "$__SELECTION__" frame $__FRAME_A__] weight mass]
        #
        for {set __FRAME_B__ $__FRAME_INITIAL__} {$__FRAME_B__ < $__FRAME_TOTAL__} {incr __FRAME_B__ $__STRIDE__} {
            # upadte atomselects to the the correct frame and calculate rgyr 
            set __RGYR_FRAME_B__ [measure rgyr [atomselect top "$__SELECTION__" frame $__FRAME_B__] weight mass]
            # calculate adn print the DRGYR
            set __DRGYR__ [expr abs($__RGYR_FRAME_A__ - $__RGYR_FRAME_B__)]
            puts $__OUTPUT_FILE__ "$__FRAME_A__ $__FRAME_B__ $__DRGYR__"
            # print the curent state
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__LOOPS_MAX__) *100]]
            puts "\tFC-drgyr -- $__PERCENT__ %"
        }
        #end for
    }
    #end for
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-drgyr -- END."
}
# proc END
#=====================================================================
#=====================================================================
















#===================================================================== A tester
#=====================================================================
proc FC-drgyr-molecules {__OUTPUT_FILE_NAME__ __MOLECULE_1_ID__ __MOLECULE_2_ID__ __SELECTION__ __STRIDE__} {
    #-----------------------------------------------------------------
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-DRGYRmolecules.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "frame_A frame_B DRGYR"
    #
    # update variables from global
    upvar 1 __FRAME_INITIAL__ __FRAME_INITIAL__
    #
    set __FRAME_TOTAL_MOLECULE_1__ [molinfo $__MOLECULE_1_ID__ get numframes]
    set __FRAME_TOTAL_MOLECULE_2__ [molinfo $__MOLECULE_2_ID__ get numframes]
    #
    set __LOOPS_MAX__ [expr (($__FRAME_TOTAL_MOLECULE_1__ / $__STRIDE__)+1) * (($__FRAME_TOTAL_MOLECULE_2__ / $__STRIDE__)+1)]
        # the +1 is mandatory to round up, otherwise the calculations are wrong
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    if {$__FRAME_TOTAL_MOLECULE_1__ == $__FRAME_TOTAL_MOLECULE_2__} {
    puts "FC-drgyr-molecules -- START."
    #
    set __FRAME_TOTAL__ $__FRAME_TOTAL_MOLECULE_1__
    #
    for {set __FRAME_A__ $__FRAME_INITIAL__} {$__FRAME_A__ < $__FRAME_TOTAL__} {incr __FRAME_A__ $__STRIDE__} {
        set __RGYR_MOLECULE_1__ [measure rgyr [atomselect $__MOLECULE_1_ID__ "$__SELECTION__" frame $__FRAME_A__] weight mass]
        for {set __FRAME_B__ $__FRAME_INITIAL__} {$__FRAME_B__ < $__FRAME_TOTAL__} {incr __FRAME_B__ $__STRIDE__} {
            # upadte atomselects to the the correct frame and calculate rgyr 
            set __RGYR_MOLECULE_2__ [measure rgyr [atomselect top "$__SELECTION__" frame $__FRAME_B__] weight mass]
            # calculate adn print the DRGYR
            set __DRGYR__ [expr abs($__RGYR_MOLECULE_1__ - $__RGYR_MOLECULE_2__)]
            puts $__OUTPUT_FILE__ "$__FRAME_A__ $__FRAME_B__ $__DRGYR__"
            # print the curent state
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__LOOPS_MAX__) *100]]
            puts "\tFC-drgyr-molecules -- $__PERCENT__ %"
        } 
        #end for
    }
    #end for
    } else {
        puts "FC-drgyr-molecules -- ERROR : trajectories for molecules $__MOLECULE_1_ID__ and $__MOLECULE_2_ID__ must have the same frame number !"
    }
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "FC-drgyr-molecules -- END."
}
# proc END
#=====================================================================
#=====================================================================














#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-map-contacts {__OUTPUT_FILE_NAME__ __SELECT__OBJECT_1___ __SELECT__OBJECT_2___} {
    #-----------------------------------------------------------------
    set __OUTPUT_FILE_NAME__ "$__OUTPUT_FILE_NAME__-MapContacts.dat"
    set __OUTPUT_FILE__ [open "$__OUTPUT_FILE_NAME__" w]
    puts $__OUTPUT_FILE__ "selection_1 selection_2 distance_frame_1 distance_frame_2 ..."
    #
    set __OBJECT_1__ [[atomselect top "$__SELECT__OBJECT_1___"] list]
    #
    set __OBJECT_2__ [[atomselect top "$__SELECT__OBJECT_2___"] list]
    #
    #----------------------------------------------------------------
    #
    set __LENGTH__OBJECT_1__ [llength $__OBJECT_1__]
    set __LENGTH__OBJECT_2__ [llength $__OBJECT_2__]
    #
    set __LOOPS_MAX__ [expr $__LENGTH__OBJECT_1__ * $__LENGTH__OBJECT_2__ ]
    set __COUNTER__ 0
    #
    #----------------------------------------------------------------
    #
    puts "FC-map-contacts -- START."
    puts "\n \t Info) This step may cause the software to stop responding. Just wait ! "
    #
    for {set __LIST_NUMBER_OBJECT_1__ 0} {$__LIST_NUMBER_OBJECT_1__ < $__LENGTH__OBJECT_1__} {incr __LIST_NUMBER_OBJECT_1__} {
        #
        set __ATOM__OBJECT_1__ [lindex $__OBJECT_1__ $__LIST_NUMBER_OBJECT_1__]
        #
        for {set __LIST_NUMBER_OBJECT_2__ 0} {$__LIST_NUMBER_OBJECT_2__ < $__LENGTH__OBJECT_2__} {incr __LIST_NUMBER_OBJECT_2__} {
            #
            set __ATOM__OBJECT_2__ [lindex $__OBJECT_2__ $__LIST_NUMBER_OBJECT_2__]
            #
            if {$__ATOM__OBJECT_1__ != $__ATOM__OBJECT_2__} {
                # measure the distance between Atom_object_1 and Atom_object_2
                set __MEASURE_RESULTS__ [measure bond [list $__ATOM__OBJECT_1__ $__ATOM__OBJECT_2__] frame all]
                # write the result in the outfile
                puts $__OUTPUT_FILE__ "[expr $__LIST_NUMBER_OBJECT_1__ +1] [expr $__LIST_NUMBER_OBJECT_2__ +1] $__MEASURE_RESULTS__"
                #    +1 because _list_number_ start at 0
                # print the curent state
            } else {
                puts $__OUTPUT_FILE__ "[expr $__LIST_NUMBER_OBJECT_1__ +1] [expr $__LIST_NUMBER_OBJECT_2__ +1] 0"
                #    +1 because _list_number_ start at 0
                # if __ATOM__OBJECT_1__ = __ATOM__OBJECT_2__ it's tha same atom, so the distance is 0, but VMD cannot do it by itself
            }
            #end if
            set __COUNTER__ [expr $__COUNTER__ +1.0]
            set __PERCENT__ [format "%.2f" [expr ($__COUNTER__ / $__LOOPS_MAX__) *100]]
            puts "FC-map-contacts -- $__PERCENT__ %"
        }
    }
    #
    #----------------------------------------------------------------
    #
    close $__OUTPUT_FILE__
    foreach variables [ info vars ] { if {[string first "atomselect" $variables] != -1} {unset $variables} }
    puts "\n FC-map-contacts -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-clustering {__FOLDER_NAME__ __MASK__ __START__ __STOP__ __STRIDE__ __SEL__} {
    #----------------------------------------------------------------- fonctionnel
    #
    puts "FC-clustering -- START."
    upvar 1 __TERMINAL__ __TERMINAL__
    # create the working folder for CPPTRAJ
    set __HERE__ [exec pwd]
    set __FOLDER_NAME__ "$__FOLDER_NAME__-cluster"
    exec mkdir $__FOLDER_NAME__
    puts "\tFC-clustering -- Folder created : $__FOLDER_NAME__"
    cd $__FOLDER_NAME__
    # create the CPPTRAJ input file
    set __OUTPUT_FILE__ [open "clustering.cpptraj" w]
    #
    #----------------------------------------------------------------- fonctionnel
    #
    # check info on parameter and trajectories files
    set __MOLINFO__ [molinfo top get filename]
    set __PARM_FILE__ [lindex [join $__MOLINFO__] 0]
    #
    # chek if the path is contained on the __PARM_FILE__ variable
    if {[string first "/" $__PARM_FILE__] == -1} {
        # if the __PARM_FILE__ don't contain "/" symbol, it means that the path is not already written 
        puts $__OUTPUT_FILE__ "parm $__HERE__/$__PARM_FILE__\n"
    } elseif {[string first "/" $__PARM_FILE__] != -1} {
        # if the __PARM_FILE__ contain "/" symbol, it means that the path is already written 
        puts $__OUTPUT_FILE__ "parm $__PARM_FILE__\n"
    }
    # end if
    #
    if { [expr {"$__SEL__" eq "first"}] } {
        set __FILE_1__ [lindex [join $__MOLINFO__] 1]
        # chek if the path is contained on the __PARM_FILE__ variable
        if {[string first "/" $__FILE_1__] == -1} {
            # if the __FILE_1__ don't contain "/" symbol, it means that the path is not already written 
            puts $__OUTPUT_FILE__ "trajin $__HERE__/$__FILE_1__ $__START__ $__STOP__ $__STRIDE__"
        } elseif {[string first "/" $__FILE_1__] != -1} {
            # if the __FILE_1__ contain "/" symbol, it means that the path is already written 
            puts $__OUTPUT_FILE__ "trajin $__FILE_1__ $__START__ $__STOP__ $__STRIDE__"
        }
        # end if
        #
        set __ITEM_NB__ 2
        set __READING__ "1 last $__STRIDE__"
    } elseif { [expr {"$__SEL__" eq "all"}] } {
        set __ITEM_NB__ 1
        set __READING__ "$__START__ $__STOP__ $__STRIDE__"
    }
    #end if
    #
    set __TRAJ_FILE__ [lrange [join $__MOLINFO__] $__ITEM_NB__ end]
    #
    foreach __FILE__ $__TRAJ_FILE__ {
        # chek if the path is contained on the __PARM_FILE__ variable
        if {[string first "/" $__FILE__] == -1} {
            # if the __FILE_1__ don't contain "/" symbol, it means that the path is not already written 
            puts $__OUTPUT_FILE__ "trajin $__HERE__/$__FILE__ $__READING__"
        } elseif {[string first "/" $__FILE__] != -1} {
            # if the __FILE_1__ contain "/" symbol, it means that the path is already written 
            puts $__OUTPUT_FILE__ "trajin $__FILE__ $__READING__"
        }
        # end if
    }
    # end foreach
    #
    puts $__OUTPUT_FILE__ "\nautoimage"
    #
    puts $__OUTPUT_FILE__ "\ncluster hieragglo averagelinkage rms $__MASK__ out cnumvtime.dat summaryhalf summary.dat cpopvtime pop.dat normframe repout rep repfmt pdb"
    #
    puts "\tFC-clustering -- CPPTRAJ run file is written."
    close $__OUTPUT_FILE__
    #
    #----------------------------------------------------------------- fonctionnel
    #
    puts "\tFC-clustering -- CPPTRAJ is runing, please wait. ... It can take long time."
    exec cpptraj < clustering.cpptraj > cpptraj.log
    puts "\tFC-clustering -- CPPTRAJ is finished."
    #
    #-----------------------------------------------------------------
    #
    cd $__HERE__
    puts "FC-clustering -- END."
}
#end proc
#=====================================================================
#=====================================================================













#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-reconstruct {__SEL__ __BOXTYPE__} {
     pbc unwrap -sel "$__SEL__" -all
     pbc wrap -$__BOXTYPE__ -centersel "$__SEL__" -center com -compound residue -all
}
#end proc
#=====================================================================
#=====================================================================












#===================================================================== FONCTIONNEL
#=====================================================================
proc FC-help {} {
#
set _output_file_ [open "FC-help.txt" w]
#
puts $_output_file_ "
███████╗██╗███╗   ██╗██████╗        ██████╗ ██████╗ ███╗   ██╗████████╗ █████╗  ██████╗████████╗███████╗
██╔════╝██║████╗  ██║██╔══██╗      ██╔════╝██╔═══██╗████╗  ██║╚══██╔══╝██╔══██╗██╔════╝╚══██╔══╝██╔════╝
█████╗  ██║██╔██╗ ██║██║  ██║█████╗██║     ██║   ██║██╔██╗ ██║   ██║   ███████║██║        ██║   ███████╗
██╔══╝  ██║██║╚██╗██║██║  ██║╚════╝██║     ██║   ██║██║╚██╗██║   ██║   ██╔══██║██║        ██║   ╚════██║
██║     ██║██║ ╚████║██████╔╝      ╚██████╗╚██████╔╝██║ ╚████║   ██║   ██║  ██║╚██████╗   ██║   ███████║
╚═╝     ╚═╝╚═╝  ╚═══╝╚═════╝        ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═╝   ╚══════╝

v2.3 \t by Tom MICLOT \t :: HELP FILE ::


WHAT CAN YOU DO 
\t * Make networks and statistic tables of (Hydrogen Bonds) and (Salt Bridges) and (Hydrophobic Interactions) for protein and nucleic.
\t   But only HBonds interactions is able for other molecules.
\t * Compare interactions between multiple simulations.
\t * Draw 3D interaction network in VMD.
\t   Enhance interacting residues as 3D spheres in VMD.
\t * Calculate a contact map betwen tow selections (based on distance).
\t * Calculate the radius of gyration and the diffrerence in RGYR.
\t * Calculate the RMSF.
\t * Calculate the RMSD-2D map.
\t * Calculate the RMSD between tow molecules.
\t * Calculate SASA for a selection.
\t * Calculate SASA for (Receptor), (Ligand) and (Interface Receptor-Ligand).
\t * Make a clustering with CPPTRAJ.

INFORMATIONS
\t * Always use \" \" to write each argument.
\t * All output files are in .dat , it correspond to text file. So they can be read, plot, etc ... very easily.
\t * The script always use the TOP molecule, else for the FC-rmsd-molecules command.
\t * Sometimes the script causes the VMD window to be unresponsive. Just wait for the script to finish working.

DEPENDENCY 
\t * CPPTRAJ    -- required only for clustering.

USAGE
\t * FC-help                     -- display this help.
\t * FC-receptor-ligand          -- interactions networks and tables of receptor/ligand, along the simulation.
\t * FC-protein-changes          -- interactions networks and tables of protein changes, along the simulation.
\t * FC-network-trace-bonds      -- draw the interaction network bonds.
\t * FC-network-trace-resids     -- darw resid in interaction network as spheres.
\t * FC-interactions-similitudes -- find the conserved interactions between multiples simulations, from network data.
\t * FC-neighbors                -- give the frequence of residues in selection_1 within CutOff of selection_2.
\t * FC-comcom                   -- measure the COM-COM distance of two selection, along the simulation.
\t * FC-polrot                   -- gives the polar coordinates of selection 2 with respect to selection 1
\t * FC-map-contacts             -- calculate the contact map between tow selections, along the simulation.
\t * FC-sasa                     -- calculate SASA of a selection, along the simulation.
\t * FC-sasa-rl                  -- calculate SASA and Interface for receptor/ligand, along the simulation.
\t * FC-rmsf                     -- calculate the RMSF of a selection, along the simulation.
\t * FC-rmsd                     -- calculate the RMSD of a selection, along the simulation.
\t * FC-2drmsd                   -- calculate the 2D-RMSD of a selection.
\t * FC-2drmsd-molecules         -- calculate the RMSD between tow molecules, along their simulations.
\t * FC-rgyr                     -- calculate the radius of gyration of a selection, along the simulation.
\t * FC-drgyr                    -- calculate the difference in radius of gyration of a selection.
\t * FC-drgyr-molecules          -- calculate the DRGYR between tow molecules, along their simulations.
\t * FC-clustering               -- make a clustering of the simulation. (CPPTRAJ)
\t * FC-renumber                 -- renumer resid of a selection from a given number.
\t * FC-reconstruct              -- center the selection and remove atom coordinate jumps.

VARIABLES (can be used by user)
\t Info) Find-Contacts use equivalent selection as GetContacts : getcontacts.github.io
\t * protein_SB_Pos -- ((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) or ((resname LYS) and (name NZ)) or ((resname ARG) and (name NH1 NH2))
\t * protein_SB_Neg -- ((resname ASP) and (name OD1 OD2)) or ((resname GLU) and (name OE1 OE2))
\t * protein_SB_All -- ((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) or ((resname LYS) and (name NZ)) or ((resname ARG) and (name NH1 NH2)) or ((resname ASP) and (name OD1 OD2)) or ((resname GLU) and (name OE1 OE2))
\t * protein_Hyd    -- hydrophobic and not backbone and (type C C1 C2 CA CB CC CE CI CK CQ CR CT CW)
\t * nucleic_SB     -- (name OP1 OP2)
\t For other selection type, please visit     https://www.ks.uiuc.edu/Research/vmd/current/ug/node90.html



======================================================================
::: FC-receptor-ligand :::

This command compute interaction search and analysis of H-Bonds, Salt bridges and Hydrophobic interactions between Receptor and Ligand, along the simulation. At the end you will get files to draw the interaction-type network, and others corresponding to interaction-type statistics for the min_frequence. 
The min_frequence parameter select residue pairs with a frequence >= min_frequence. The value must be between 0 and 1, and can take three decimal places.

USAGE
FC-receptor-ligand \"OutputFileName\" \"receptor\" \"ligand\" \"min_frequence\"

OUTFILE
\t * The outfiles will be OutputFileName-HBonds-Network-Frequences.dat , OutputFileName-HBonds-TableStats.dat , OutputFileName-Hydrophobic-Network-Frequences.dat , OutputFileName-Hydrophobic-TableStats.dat , OutputFileName-SaltBridges-Network-Frequences.dat , OutputFileName-SaltBridges-TableStats.dat .
\t * Type Network is use to draw network bonds and residues, but also to get interactions similitudes between multiple simulations.
\t   Type TableStats is a table with statistical analysis of interactions.
\t * All outfiles are column based files.
\t * In the TableStats files, interaction are selected by the frequence (>= minimum frequency) and by the mean of the atomic distance (<= 4 angstroms).


WARNINGS
\t * It may happen that the statistical table files are empty.
\t   This means that the frequency of interactions are all lower than the minimum frequency chosen, or the mean of distance is > 4 angstroms.
\t   In this case you should reduce the value of the minimum frequency.
\t * For protein-X interaction search, please always use \"receptor\" for your (protein/peptide/amono acids) selection.
\t * protein-protein : H-Bonds, Salt bridges and Hydrophobic available.
\t   protein-nucleic : H-Bonds and Salt bridges available.
\t   nucleic-nucleic : only H-Bonds available.
\t   for others      : only H-Bonds available.



======================================================================
::: FC-protein-changes :::

This command compute interaction search and analysis of H-Bonds, Salt bridges and Hydrophobic interactions of protein changes along the simulation. At the end you will get files to draw the interaction-type network, and others corresponding to interaction-type statistics for the min_frequence. 
The min_frequence parameter select residue pairs with a frequence >= min_frequence. The value must be between 0 and 1, and can take three decimal places.

USAGE
FC-protein-changes \"OutputFileName\" \"min_frequence\"

OUTFILE
\t * The outfiles will be OutputFileName-HBonds-Network-Frequences.dat , OutputFileName-HBonds-TableStats.dat , OutputFileName-Hydrophobic-Network-Frequences.dat , OutputFileName-Hydrophobic-TableStats.dat , OutputFileName-SaltBridges-Network-Frequences.dat , OutputFileName-SaltBridges-TableStats.dat .
\t * Type Network is use to draw network bonds and residues, but also to get interactions similitudes between multiple simulations.
\t   Type TableStats is a table with statistical analysis of interactions.
\t * All outfiles are column based files.
\t * In the TableStats files, interaction are selected by the frequence (>= minimum frequency) and by the mean of the atomic distance (<= 4 angstroms).

WARNING
\t * It may happen that the statistical table files are empty.
\t   This means that the frequency of interactions are all lower than the minimum frequency chosen, or the mean of distance is > 4 angstroms.
\t   In this case you should reduce the value of the minimum frequency.



======================================================================
::: FC-network-trace-bonds :::

This command draw the network from the (Input File), at the curent frame.
It take 4 parameters :
\t * InputFileName -- name of the file to be read. (Ex: xxx-HBonds-TableStat-Strong.dat)
\t * expansion     -- used to grow, or reduce, radius of network cylinder.
\t * MinRadius     -- select minimum radius you want to display
\t * color         -- color of the networ cylinder. All VMD color name are are available. (red , green , blue ,...).

USAGE
FC-network-trace-bonds \"InputFileName\" \"MinRadius\" \"expansion\" \"color\"

You can clear all the drawed object with the command :
     draw delete all



======================================================================
::: FC-network-trace-resids :::

This command draw residues involved in the network from the (Input File) as spheres, at the curent frame.
It take 5 parameters :
\t * InputFileName -- name of the file to be read. (Ex: xxx-HBonds-TableStat-Strong.dat)
\t * expansion     -- used to grow, or reduce, radius of network cylinder.
\t * MinRadius     -- select minimum radius you want to display
\t * color         -- color of the networ cylinder. All VMD color name are are available. (red , green , blue ,...).
\t * 1/2           -- chose 1 or 2. In the InputFile residues are written in (resid in column 1) and (resid in column 2) format.

USAGE
FC-network-trace-resids \"InputFileName\" \"MinRadius\" \"expansion\" \"color\" \"1/2\"

You can clear all the drawed object with the command :
     draw delete all




======================================================================
::: FC-comcom :::

This command measure the COM-COM distance evolution along the simulation.

USAGE
FC-rmsd \"OutFileName\" \"selection_1\" \"selection_2\"

OUTFILE
\t * The outfile will be OutFileName-COMCOM.dat
\t   So if you write MyFile as file name, the file will be MyFile-COMCOM.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame and COM-COM-distance




======================================================================
::: FC-polrot :::

This command measure the polar coordinates of selection 2 with respect to selection 1. This the polar

USAGE
FC-polrot \"OutFileName\" \"selection_1\" \"selection_2\"

OUTFILE
\t * The outfile will be OutFileName-POLROT.dat
\t   So if you write MyFile as file name, the file will be MyFile-POLROT.dat
\t * The outfile is column based file.
\t   Columns are respectively Frame, R, Theta angle (polar), Psi angle (polar), Theta angle (degree), Psi angle (degree)



======================================================================
::: FC-sasa :::

This command calculate the SASA of the selection. You can aslo specify a restriction parameter, but it's not necessary. It's similar to the sasa.tcl script available on the VMD website.
<< The restrict option can be used to prevent internal protein voids or pockets from affecting the surface area results. >>

USAGE
FC-sasa \"OutFileName\" \"selection\" \"restriction\"

OUTFILE
\t * The outfile will be OutFileName-SASA.dat
\t   So if you write MyFile as file name, the file will be MyFile-SASA.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame and SASA



======================================================================
::: FC-sasa-rl :::

This command calculate the SASA for receptor and ligand interaction.
It compute (SASA Receptor) , (SASA Ligand) , (SASA Receptor/Ligand complex) and (Interface of Receptor/Ligand)

USAGE
FC-sasa-rl \"OutFileName\" \"receptor\" \"ligand\"

OUTFILE
\t * The outfiles will be OutFileName-Receptor-SASA.dat, OutFileName-Ligand-SASA.dat, OutFileName-Complex-SASA.dat, OutFileName-RLInterface.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame and SASA



======================================================================
::: FC-rmsf :::

This command calculate the RMSF for each atom in the selection. The alignment is automaticaly done before RMSF calculations.

SELECTION EXEMPLES
\t * Protein selection - protein and name CA
\t * Nucleic selection - nucleic and name N1

USAGE
FC-rmsf \"OutFileName\" \"selection\"

OUTFILE
\t * The outfile will be OutFileName-RMSF.dat
\t   So if you write MyFile as file name, the file will be MyFile-RMSF.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively AtomNumber and RMSF




======================================================================
::: FC-rmsd :::

This command calculate the RMSD of the selection. The alignment is automaticaly done during RMSD calculations.

USAGE
FC-rmsd \"OutFileName\" \"selection\"

OUTFILE
\t * The outfile will be OutFileName-RMSD.dat
\t   So if you write MyFile as file name, the file will be MyFile-RMSD.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame and RMSD



======================================================================
::: FC-2drmsd :::

This command calculate the 2D-RMSD of the selection. The alignment is automaticaly done during 2D-RMSD calculations.
You can use a graphing utility like Gnuplot to display an heatmap.
For exemple, in Gnuplot you can wrote:     plot 'file' using 1:2:3 with image

USAGE
FC-2drmsd \"OutFileName\" \"selection\" \"stride\"

OUTFILE
\t * The outfile will be OutFileName-2dRMSD.dat
\t   So if you write MyFile as file name, the file will be MyFile-2dRMSD.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame_A and Frame_B.
\t   Column 3 is the RMSD value.



======================================================================
::: FC-rmsd-molecules :::

This command calculate the RMSD between the selection into tow trajectories. It's an 2D-RMSD like process.
The alignment of the tow molecule is automaticaly done during the calculations
For exemple,y ou can plot an heatmap or make a histogram distribution.

USAGE
FC-rmsd-molecules \"OutFileName\" \"Molecule_1_ID\" \"Molecule_2_ID\" \"selection\" \"stride\"

OUTFILE
\t * The outfile will be OutFileName-RMSDmolecules.dat
\t   So if you write MyFile as file name, the file will be MyFile-RMSDmolecules.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame_A and Frame_B.
\t   Column 3 is the RMSD value.




======================================================================
::: FC-rgyr :::

This command calculate the radius of gyration of the selection.

USAGE
FC-rgyr \"OutFileName\" \"selection\"

OUTFILE
\t * The outfile will be OutFileName-RMSD.dat
\t   So if you write MyFile as file name, the file will be MyFile-RMSD.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame and RGYR



======================================================================
::: FC-dgyr :::

This command calculate the difference in radius of gyration of the selection.
You can use a graphing utility like Gnuplot to display an heatmap.
For exemple, in Gnuplot you can wrote:     plot 'file' using 1:2:3 with image

USAGE
FC-drgyr \"OutFileName\" \"selection\" \"stride\"

OUTFILE
\t * The outfile will be OutFileName-2dRMSD.dat
\t   So if you write MyFile as file name, the file will be MyFile-DRGYR.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame_A and Frame_B.
\t   Column 3 is the DRGYR value.



======================================================================
::: FC-drgyr-molecules :::

This command calculate the DRGYR between the selection into tow trajectories. It's an 2D-RGYR like process.
For exemple,y ou can plot an heatmap or make a histogram distribution.

USAGE
FC-drgyr-molecules \"OutFileName\" \"Molecule_1_ID\" \"Molecule_2_ID\" \"selection\" \"stride\"

OUTFILE
\t * The outfile will be OutFileName-DRGYRmolecules.dat
\t   So if you write MyFile as file name, the file will be MyFile-DRGYRmolecules.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Frame_A and Frame_B.
\t   Column 3 is the DRGYR value.



======================================================================
::: FC-map-contacts :::

This command measure the distance between atoms in (selection 1) and (selection2) for all frames in the MD.
You can use a graphing utility like Gnuplot to display an heatmap, or make an animated of heatmap.
For exemple, in Gnuplot you can wrote:     plot 'file' using 1:2:345 with image

USAGE
FC-map-contacts \"OutFileName\" \"selection 1\" \"selection 2\"

SELECTION EXEMPLES
\t * Protein selection - protein and name CA
\t                       resid 26 to 251 and name CA
\t * Nucleic selection - nucleic and name N1
\t                       resid 1 to 24 and name N1
\t * Ligand  selection - resname LIG
\t                       resid 124

OUTFILE
\t * The outfile will be OutFileName-MapContacts.dat
\t   So if you write MyFile as file name, the file will be MyFile-MapContacts.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are Resid (for protein or nucleic) or Atom Number (for ligand).
\t   Column >= 3 are distance value for one frame. Ex: column 3 is all distance for frame 0, column 4 for frame 1, etc.

WARNINGS
\t * It always use all defined atoms in the selection.
\t   So, for ligand molecule you will get all atoms. In that case, atoms are numbered according to their order in the selection, but not according to their index.
\t * For protein selection, always use : CA .
\t * For nucleic selection, always use : N1 or C1' .



======================================================================
::: FC-clustering :::

This command use CPPTRAJ to calculate a cluster of structures in the simulation. So, it can take a long time.
The command automatically recognizes the parameter file and the trajectory files loaded in VMD. 
Ensure you have load Amber, or AmberTools, in the terminal were you run VMD.
It take 6 drawing parameters :
\t * FolderName -- name of the CPTRAJ work dir. The final name will be FolderName-cluster .
\t * mask       -- mask selection, in CPPTRAJ format. See  https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/
\t * start      -- frame to begin reading.
\t * stop       -- frame to stop reading at.
\t * stride     -- stride for reading in trajectory frames. Stride is applied to all trajectory files loaded (the first and the others).
\t * first/all  -- \"first\" applies (start_value stop_value stride_value) to the first trajectory file loaded, and (1 last stride_value) to the others.
\t              -- \"all\" applies (start_value stop_value stride_value) to the all trajectory files loaded.

USAGE
FC-clustering \"FolderName\" \"mask\" \"start\" \"stop\" \"stride\" \"first/all\"

OUTFILE
\t * The outfiles are the file to run CPPTRAJ, and all output file of this software.

WARNING
\t * Please take into account that VMD has a command (measure cluster), but which is less efficient than the analysis by CPPTRAJ.



======================================================================
::: FC-neighbors :::

This command give the frequence of residues in (selection_1) within (CutOff) of (selection_2).

USAGE
FC-neighbors \"OutFileName\" \"selection_1\" \"selection_2\" \"CutOff\"

OUTFILE
\t * The outfile will be OutFileName-NEIGHBORS.dat
\t   So if you write MyFile as file name, the file will be MyFile-NEIGHBORS.dat
\t * The outfile is column based file.
\t   Column 1 and 2 are respectively Residue , Count in the MD, Frequence



======================================================================
::: FC-renumber :::

This command renumber the IDs of the residues of the (selection), starting from the (startIDnumber).

USAGE
FC-renumber \"selection\" \"startIDnumber\"

WARNING
\t * It is possible that the Sequence Viewer in VMD does not show the correct numbers after using this command.
\t   So please use :





======================================================================
::: FC-reconstruct :::

This command uses a wrap/unwrap protocol to cleanly remove large jumps in atom coordinates of the selection (and the errors this causes in the bonds), and put the selection in the center of the cell.
You can chose the parallelepiped or orthorhombic option to << wrap the atoms into the unitcell parallelepiped or the corresponding orthorhombic box with the same volume and center as the (non-orthrhombic) unitcell. >>

USAGE
FC-reconstruct \"selection\" \"parallelepiped | orthorhombic\"

WARNING
\t * If you don't know which option to choose, you can use parallelepiped as defult.
\t * Using FC-reconstruct command may cause deformation of the water box.
"
#
#
close $_output_file_
upvar 1 __TERMINAL__ __TERMINAL__
exec $__TERMINAL__ -e "vi FC-help.txt ; bash"
after 1000
exec rm -f FC-help.txt
}
#end proc 
#=====================================================================
#=====================================================================












#=====================================================================
#=====================================================================
puts "\t Info) Find-Contacts is ready !"
puts "\n \t Info) Please don't forget to unwrap your system (or use FC-reconstruct) before using Find-Contacts !"
puts "\n\t Info) Write  FC-help  to see the help file. \n"
#=====================================================================
#=====================================================================
