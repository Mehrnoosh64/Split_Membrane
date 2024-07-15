# Split Membrane

***
How to make a split bilayer:
You need to run "Splitter.sh".
To run the script, you need to have the following files in the working folder: 
1. input.gro 
2. input.tpr
    (.gro and .tpr files of your desired membrane made by charmm36 force field. 
    For example, you can use step7_10.gro and step7_10.tpr obtained from charmm-gui website.)                                      
3. initial_topol.top                                                                                      
4. water_deletor.pl (Written by Justin Lemkul, jalemkul@vt.edu)                                           
5. Split_Slipids_2020.ff folder
                                                                      
The current version supports systems containing PSM, POPC, POPS, POPE, POPA, and Cholesterol.
The model is easily extendable to other lipid types.


# Reuniting the split membrane

***
If you want to reunite the split membrane after your simulation, you will need to use the merger.py code. If you did not use any restrictions on the head and tail groups to prevent flip-flopping during the simulation, it is likely that some head/tail groups moved to another leaflet. I have provided two different Python codes, one of which may work for your system. 
After reuniting the head and tails, I suggest a minimization with loose restraints on the head group P atoms. For example, for POPC, you can add the following lines to the .itp file

#ifdef POSRES_POPC
[ position_restraints ]
20 1 10 10 10
#endif

Then you have to equilibrate the system. I turned off the restraints for the z direction while keeping them for the xy direction during the equilibration. This way you don't let the head groups move much from their positions.

Note: If you have proteins in your systems, you will need to add them to the merged.gro manually.


***
If you are using the script, please cite us:                                            
                                                                                                       
Written by Mehrnoosh Kh. Hazrati, Ph.D.                                                                
Robert Vácha Lab, CEITEC – Central European Institute of Technology, Masaryk University     
