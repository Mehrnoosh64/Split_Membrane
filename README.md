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
5. Split_Splipids_2020.ff folder
                                                                      
The current version supports systems containing PSM, POPC, POPS, POPE, POPA, and Cholesterol.
The model is easily extendable to other lipid types.

***
If you are using the script, please cite us:                                            
                                                                                                       
Written by Mehrnoosh Kh. Hazrati, Ph.D.                                                                
Robert Vácha Lab, CEITEC – Central European Institute of Technology, Masaryk University     
