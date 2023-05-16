# Split_Membrane

***
How to make a split bilayer:
You need to run "Splitter.sh".
To run the script you need to have following files in the working folder: 
1. step7_10.gro 
2. step7_10.tpr
    (.gro and .tpr files of your desired membrane made by charmm36 force field. 
    Default files are step7_10.gro and step7_10.tpr obtained from charmm-gui website. 
    If you use other .gro and .tpr files, edit the name in the script)                                      
3. initial_topol.top                                                                                      
4. water_deletor.pl (Written by Justin Lemkul, jalemkul@vt.edu)                                           
5. Split_Splipids_2020.ff folder
                                                                      
Current version supports systems containing PSM, POPC, POPS, POPE, POPA, and Cholesterol.
The model is easily extendable to other lipid types.

***
If you are using the script, please cite this paper:      

???                                             
                                                                                                       
Written by Mehrnoosh Kh. Hazrati, Ph.D.                                                                
Robert Vácha Lab, CEITEC – Central European Institute of Technology, Masaryk University     
