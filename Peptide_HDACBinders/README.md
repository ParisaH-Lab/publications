# General requirements
## Rosetta Software
Rosetta software is free for academic use. You can get the license from the following link:     
https://www.rosettacommons.org/software/license-and-download
The code provided in this paper have been updated to run on Rosetta version before February 2019. If you want to run on versions after that, you need to update resfile usage to use `PackerPalette`.

To install and build, follow the protocols in the link below: https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build
Please note that with 1-4 cores, it will take *hours* to fully build Rosetta on a desktop.

## Pyrosetta
You may download and build pyrosetta after obtaining academic license (see above) by following protocols in this page:
http://www.pyrosetta.org/dow
The code provided in this script has been tested with python 3.6 and pyrosetta 2018.42. 

Installation of pyrosetta is rather fast.

## Other packages and softwares
1. Python (version 3.6 works for most of the code. Code in method 1 is older and was tested on python 2.7)
2. Numpy and Pandas
3. Jupyter Notebook (only required to run the demo for Docking analysis)

# Protocols
This folder contains subdirectories with code and detailed isntruction on how to run each of four design methods. Below is a brief description of what each folder contains

## general_files
Includes all the files that are used in more than one folder:
1. resfiles (with and without ncAA)
2. params file for SHA
3. composition files

## Method1
Includes:
1. `README.md` with details of how to run analysis
2. **test** folder with sample input files
3. Codes necessary to run the docking part (`main.py`, `pdbfile.py`, `translate_rotate.py`, `merge_8.xml`) as well as samples on how to run this in large scale (`merge_all.sh`). This step often takes *seconds to minutes* for one scaffold and one target. However, it needs to be performed on many different orientations of the anchor and different scaffolds. Depending on your scaffold and targets, it can take hours to go through everything.
4. The design script, `Method1_Design_7mers.xml`, and a script to batch submit multiple design jobs, `Method1_Design_Loop_7mers.sh`. The design steps takes *few minutes*

## Method2
Includes:
1. `README.md` with details of how to run analysis
2. **test** folder with sample input files
3. Design script `Method2_inSitu_KIC.xml`. This script both extends the peptide chain and designs. If the extension doesn't result in a closed solution, it will fail within *seconds*. A successful extension and design round often takes *several minutes*

## Method3_4
Includes:
1. `README.md` with details of how to run analysis
2. **test** folder with sample input files
3. The code required to exhaustively sample the torsions and residues in  residues immediately before and after anchor, `Method3-4_Exhaust.py`. Because the code goes through a large number of potential phi and psi angles, we usually run this *overnight*. The output of this stage is a csv file with the score of the pose, shape complementarity, computed ∆∆G (with or without repacking), number of unsatisfied hydrogen bonds at the interface, and peptide energy.
4. The code required for exhaustive sampling (step 3) as well as design in method 4 (`Method4_Design.py`, `Method4_*.xml`). This design method goes through three different design runs, thus it takes *hours* for it to finish 

# Output
Because the processes of design and packing are stochastic, the output files can vary from run to run. Thus, we did not provide expected outputs. However, for each method, we clearly mentioned what the output should be if it is not a PDB or a silent file.
