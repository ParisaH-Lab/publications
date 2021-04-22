## Requirements
1. Python 3 or above (tested up tp 3.6)
2. pyrosetta 2018.4
3. Rosetta (any version from 01/2018-02/2019)
4. Certain files provided in **general_files** folder as well as privded residue_4.csv file

## Method 3 protocol
1. Start by exhaustive sampling of residues around native using `Method3-4_Exhaust.py`. You need to run this script once for the residue after the anchor and once for the residue before.
2. The outputs from step 1 will then be used in this step. After runnung `Method3-4_Exhaust.py`,  a csv file will be returned with scores, interface metrics, and torsion value information. The torsion values for the best poses will be used as variables in the design step. You can design by running `Method3_all_in_one.xml`. The script provided is currently hardcoded for generating peptides with 9 residues. It reads in `Method3_GeneralDesign.xml`, so you need that file as well. You can run the script using: 

`<path_to_Rosetta>/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-parser:protocol Method3_all_in_one.xml \
-in:file:s <input_file> \
-score:weights <desired_score> \
-nstruct <number_of_runs> \
-out:file:silent <output_name>`

Both scripts are fully annotated.

## Method 4 protocol
1. Start by exhaustive sampling using `Method3-4_Exhaust.py` and then running multiple rounds of minimization on peptide chain. The best scoring outcomes where selected later for design.
2. The best outputs from last step will be used as inputs for this step. Design using `Method4_Design.py`. In order for the deisgn to work you will need all the files labeled `Method4_*.xml`. The files describe different modes of design for residues at different types of interfaces.

Both scripts are fully annotated. The arguments are detailed in python scripts.

## Example
The example input file with SHA is provided in **test** folder.

Note: these scripts are hard-coded for specific SHA and HDAC2 case.
Note: the script currently only works if the peptide chain is chain 1 and the pose has two chains.
