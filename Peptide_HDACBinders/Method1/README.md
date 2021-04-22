## Requirments
1. The code in this step is updated to run with python2.7 and Rosetta versions  before 02/2019. If you want to run them on a more recent version, you may need to make some modifications.
2. Some of the necessary files are provided in **general_file** folder.

## Step1. Docking:
Docking script is `main.py`. It aligns a residue based on C, CA, N, O position and currently does not work for Glycine. It goes through all residues in the cyclic peptide and places each residue such that the backbone atoms are aligned with SHA backbones. It then outputs the aligned pose.

You will need `translate_rotate.py` and `pdbfile.py` codes for main to function.

You run the code as:
`main.py -f <target-protein> --targetresidues <ALA #1>" -q <query peptide> --queryresidue "<ALA #2>"`

This rotates the peptide such that #2 residue in the peptide is aligned on #1 residue in the target. Please note that all peptide residues and the anchor should be changed to Ala for the script to work for simplicity.

After this, one needs to merge the transposed peptide with the target empty pose. The provided `merge_8.xml` script represents an example on how to do it for 8mer peptides.

`merge_all.sh` shows how to perform a batch alignment from a directory of potential scaffolds onto protein templates.

## Step2. Design:
Given a docked scaffold, you can run `Method1_Design_7mers` file to design. Please note that this design is for 7mers. In order to design for other sizes, the residue numbers where the bond is declared needs to get changed. In newer versions of Rosetta, you can use the mover `PeptideCyclizeMover` without requiring to mention residue names. Examples of this mover are presented in other design methods.

The `Method1_Design_Loop_7mers` file describes the details of how you can run this on a large dataset of docked conformations. 
It will generate the **%%num%%** and **%%k_n%%** variables in the xml.
Please note that this file is hard-coded for the way the docking step generates outputs.

Both scripts are annotated.

## Example
The folder **test** includes a test scaffold, and a test target protein. 
