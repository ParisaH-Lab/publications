## Requirements:
1. Rosetta software, any versions from 01/2018-02/2019. In order for it to work for versions after that, you will need to check `PackerPalette` in Rosetta.
2. Some of the files required for running this script are provided in **general_files** folder.

## Running
Note: This code is designed for the specific case of this paper. The numberings in the code need to be modified according to your input.

The code for expanding the SHA-Trp and design is `Method2_inSitu_KIC.xml`. The xml document is fully annotated.

To run xml:

`<path_to_Rosetta>/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-parser:protocol Method2_inSitu_KIC.xml \
-in:file:s <input_file> \
-score:weights <desired_score> \
-nstruct <number_of_runs> \
-ignore_waters false \
-auto_setup_metals true \
-out:file:silent <output_name> `

## Example
The start input file can be found in **test** folder. You need to prepare the example first by running this command:

`<path-to-Rosetta>/main/bin/relax.default.linuxgccrelease -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -score:weights ref2015.wts -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -auto_setup_metals true -s <input_pdb> -ignore_waters false -auto_setup_metals true -extra_res_fa SHA.params`

Sometimes, the addition of Zn ion and waters may break the code. In that case, you can simply remove the zinc ion and waters from the input file and the code. The xml is written such that it maintain the orientation of SHA and Trp even without waters and zinc being in the input.
