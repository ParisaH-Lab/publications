#!/bin/bash

# Defining global variables
## simple counter
COUNTER2=0
## Number of SHA residue in the pocket
COUNTER1=371
## User-input variable of a given comformation
list_num=$1

# Generating clean directories to output results in
mkdir des_${list_num}
cd des_${list_num}/


# Reading through the scaffolds that are docked
while read b
do

    # Gets the size of peptide
    ## Defines this way because of the way the names are originally set
    size=$(echo $b | awk 'BEGIN {FS= "_"} {print $1}')

    # Gets the name of pdb
    pdb_name=$(echo $b | cut -f 1 -d '_' --complement)

    # Loops through different docks of a given scaffold
    for t in dock/out_7/round1/${b}/*pdb; do

        # Laborious way to generate resfiles for each PDB
        while read i
        do

            resn=$(echo $i | awk '{print $4}') # residue name column
            resi=$(echo $i | awk '{print $6}') # residue ID column
            B=$(echo $resi $COUNTER2 | awk '{print $1 - $2}')

            if [ $B != '0' ] #to make sure we do this once for each residue position
            then
                if [ $COUNTER2 == 0 ] # initial parts of resfile
                then
                    echo "NATRO"
                    echo
                    echo "start"
                    echo
                fi

                let COUNTER2=COUNTER2+1

                new=$(echo $resi $COUNTER1 | awk '{print $1 + $2}')

                # Ensures that the positions that are D stay D and positions that are L stay L.
                ## Important for keeping the structure.

                if [[ $resn == D* ]] # No canonical AA starts with D!
                then
                    if [ $resn == "DPR" ]
                    then
                        echo $new" C EMPTY NC DPR" #keeps DPRO the way it is
                    else
                        echo $new" C EMPTY NC DAL NC DAS NC DCS NC DPH NC DGU NC DHI NC DIL NC DLY NC DLE NC DME NC DAN NC DGN NC DAR NC DSE NC DTH NC DVA NC DTR NC DTY"
                    fi

                elif [ $resn == "PRO" ]
                then
                    echo $new" C PIKAA P" # keeps PRO the way it is

                elif (($resi == $size)) #last resdiude is SHA based on dock
                then
                    echo $new" C EMPTY NC SHA"
                else
                    echo $new" C NOTAA CGP" # only restriction for L-AA
                fi

            fi

        done < 7mer/${pdb_name}.pdb > des_${list_num}/resfile_7.txt

        let COUNTER2=0 #prepare for next round

        <Path-To-Rosetta-bin>/rosetta_scripts.default.linuxgccrelease -in:file:s $t  -in:file:fullatom -parser:view -parser:protocol design_7aa.xml -parser:script_vars num=$((${size}+371)) -nstruct 1 -extra_res_fa SHA.params

        rm -rf resfile_7.txt

    done

done < list_${list_num}.txt
