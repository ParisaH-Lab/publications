s=$1

for b in $(cat binder2.list)
do
    for ((i=1; i<=12; i++))
    do
        echo $b $i
        cp scaffolds/$b .
        cp template/$s .;
        python main.py -f $s --targetresidues "ALA 370" -q $b --queryresidue "ALA $i"
        rm -rf $b
        rm -rf $s
        mv realign${b} ${s}_realign${b}.res${i}.pdb
        cp template/${s}_delete.pdb ${s}_realign${b}.res${i}_input.pdb; Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s ${s}_realign${b}.res${i}_input.pdb -in:file:fullatom -parser:view -parser:protocol merge.xml -parser:script_vars alnpdb=${s}_realign${b}.res${i}.pdb -nstruct 1 -mute all -overwrite; rm ${s}_realign${b}.res${i}_input.pdb -rf; rm ${s}_realign${b}.res${i}.pdb -rf;
    done
done
