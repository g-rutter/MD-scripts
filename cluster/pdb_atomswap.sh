#!/bin/sh

PDB=$1

lines_per_entry=`awk '/ENDMDL/{print NR; exit}' $PDB`
atoms=$((lines_per_entry-2))
molecules=$((`awk '/ATOM *'$atoms' /{print $5; exit}' $PDB`+1))
residues=$((`awk '/ATOM *'$atoms' /{print $6; exit}' $PDB`))

atoms_per_molecule=$((atoms/molecules))
residues_per_molecule=$((residues/molecules))

n_lines=`wc -l $PDB | awk '{print $1}'`
n_entries=$((n_lines/lines_per_entry))

if [[ $molecules -eq 2 ]]; then
    swaps=(1 -1)
fi

awk '/^MODEL /{
                $2+='$n_entries'
                header=$0
            }

    /^ATOM/{
                if ($5 == 0){
                    $5=1;
                    $2+='$atoms_per_molecule'
                    $6+='$residues_per_molecule'
                    $12+=1
                }
                else if ($5 = 1){
                    $5=0;
                    $2-='$atoms_per_molecule'
                    $6-='$residues_per_molecule'
                    $12-=1
                }
                lines[$2]=sprintf("ATOM  %5i  %-3s %-3s %-2i %-6i%8.3f%8.3f%8.3f%6.2f%6.2f      %-4i  ", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12);
            }

    /^ENDMDL/{
                print header
                for (x=1; x<='$atoms'; x++)
                    print lines[x]
                print
            }
            ' $PDB  >$PDB.swapped
