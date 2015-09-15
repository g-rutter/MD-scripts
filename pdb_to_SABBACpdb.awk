#!/usr/bin/awk -f

#Script to make my PDB files be recognised by SABBAC
#The steps are:-
    #- Remove all  atoms except CA
    #- Write resname in 3 letter format
    #- (Maybe not necessary) Add one-letter element in col 78

/^ATOM.* CA /{
    #Get element field for col78
    name=$3
    element=substr(name, 0, 1)

    oneletter=$4
    if (oneletter ~/A/) threeletter="ALA"
    if (oneletter ~/R/) threeletter="ARG"
    if (oneletter ~/N/) threeletter="ASN"
    if (oneletter ~/D/) threeletter="ASP"
    if (oneletter ~/C/) threeletter="CYS"
    if (oneletter ~/Q/) threeletter="GLN"
    if (oneletter ~/E/) threeletter="GLU"
    if (oneletter ~/G/) threeletter="GLY"
    if (oneletter ~/H/) threeletter="HIS"
    if (oneletter ~/I/) threeletter="ILE"
    if (oneletter ~/L/) threeletter="LEU"
    if (oneletter ~/K/) threeletter="LYS"
    if (oneletter ~/M/) threeletter="MET"
    if (oneletter ~/F/) threeletter="PHE"
    if (oneletter ~/P/) threeletter="PRO"
    if (oneletter ~/S/) threeletter="SER"
    if (oneletter ~/T/) threeletter="THR"
    if (oneletter ~/W/) threeletter="TRP"
    if (oneletter ~/Y/) threeletter="TYR"
    if (oneletter ~/V/) threeletter="VAL"

    #Print out fixed line
    printf("ATOM  %5i  %-3s %-3s %-2i %-6i%8.3f%8.3f%8.3f%6.2f%6.2f      %-4i %s\n", $2, name, threeletter, $5, $6, $7, $8, $9, $10, $11, $12, element);
}

