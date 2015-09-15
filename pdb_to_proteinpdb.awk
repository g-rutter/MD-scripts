#!/usr/bin/awk -f

#Script to make my PDB files be recognised as the correct number of contiguous proteins
#The steps are:-
    #- [A-Z]SC -> O  
    #- Give gly an O atom
    #- Add one-letter element in col 78

!/^ATOM/{
    print
}

/^ATOM/{
    #Set ASC atom names to be fake oxygens
    if ($3 ~ "[A-Z]SC") { name="O" }
    else { name=$3 }
    #Get element field for col78
    element=substr(name, 0, 1)

    #Print out fixed line
    printf("ATOM  %5i  %-3s %-3s %-2i %-6i%8.3f%8.3f%8.3f%6.2f%6.2f      %-4i %s\n", $2, name, $4, $5, $6, $7, $8, $9, $10, $11, $12, element);

    #Give gly an oxygen too
    if ($4 == "G" && $3 == "C"){
        printf("ATOM  %5i  %-3s %-3s %-2i %-6i%8.3f%8.3f%8.3f%6.2f%6.2f      %-4i %s\n", $2, "O", $4, $5, $6, $7, $8, $9, $10, $11, $12, "O");
    }
}
