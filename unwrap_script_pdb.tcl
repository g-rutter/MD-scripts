#example usage:
#vmd out.ID11.pdb -dispdev text -e /home/theory/phrlaq/scripts/MD/unwrap_script_pdb.tcl -args 11
#produces T11-unwrapped.dcd
#set pbc correctly!

set dir ~/.dots/vmd/pbctools
source $dir/pkgIndex.tcl
package require pbctools

set temperature [lindex $argv 0]
set unwrapped_filename T$temperature-unwrapped.dcd

pbc set {190.0 190.0 190.0} -all
pbc join chain -all -bondlist
animate write dcd $unwrapped_filename
exit
