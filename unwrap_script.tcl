#example usage:
#for i in {1..15}; do vmd T${i}.dcd -psf ../*.psf  -e unwrap_script.tcl; mv traj_unwrapped.dcd T${i}.dcd; done
set dir ~/.dots/vmd/pbctools
source $dir/pkgIndex.tcl
package require pbctools

set temperature [lindex $argv 0]
set unwrapped_filename T$temperature-unwrapped.dcd

pbc join chain -all -bondlist
animate write dcd $unwrapped_filename
exit
