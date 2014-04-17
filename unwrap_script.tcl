#example usage:
#for i in {1..15}; do vmd T${i}.dcd -psf ../*.psf  -e unwrap_script.tcl; mv traj_unwrapped.dcd T${i}.dcd; done
set dir ~/.dots/vmd/pbctools
source $dir/pkgIndex.tcl
package require pbctools

pbc join chain -all -bondlist
animate write dcd {traj_unwrapped.dcd}
exit
