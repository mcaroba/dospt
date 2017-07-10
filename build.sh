cd src



########################################################################
# This compiles dfftpack
#
cd dfftpack
#make
cd ../
#
########################################################################



########################################################################
# This compiles and locally installs Voro++
#
cd voro++/voro++-0.4.6
#make
#make install
cd ../../
#
########################################################################



########################################################################
# Compile the Voro++ functions used as external C code that can be
# handled by Fortran
#
g++ -c voronoi.cc -Ivoro++/voro++-0.4.6/src/
#
########################################################################



########################################################################
# Define directories for the different libraries used
#
# Voro++ library
#
vorolibs="-Lvoro++/install/lib/ -Lvoro++/voro++-0.4.6/src/ -lvoro++"
#
# dfftpack library
#
dfftlibs="-Ldfftpack/ -ldfftpack"
#
# Generic stuff
#
otherlibs="-lstdc++"
#
# All libraries
#
libs="$otherlibs $dfftlibs $vorolibs"
########################################################################



########################################################################
# Flags
#
flags="-fopenmp -std=legacy"
#
########################################################################



########################################################################
# Source code files
#
# Original DoSPT files written (mostly) by Miguel Caro
#
main="DoSPT.f90"
subroutines="misc.f90 volume.f90 sort_supergroups.f90"
src="$main $subroutines"
modules="constants.f90 read_input.f90 good_bye.f90 read_trajectory.f90 \
         rebuild_topology.f90 calc_dos.f90 misc2.f90 partition.f90 \
         fluidicity.f90 thermodynamic.f90"
#debug="-g -fcheck=all -Wall"
#debug="-g -fcheck=all"
#
# Add other files
#
other="voronoi.o lowess.f dfftpack/*.o"
src="$src $other"
#
########################################################################



########################################################################
# Compile modules
#
gfortran $debug $flags -c $modules $subroutines
modules=`echo $modules | sed 's/f90/o/g'`
#
########################################################################



########################################################################
# DoSPT output binary file
#
outdir="../bin/"
mkdir -p $outdir
#
out="$outdir/DoSPT"
#
########################################################################



########################################################################
# Compile DoSPT
#
gfortran $debug $flags $src $libs $modules -o $out
#
########################################################################
