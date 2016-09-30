cd src



########################################################################
# This compiles dfftpack
#
cd dfftpack
make
cd ../
#
########################################################################



########################################################################
# This compiles and locally installs Voro++
#
cd voro++/voro++-0.4.6
make
make install
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
src="DoSPT.f90 misc.f90 read_trajectory.f90 fluidicity.f90 volume.f90"
#
# Add other files
#
src="$src voronoi.o lowess.f"
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
gfortran $flags $src $libs -o $out
#
########################################################################
