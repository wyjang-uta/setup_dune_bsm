echo "Enabling LDM environment ..."
echo "Initialize ..."
#unset GSLDIR
#unset GSLBIN
#unset GSLLIB
#unset GLOBESDIR
#unset GLOBESBIN
#unset GLOBESLIB
#unset GENIEDIR
#unset GENIEBIN
#unset GENIELIB
#unset G4DIR
#unset G4BIN
#unset G4LIB

echo "Setting up GSL2.5 ..."
setup gsl v2_5 -q prof

#export GSLDIR=/dune/app/users/wyjang/dune_bsm/gsl/2.5
#echo "GSLDIR=$GSLDIR"
#export GSLBIN=$GSLDIR/bin
#export PATH=$PATH:$GSLBIN
#echo "GSLBIN=$GSLBIN"
#export GSLLIB=$GSLDIR/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSLLIB
#echo "GSLLIB=$GSLLIB"
#export C_INCLUDE_PATH=$GSLDIR/include

echo "Setting up GSL2.5 ... DONE"

#echo "Setting up GLoBES 3.0.11 ..."
#
#export GLOBESDIR=/dune/app/users/wyjang/dune_bsm/globes/3.0.11
#echo "GLOBESDIR=$GLOBESDIR"
#export GLOBESBIN=$GLOBESDIR/bin
#echo "GLOBESBIN=$GLOBESBIN"
#export PATH=$PATH:$GLOBESBIN
#export GLOBESLIB=$GLOBESDIR/lib
#echo "GLOBESLIB=$GLOBESLIB"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GLOBESLIB
#
#echo "Setting up GLoBES 3.0.11 ... DONE"

echo "Setting up GEANT4.10.6 ..."
setup geant4 v4_10_6_p01d -q e20:prof

#export G4DIR=/dune/app/users/wyjang/geant4/4.10.07/x86_64-sl7-gnu
#echo "G4DIR=$G4DIR"
#export G4BIN=$G4DIR/bin
#export PATH=$PATH:$G4BIN
#source $G4DIR/bin/geant4.sh
#export G4LIB=$G4DIR/lib64
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB

echo "Setting up GEANT4 ... DONE"

echo "Setting up GENIE v3.0.6 ..."

setup genie v3_00_06p -q c7:prof
#export GENIEDIR=/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06i/Linux64bit+3.10-2.17-c7-prof
#echo "GENIEDIR=$GENIEDIR"
#export GENIEBIN=$GENIEDIR/bin
#export PATH=$PATH:$GENIEBIN
#echo "GENIEBIN=$GENIEBIN"
#export GENIELIB=$GENIEDIR/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENIELIB
#echo "GENIELIB=$GENIELIB"

echo "Setting up GENIE v3.0.6 ... DONE"

