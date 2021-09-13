function ldm {
  echo "Enabling LDM environment ..."
  echo "Initialize ..."
  unset GENIEDIR
  unset GENIEBIN
  unset GENIELIB
  unset G4DIR

  echo "Setting up GSL2.5 ..."
  export GSLDIR=/cvmfs/larsoft.opensciencegrid.org/products/gsl/v2_5/Linux64bit+3.10-2.17-prof
  echo "GSLDIR=$GSLDIR"
  export GSLBIN=$GSLDIR/bin
  echo "GSLBIN=$GSLBIN"
  export GSLLIB=$GSLDIR/lib
  echo "GSLLIB=$GSLLIB"

  echo "Setting up GENIE ..."
  export GENIEDIR=/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_00_06i/Linux64bit+3.10-2.17-c7-prof
  echo "GENIEDIR=$GENIEDIR"
  export GENIEBIN=$GENIEDIR/bin
  echo "GENIEBIN=$GENIEBIN"
  export GENIELIB=$GENIEDIR/lib
  echo "GENIELIB=$GENIELIB"
  echo "Setting up GENIE ... FINISHED"

  echo "Setting up GEANT4 ..."
  export G4DIR=/dune/app/users/wyjang/geant4/4.10.07/x86_64-sl7-gnu
  echo "G4DIR=$G4DIR"
  source $G4DIR/bin/geant4.sh
  echo "Setting up GEANT4 ... FINISHED"
}
