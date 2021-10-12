# HNL analysis requirement summary
# Geant4 version 4.10.6
# GENIE version 3.0.6
# GSL version 2.5

echo "Enabling Heavy Neutral Lepton environment ..."

unset GSLDIR
unset GSLBIN
unset GSLLIB
unset GENIEDIR
unset GENIEBIN
unset GENIELIB
unset G4DIR
unset G4BIN
unset G4LIB

echo "Setting GSL2.5 ..."
setup gsl v2_5 -q prof
echo "Setting GSL2.5 ... DONE"

echo "Setting GENIE ..."
setup genie v3_00_06p -q c7:prof
echo "Setting up GENIE ... DONE"

echo "Setting up GEANT4 ..."
setup geant4 v4_10_6_p01d -q e20:prof
echo "Setting up GEANT4 ... DONE"

