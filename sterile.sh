# LED/Sterile neutrino analysis requirement summary
# (1) GSL v2.5
# (2) GLoBES v3.2.18

echo "Enabling Sterile Neutrino environment ..."
echo "Initialize ..."
unset GSLDIR
unset GSLBIN
unset GSLLIB
unset GLOBESDIR
unset GLOBESBIN
unset GLOBESLIB

echo "Setting up GSL2.5 ..."
setup gsl v2_5 -q prof
echo "Setting up GSL2.5 ... DONE"

echo "Setting up GLoBES 3.0.11 ..."

export GLOBESDIR=/dune/app/users/wyjang/dune_bsm/globes/3.0.11
echo "GLOBESDIR=$GLOBESDIR"
export GLOBESBIN=$GLOBESDIR/bin
echo "GLOBESBIN=$GLOBESBIN"
export PATH=$PATH:$GLOBESBIN
export GLOBESLIB=$GLOBESDIR/lib
echo "GLOBESLIB=$GLOBESLIB"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GLOBESLIB

echo "Setting up GLoBES 3.0.11 ... DONE"
