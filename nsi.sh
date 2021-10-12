# NSI software requirement summary
# (1) GLoBES 3.2.18

echo "Enabling Non-standard Interaction environment ..."
echo "Initialize ..."
unset GLOBESDIR
unset GLOBESBIN
unset GLOBESLIB

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
