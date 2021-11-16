# Standard set up script for DUNE BSM CAFAna
# M Wallbank, University of Cincinnati <wallbank@fnal.gov>, November 2021

# As far as I know, there is no 'standard' CAFAna build available on the DUNE VMs;
# this is therefore specifically to set up my build.
# If you check out your own version of the repo, change the first line!

BASE=/dune/app/users/wallbank
DIRECTORY=${BASE}/lblpwgtools

source /grid/fermiapp/products/dune/setup_dune.sh
source ${DIRECTORY}/CAFAna/build/Linux/CAFAnaEnv.sh
