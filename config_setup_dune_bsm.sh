#!/bin/bash
if [[ ":$PATH:" == *":$PWD:"* ]]; then
	echo "setup_dune_bsm is already configured."
else
echo '
# automatically added by setup_dune_bsm - begin
PATH='"$PWD"':$PATH
# automatically added by setup_dune_bsm - end' >> $HOME/.bash_profile
fi
