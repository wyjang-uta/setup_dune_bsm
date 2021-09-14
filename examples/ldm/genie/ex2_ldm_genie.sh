
#!/bin/bash
# Script to run DM interaction generation - wyjang Dec. 20. 2020.
# Scalar DM run. !Use GDM10_00a_00_000 for Ferminoic DM run

# -m : mass of dark matter
# -z : the ratio of the mediator mass to the dark matter mass. Default: 0.5
# -g : the Z' gauge coupling. Default: taken from $GENIE/config/CommonParam.xml -> 0.1

gevgen_dm -t 1000180400 \
  -m 2 \
  -z 0.1 \
  -g 0.1 \
  -n 100000 \
  -e 10 \
  --tune GDM18_00b_00_000 \
  --run 1000 \
  --seed 172188 \
  --event-generator-list DME \
  -o test_dme_m2.root \
  --cross-sections ldm_xsec_n500_e120_only_lar.xml
