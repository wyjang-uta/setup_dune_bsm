#!/bin/bash
gmkspl_dm -m 1 2 5 10 \
	-f ../geom/nd_hall_only_lar.gdml \
	-n 500 \
	-e 120 \
	-g 1 \
	-z 0.1 \
	--event-generator-list DM \
	-o ldm_xsec_mass_1_2_5_10_n500_e120_only_lar.xml
