#!/bin/bash
# Identify population-specific SNP sites from an input set of sites.
python3 ~/aaylward/PopSites/pop_sites.py \
--buff 500 \
--fdr 0.05 \
--inpath ~/aaylward/Data/test.txt \
--pop CEU \
--rep 100
exit
