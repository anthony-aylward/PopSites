#!/bin/bash
# Identify population-specific SNP sites from an input set of sites.
python3 ~/PopSites/pop_sites.py \
--buff 500 \
--fdr 0.05 \
--inpath ~/PopSites/Data/test.txt \
--pop EUR \
--rep 100
exit
