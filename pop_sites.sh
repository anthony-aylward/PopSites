#!/bin/bash
# Identify population-specific SNP sites from an input set of sites.
opt/python/bin/python3 ~/PopSites/pop_sites.py \
--buff 500 \
--fdr 0.05 \
--inpath ~/Data/test.txt \
--pop CEU \
--rep 100
exit
