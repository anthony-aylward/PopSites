#!/bin/bash
# Identify population-specific SNP sites from an input set of sites.
python3 ~/PopSites/pop_sites.py \
--buff 1000 \
--fdr 0.05 \
--gtf ~/GTF/hg19_genes.gtf \
--inpath ~/Data/HPG_genes_with_sites_after_constraint9.txt \
--pop CEU \
--rep 100
exit
