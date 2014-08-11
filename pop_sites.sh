#!/bin/bash
export PATH=/opt/python/bin:$PATH
export PYTHONPATH=/home/jpg/opt/python/lib/python2.7/site-packages:$PYTHONPATH
python ~/PopSites/pop_sites.py \
--buff 500 \
--fdr 0.05 \
--inpath ~/Data/test.txt \
--pop CEU \
--rep 100
exit
