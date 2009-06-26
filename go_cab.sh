#!/bin/bash
for t in `seq 1 150`; do 
/usr/bin/qsub -q medium@d0cabsrv1 -l nodes=1 -k oe -m ae aa_bs.sh 
/usr/bin/qsub -q medium@d0cabsrv2 -l nodes=1 -k oe -m ae aa_bs.sh 
done
