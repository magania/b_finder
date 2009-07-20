#!/bin/bash
for t in `seq 1 250`; do 
/usr/bin/qsub -q medium@d0cabsrv1 -l nodes=1 -k oe -m ae a2_jpsi.sh 
/usr/bin/qsub -q medium@d0cabsrv2 -l nodes=1 -k oe -m ae a2_jpsi.sh 
done
