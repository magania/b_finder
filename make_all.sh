#!/bin/bash

echo "AA=/home/magania/Bs/aatrack/P17/AA" > Makefile
echo "DFLAGS=MC" >> Makefile
cat Makefile.template >> Makefile
make -j4
make clean
echo "AA=/home/magania/Bs/aatrack/P17/AA" > Makefile
echo "DFLAGS=P17" >> Makefile
cat Makefile.template >> Makefile
make -j4
make clean
echo "AA=/home/magania/Bs/aatrack/P20/AA" > Makefile
echo "DFLAGS=P20" >> Makefile
cat Makefile.template >> Makefile
make -j4
make clean
echo "AA=/home/magania/Bs/aatrack/P21/AA" > Makefile
echo "DFLAGS=P21" >> Makefile
cat Makefile.template >> Makefile
make -j4
make clean
ln -s bs_finder_P17 bs_finder_P14

