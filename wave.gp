#!/usr/bin/gnuplot --persist

# csv must be , and space seperated
#gnuplot -c wave.gp 'out_put.png' 'datafile.dat'
set term png
set output ARG1
splot ARG2 matrix


