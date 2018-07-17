#!/usr/bin/gnuplot --persist

# csv must be , and space seperated
#gnuplot -c animate.gp 'out_put.gif' 'datafile.dat'
set term gif animate delay 1
set output ARG1
stats ARG2 nooutput
do for [i=1:int(STATS_blocks)] {
    splot ARG2 matrix index (i-1)
}


