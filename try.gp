#!/usr/bin/gnuplot --persist

# csv must be , and space seperated
#gnuplot -c animate.gp 'out_put.gif' 'datafile.dat' 'min/max scale'
set term gif animate delay 1
set output ARG1
stats ARG2 nooutput
set zrange [-ARG3:ARG3]
set dgrid3d 60,60 gauss 1
unset  key
set view 60, 75, 1, 1
set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set style data lines
DEBUG_TERM_HTIC = 119
DEBUG_TERM_VTIC = 119

do for [i=1:int(STATS_blocks)] {
    splot ARG2 matrix index (i-1) u 1:2:3 w l
}


