# RippleTank


in its current state, just run make and it should build.
to run on the cluster, make sure to set the $RIPPLE_PROJ env var to the location of your RippleTank dir 

running as openmp, it should not put output to result.* but instead to source.dat


To plot a simuation(not using mpi), do the following
gnuplot -c animate.gp output.gif source.dat 10


once mpi and halo exchange is working we can figure out a new way to plot as well as put data in one file and in order


Make sure to create your own branch for changes so we dont loose work!!
