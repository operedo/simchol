# simchol

Large-scale stochastic simulation of 3D random functions using LU decomposition method with MPI.

Compile/run "normal" binary
-----------------------

1) make clean; make
2) mnsubmit run.sh

Compile/run "static trace" binary
-----------------------

1) make -f makefile_static
2) edit extrae.xml to set the destination folder
3) edit run_trace_static.sh to set the destination/auxiliary folders
4) mnsubmit run_trace_static.sh



