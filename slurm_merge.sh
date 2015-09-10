#!/bin/bash
# @ initialdir = .
# @ output = seq_merge.out
# @ error =  seq_merge.err
# @ total_tasks = 1
# @ cpus_per_task = 1
# @ tasks_per_node = 4
# @ wall_clock_limit = 01:00:00

export EXTRAE_HOME=/gpfs/apps/NVIDIA/CEPBATOOLS/extrae/2.2.1/64

${EXTRAE_HOME}/bin/mpi2prv -syn -f /gpfs/scratch/bsc21/bsc21021/traces/TRACE.mpits -e ./simchol -o case.prv


