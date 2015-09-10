#!/bin/bash
# @ initialdir = .
# @ output = job_trace_%j.out
# @ error =  job_trace_%j.err
# @ total_tasks = 12
# @ cpus_per_task = 1
# @ tasks_per_node = 12
# @ tracing =  1
# @ wall_clock_limit = 00:15:00


SIMS=1

#MATSIZE=31127
#BLOCKSIZE=512

#MATSIZE=117527
#BLOCKSIZE=512

MATSIZE=9527
BLOCKSIZE=32



export OMP_NUM_THREADS=1
export OMP_SCHEDULE="dynamic"

mpirun -mca mpi_yield_when_idle 0 ./trace_static.sh ./simchol_trace_static muestras.dat >> salida_trace_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 2>&1  

	

mv /gpfs/scratch/bsc21/bsc21021/traces/case.prv /gpfs/scratch/bsc21/bsc21021/traces/caseNP${SLURM_NPROCS}NT${OMP_NUM_THREADS}.prv
mv /gpfs/scratch/bsc21/bsc21021/traces/case.pcf /gpfs/scratch/bsc21/bsc21021/traces/caseNP${SLURM_NPROCS}NT${OMP_NUM_THREADS}.pcf
mv /gpfs/scratch/bsc21/bsc21021/traces/case.row /gpfs/scratch/bsc21/bsc21021/traces/caseNP${SLURM_NPROCS}NT${OMP_NUM_THREADS}.row

rm /gpfs/scratch/bsc21/bsc21021/traces/TRACE* /gpfs/scratch/bsc21/bsc21021/traces/set-0/*

