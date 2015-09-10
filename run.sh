#!/bin/bash
# @ job_name = my_job
# @ initialdir = .
# @ output = job_%j.out
# @ error =  job_%j.err
# @ wall_clock_limit = 5:00:00
# @ total_tasks = 30
# @ cpus_per_task = 12
# @ tasks_per_node = 1


SIMS=1

#MATSIZE=31127
#BLOCKSIZE=512

MATSIZE=117527
BLOCKSIZE=512

#MATSIZE=9527
#BLOCKSIZE=32

sl_get_machine_list -u > machinelist.txt
cat machinelist.txt | while read a ; do
ssh -n $a "$PWD/check.sh 5" &
done

#srun -N $SLURM_NNODES -n $SLURM_NNODES $PWD/check.sh 5 &
 
      export OMP_NUM_THREADS=12
      export OMP_SCHEDULE="dynamic"
      echo "---------------------------------------------" >> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 
      echo 'START: JOB SUBMITTED WITH SRUN AT: ' `date` >> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 
      echo "SLURM_NPROCS=$SLURM_NPROCS: ">> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 
      echo "OMP_NUM_THREADS=$OMP_NUM_THREADS: ">> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 
      #/usr/bin/time srun ./LUB_openmp_cova muestras.dat >> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS} 2>&1
      /usr/bin/time srun ./simchol muestras-nscore.out_deletedRows >> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 2>&1
      echo 'END: JOB SUBMITTED WITH SRUN AT: ' `date` >> salida_NP_${SLURM_NPROCS}_NT_${OMP_NUM_THREADS}_${SIMS}_${MATSIZE}_${BLOCKSIZE} 


cat machinelist.txt | while read a; do
ssh -n $a "killall check.sh" &
ssh -n $a "killall top" &
ssh -n $a "pkill -u bsc21021" &
done

