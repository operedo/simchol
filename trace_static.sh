#!/bin/bash
#export MPTRACE_CONFIG_FILE=./mpitrace.xml
#export MPITRACE_ON=1
#export LD_LIBRARY_PATH=/apps/PAPI/4.4.0/lib:/gpfs/apps/NVIDIA/CEPBATOOLS/extrae/latest/64/lib
#export LD_PRELOAD=libmpitrace.so
#$@


export EXTRAE_HOME=/gpfs/apps/NVIDIA/CEPBATOOLS/extrae/2.2.1/64
export EXTRAE_CONFIG_FILE=./extrae.xml
source ${EXTRAE_HOME}/etc/extrae.sh

## Run the desired program
$*
