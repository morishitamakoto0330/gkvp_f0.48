#!/bin/csh

###  NOTE  ###
###    Process number x Thread number =< 32 Cores x Node numbers
###    Available resource groups for Plasma Simulator:
###      small        #     1 - 12   nodes, 15 min., 32GB/nodes  2run/1submit
###      medium       #    25 - 1152 nodes,  5 hour, 32GB/nodes  4run/8submit
###      large1h      #  1921 - 2160 nodes,  1 hour, 32GB/nodes  1run/2submit
###      large        #  1921 - 2160 nodes, 10 hour, 32GB/nodes  1run/2submit
###    Available resource groups for LHD analysis server:
###      X24          #     1 - 24   nodes, 24 hour, 32GB/nodes 16run/32submit
###
###   Mon/Tue/Thu/Fri # 9:00  - 21:00 # medium = small > large1h > laege
###                   # 21:00 - 22:00 # large1h > large > medium = small
###   Wed/Sat/Sun     # 0:00  - 24:00 # large > large1h > medium = small
###              
###    To submit a interactive job, "pjsub --interact sub.q"
###                    a batch job, "pjsub sub.q"
###            To check job status, "pjstat" for step jobs "pjstat -E"
###                  To delete job, "pjdel JOBID"
###     To show budget information, "pstime"
##############


#PJM -L "rscunit=fx"
#PJM -L "rscgrp=X24"
#PJM -L "node=8"
#PJM -L "elapse=00:10:00"
#PJM -j
#PJM --mpi "proc=32"
#### --mpi "rank-map-hostfile=myrankmap"
#PJM -g 16391

## Note that Max. core num. per 1 node on PS is 32. 

setenv PARALLEL 8          # Thread number for automatic parallelization
setenv OMP_NUM_THREADS 8   # Thread number for Open MP


set DIR=%%DIR%%
set LDM=gkvp_mpifft.exe
set NL=gkvp_f0.48_namelist.%%%

### Run
cd ${DIR}
setenv fu05 ${DIR}/${NL}

module load fftw-fx/3.3.4
#module load fftw-fx/3.3.4-4simd

date
#rm -rf ${DIR}/Fprofd*
#mkdir ${DIR}/Fprofd_Stati
#fapp -C -d ${DIR}/Fprofd_Stati -Impi,hwm -L1 -Hevent=Statistics mpiexec ${DIR}/${LDM}
#fipp -C -d ${DIR}/Fprofd_Stati -Ihwm -Srange mpiexec ${DIR}/${LDM}
mpiexec ${DIR}/${LDM}
date 
touch complete

