#!/bin/sh
#PJM -N "opt_mma_helmholtzfilter"
#PJM -L rscgrp=cx-single
#PJM -L node=1
#PJM -L elapse=1:00:00
#PJM -j
#PJM -o log.lst

#module load gcc/8.4.0 cuda/11.2.1
module load intel
./opt_mma_helmholtzfilter.out