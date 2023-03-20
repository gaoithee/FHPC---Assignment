#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=dp-mkl-cl
#SBATCH --nodes=1
#SBATCH --nodelist=epyc[007]
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load mkl

export OMP_PLACES=cores
export OMP_PROC_BIND=close

for numcores in {1..64}
        do export OMP_NUM_THREADS=$numcores
        echo $numcores >> weak-dp-mkl-cl.dat
        for i in {1..10}
                do ./gemm_mkl.x 12000 12000 12000 | grep GFLOPS >> weak-dp-mkl-cl.dat
        done
done
